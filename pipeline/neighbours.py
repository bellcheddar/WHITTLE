"""Nearest experimental neighbour in SAbDab-nano. Spec 6.2 items 7 and 8.

Three things come out of this stage:

1. **Per-region identity** to the closest experimentally determined nanobody.
   Reported per region, not whole-chain, because framework identity is above
   90 percent for essentially every VHH and tells you nothing. The CDR identity
   is the number that carries information.
2. **TM-score** between this candidate's predicted model and that experimental
   structure, via TM-align (`tmtools`).
3. **The experimental-structure badge**: a near-identical sequence has a PDB
   entry. **Positive signal only** (spec 6.2 item 8). Its absence says nothing
   about expressibility and never reduces a score; nothing in
   `pipeline/score.py` reads this column.

Spec 5.3: SAbDab-nano is **not** mirrored into the repository or onto the
droplet. The summary is downloaded into the gitignored `data/sabdab/`, only
derived columns are stored (nearest PDB id, TM-score, identities, same-target
flag), and the snapshot date goes into the provenance table.

    pixi run python -m pipeline.neighbours
    pixi run python -m pipeline.neighbours --target EGFR --no-structures
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import io
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from pipeline.common import (DATA, STRUCTURES, TARGETS, connect, http_get, log,
                             manifest_done, manifest_write, write_json, RateLimiter)

SABDAB = DATA / "sabdab"
SUMMARY = SABDAB / "sabdab_nano_summary.tsv"
PDB_CACHE = SABDAB / "pdb"
SUMMARY_URL = ("https://sabdab.opig.stats.ox.ac.uk/api/download/"
               "all-single-domain-summary")

# Which SAbDab antigen names count as a target-matched reference. Only EGFR and
# TNFa have known nanobodies (spec 5.1 item 4); the other three are left empty
# rather than matched loosely.
TARGET_ANTIGEN_PATTERNS = {
    "EGFR": ("epidermal growth factor receptor", "egfr", "erbb-1", "erbb1"),
    "TNFa": ("tumor necrosis factor", "tumour necrosis factor", "tnf-alpha",
             "tnf alpha", "tnfa"),
}

BADGE_IDENTITY = 0.95
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


# ------------------------------------------------------------------ reference

def download_summary(force: bool = False) -> Path:
    SABDAB.mkdir(parents=True, exist_ok=True)
    if SUMMARY.exists() and not force:
        return SUMMARY
    log("neighbours: downloading the SAbDab single-domain summary")
    response = http_get(SUMMARY_URL, timeout=300)
    SUMMARY.write_bytes(response.content)
    return SUMMARY


def load_reference() -> list[dict]:
    """One row per unique VH sequence, with every PDB entry that carries it."""
    path = download_summary()
    by_sequence: dict[str, dict] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            seq = (row.get("VH") or "").strip().upper()
            if len(seq) < 90 or len(seq) > 200:
                continue
            if not seq.isalpha():
                continue
            entry = by_sequence.setdefault(seq, {
                "sequence": seq,
                "cdr1": (row.get("CDR-H1") or "").strip(),
                "cdr2": (row.get("CDR-H2") or "").strip(),
                "cdr3": (row.get("CDR-H3") or "").strip(),
                "pdbs": [], "chains": [], "antigens": set(), "dates": [],
            })
            pdb = (row.get("PDB") or "").strip()
            # SAbDab2 ids look like pdb_000010nd; the last four characters are
            # the classic PDB code.
            short = pdb.replace("pdb_0000", "").lower() if pdb.startswith("pdb_") else pdb.lower()
            if short and short not in entry["pdbs"]:
                entry["pdbs"].append(short)
                entry["chains"].append((row.get("Hchain") or "").strip())
            antigen = (row.get("antigen_name") or "").strip().lower()
            if antigen and antigen not in ("na", "none"):
                entry["antigens"].add(antigen)
            date = (row.get("date") or "").strip()
            if date:
                entry["dates"].append(date)
    out = []
    for entry in by_sequence.values():
        entry["antigens"] = sorted(entry["antigens"])
        entry["pdb"] = entry["pdbs"][0] if entry["pdbs"] else ""
        entry["chain"] = entry["chains"][0] if entry["chains"] else "H"
        out.append(entry)
    log(f"neighbours: SAbDab-nano reference set, {len(out)} unique VH sequences "
        f"from {sum(len(e['pdbs']) for e in out)} structure instances")
    return out


def antigen_matches(antigens: list[str], target: str) -> bool:
    patterns = TARGET_ANTIGEN_PATTERNS.get(target)
    if not patterns:
        return False
    joined = " | ".join(antigens)
    return any(p in joined for p in patterns)


# ------------------------------------------------------------------- identity

def identity(a: str, b: str) -> float:
    """Identity over the longer sequence, from a global alignment."""
    if not a or not b:
        return 0.0
    if a == b:
        return 1.0
    la, lb = len(a), len(b)
    prev = list(range(lb + 1))
    for i in range(1, la + 1):
        row = [i] + [0] * lb
        ai = a[i - 1]
        for j in range(1, lb + 1):
            row[j] = min(prev[j] + 1, row[j - 1] + 1,
                         prev[j - 1] + (ai != b[j - 1]))
        prev = row
    return 1.0 - prev[lb] / max(la, lb)


def shortlist(candidates: list[dict], reference: list[dict], top_n: int
              ) -> dict[str, list[int]]:
    """Top-N reference rows per candidate, by MMseqs2 where it is available.

    12 million exact pairwise alignments is not affordable, so MMseqs2 does the
    search and the exact region-resolved identity is computed only on what it
    returns. Without MMseqs2 the fallback is a cheap 4-mer Jaccard prefilter,
    which is coarse but keeps the stage running.
    """
    binary = shutil.which("mmseqs")
    if binary is not None:
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            (tmpdir / "q.fasta").write_text(
                "".join(f">{c['id']}\n{c['sequence']}\n" for c in candidates))
            (tmpdir / "t.fasta").write_text(
                "".join(f">r{i}\n{r['sequence']}\n" for i, r in enumerate(reference)))
            try:
                subprocess.run(
                    [binary, "easy-search", str(tmpdir / "q.fasta"),
                     str(tmpdir / "t.fasta"), str(tmpdir / "hits.m8"),
                     str(tmpdir / "work"), "--max-seqs", str(top_n * 4),
                     "-s", "6.0", "--threads", "8", "-v", "0",
                     "--format-output", "query,target,fident,alnlen,bits"],
                    capture_output=True, check=True, timeout=1800)
                hits: dict[str, list[tuple[float, int]]] = {}
                for line in (tmpdir / "hits.m8").read_text().splitlines():
                    parts = line.split("\t")
                    if len(parts) < 5:
                        continue
                    query, target, fident = parts[0], parts[1], float(parts[2])
                    hits.setdefault(query, []).append((fident, int(target[1:])))
                out = {}
                for cid, rows in hits.items():
                    rows.sort(reverse=True)
                    out[cid] = [idx for _f, idx in rows[:top_n]]
                if out:
                    log(f"neighbours: MMseqs2 shortlisted hits for {len(out)} "
                        f"of {len(candidates)} candidates")
                    return out
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired,
                    FileNotFoundError, ValueError):
                log("neighbours: MMseqs2 search failed; falling back to a k-mer "
                    "prefilter")

    def kmers(seq: str, k: int = 4) -> set[str]:
        return {seq[i:i + k] for i in range(len(seq) - k + 1)}

    ref_kmers = [kmers(r["sequence"]) for r in reference]
    out = {}
    for cand in candidates:
        qk = kmers(cand["sequence"])
        scored = sorted(
            ((len(qk & rk) / max(len(qk | rk), 1), i) for i, rk in enumerate(ref_kmers)),
            reverse=True)
        out[cand["id"]] = [i for _s, i in scored[:top_n]]
    return out


# ------------------------------------------------------------------ TM-score

def fetch_chain(pdb_id: str, chain: str, limiter: RateLimiter) -> list[tuple] | None:
    """CA coordinates and the one-letter sequence for one chain of one entry."""
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    cached = PDB_CACHE / f"{pdb_id}.cif.gz"
    if not cached.exists():
        try:
            response = http_get(f"https://files.rcsb.org/download/{pdb_id}.cif.gz",
                                limiter=limiter, timeout=90)
            if response.status_code != 200 or not response.content:
                return None
            cached.write_bytes(response.content)
        except Exception:
            return None
    try:
        import gemmi

        # gemmi reads a gzipped mmCIF straight off disk; this version has no
        # read_structure_from_string, and going via a string gains nothing.
        structure = gemmi.read_structure(str(cached))
        structure.setup_entities()
        model = structure[0]
        for ch in model:
            if ch.name != chain:
                continue
            coords, seq = [], []
            for residue in ch:
                atom = residue.find_atom("CA", "*")
                one = gemmi.find_tabulated_residue(residue.name)
                if atom is None or one is None or not one.is_amino_acid():
                    continue
                coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
                seq.append(THREE_TO_ONE.get(residue.name.upper(), "X"))
            if len(coords) >= 50:
                return [coords, "".join(seq)]
    except Exception:
        return None
    return None


def model_chain(path: Path) -> list[tuple] | None:
    try:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            coords, seq = [], []
            for line in fh:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    coords.append((float(line[30:38]), float(line[38:46]),
                                   float(line[46:54])))
                    seq.append(THREE_TO_ONE.get(line[17:20].strip().upper(), "X"))
        return [coords, "".join(seq)] if len(coords) >= 50 else None
    except Exception:
        return None


def tm_align_result(model: list, experimental: list):
    try:
        import numpy as np
        from tmtools import tm_align

        return tm_align(np.asarray(model[0], dtype=np.float64),
                        np.asarray(experimental[0], dtype=np.float64),
                        model[1], experimental[1])
    except Exception:
        return None


def tm_score(model: list, experimental: list) -> float | None:
    result = tm_align_result(model, experimental)
    if result is None:
        return None
    # Normalised by the model, which is the chain we are making a claim about.
    return round(float(result.tm_norm_chain1), 4)


def write_superposed(model_path: Path, out_path: Path, result) -> bool:
    """Write the model with every atom moved onto the experimental frame.

    DESIGN 3.5 asks the Neighbours panel to lead with the superposition, and
    spec 6.2 item 7 calls it the panel that justifies having a 3D viewer rather
    than decorating with one. Two structures loaded in their own coordinate
    frames sit hundreds of angstroms apart and are not a superposition, so the
    transform TM-align already computed is applied here, once, at build time:
    x' = u.x + t, verified to drop the centroid separation from 256 to 1.2
    angstroms on the first pair it was tried against.

    Doing it here rather than in the browser keeps the droplet free of
    request-time computation and means the viewer only has to load two files.
    """
    import numpy as np

    try:
        u = np.asarray(result.u, dtype=np.float64)
        shift = np.asarray(result.t, dtype=np.float64)
        opener = gzip.open if str(model_path).endswith(".gz") else open
        with opener(model_path, "rt") as fh:
            lines = fh.read().splitlines()
        out: list[str] = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])])
                moved = u @ xyz + shift
                out.append(line[:30]
                           + f"{moved[0]:8.3f}{moved[1]:8.3f}{moved[2]:8.3f}"
                           + line[54:])
            else:
                out.append(line)
        header = ("REMARK  SUPERPOSED ONTO THE NEAREST EXPERIMENTAL NANOBODY "
                  "BY TM-ALIGN")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(out_path, "wt", compresslevel=9) as fh:
            fh.write(header + "\n" + "\n".join(out) + "\n")
        return True
    except Exception:
        return False


# ----------------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", action="append", default=None)
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--no-structures", action="store_true",
                        help="skip the TM-score step, which fetches PDB entries")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    targets = args.target or TARGETS
    reference = load_reference()
    conn = connect()
    done = manifest_done("neighbours")

    candidates = []
    for target in targets:
        for row in conn.execute(
                "SELECT id, target, sequence, cdr1, cdr2, cdr3 FROM candidates "
                "WHERE target = ? ORDER BY id", (target,)):
            if row["id"] in done:
                continue
            candidates.append(dict(row))
            if args.limit and len(candidates) >= args.limit:
                break
        if args.limit and len(candidates) >= args.limit:
            break
    conn.close()
    if not candidates:
        log("neighbours: nothing to do")
        return 0

    log(f"neighbours: {len(candidates)} candidates against "
        f"{len(reference)} unique experimental nanobody sequences")
    shortlists = shortlist(candidates, reference, args.top_n)

    limiter = RateLimiter(6)
    structure_cache: dict[str, list | None] = {}
    started = time.time()
    with_badge = with_tm = same_target = with_superposed = 0

    for n, cand in enumerate(candidates, 1):
        indices = shortlists.get(cand["id"], [])
        best = None
        for idx in indices:
            ref = reference[idx]
            whole = identity(cand["sequence"], ref["sequence"])
            if best is None or whole > best["identity_whole"]:
                best = {
                    "ref": ref,
                    "identity_whole": round(whole, 4),
                    "identity_cdr": round(
                        (identity(cand["cdr1"] or "", ref["cdr1"])
                         + identity(cand["cdr2"] or "", ref["cdr2"])
                         + identity(cand["cdr3"] or "", ref["cdr3"])) / 3.0, 4),
                    "identity_cdr3": round(
                        identity(cand["cdr3"] or "", ref["cdr3"]), 4),
                }
        if best is None:
            manifest_write("neighbours", {"id": cand["id"], "target": cand["target"],
                                          "status": "failed:no hit"})
            continue
        ref = best["ref"]
        # Framework identity, from the chain with the three CDRs removed.
        cand_fw = cand["sequence"]
        ref_fw = ref["sequence"]
        for key in ("cdr1", "cdr2", "cdr3"):
            if cand.get(key):
                cand_fw = cand_fw.replace(cand[key], "", 1)
            if ref.get(key):
                ref_fw = ref_fw.replace(ref[key], "", 1)
        identity_fw = round(identity(cand_fw, ref_fw), 4)

        record = {
            "id": cand["id"], "target": cand["target"], "status": "ok",
            "nearest_pdb": ref["pdb"].upper(),
            "nearest_chain": ref["chain"],
            "identity_whole": best["identity_whole"],
            "identity_fw": identity_fw,
            "identity_cdr": best["identity_cdr"],
            "identity_cdr3": best["identity_cdr3"],
            "nearest_antigens": "; ".join(ref["antigens"][:3]),
            "same_target_reference": 1 if antigen_matches(ref["antigens"],
                                                          cand["target"]) else 0,
            # Positive signal only (spec 6.2 item 8). Nothing in the scoring
            # reads this column.
            "badge": 1 if best["identity_whole"] >= BADGE_IDENTITY else 0,
            "n_pdb_entries": len(ref["pdbs"]),
        }
        if record["badge"]:
            with_badge += 1
        if record["same_target_reference"]:
            same_target += 1

        if not args.no_structures and ref["pdb"]:
            key = f"{ref['pdb']}_{ref['chain']}"
            if key not in structure_cache:
                structure_cache[key] = fetch_chain(ref["pdb"], ref["chain"], limiter)
            experimental = structure_cache[key]
            model_path = STRUCTURES / cand["target"] / f"{cand['id']}.pdb.gz"
            if experimental and model_path.exists():
                aligned = tm_align_result(model_chain(model_path), experimental)
                if aligned is not None:
                    record["nearest_tm"] = round(float(aligned.tm_norm_chain1), 4)
                    record["nearest_rmsd"] = round(float(aligned.rmsd), 3)
                    superposed = (STRUCTURES / cand["target"]
                                  / f"{cand['id']}.sup.pdb.gz")
                    with_tm += 1
                    if write_superposed(model_path, superposed, aligned):
                        record["superposed"] = 1
                        with_superposed += 1

        manifest_write("neighbours", record)
        if n % 100 == 0 or n == len(candidates):
            rate = n / max(time.time() - started, 1e-9)
            log(f"neighbours: {n}/{len(candidates)}, {with_badge} with an "
                f"experimental badge, {with_tm} with a TM-score, "
                f"{with_superposed} superposed, {rate * 60:.0f} per minute")

    write_json(SABDAB / "snapshot.json", {
        "source": "SAbDab-nano single-domain summary",
        "url": SUMMARY_URL,
        "snapshot_date": dt.datetime.fromtimestamp(
            SUMMARY.stat().st_mtime).strftime("%Y-%m-%d"),
        "unique_sequences": len(reference),
        "badge_identity_threshold": BADGE_IDENTITY,
        "note": "Derived values only. The database is not mirrored into this "
                "repository or onto the droplet (spec 5.3).",
    })
    log(f"neighbours: complete. {with_badge} candidates carry an experimental "
        f"badge, {same_target} have a target-matched reference, "
        f"{with_tm} have a TM-score and {with_superposed} have a superposed "
        f"model written for the drawer")
    return 0


if __name__ == "__main__":
    sys.exit(main())
