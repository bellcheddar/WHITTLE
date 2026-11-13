"""Antigen constructs for Tier 2. Spec 5.2 and 6.4 item 1.

Fetches each target's canonical UniProt sequence once, slices the construct to
the boundaries in `config/targets.yaml`, and caches it under
`data/antigens/`. The boundaries are authoritative and are not re-derived here:
they are UniProt numbering for the full-length precursor, accessed 2026-11-13,
and spec 5.2 says not to round them or substitute a PDB construct without a
DECISIONS.md line.

The md5 of each construct goes into the `targets` table, so a later run that
silently fetched a different isoform would show up as a changed hash rather
than as quietly different co-folds.

    pixi run python -m pipeline.antigens
"""
from __future__ import annotations

import hashlib
import sys

from pipeline.common import (ANTIGENS, RateLimiter, connect, http_get, log,
                             manifest_write, targets_config, write_json)

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"


def fetch_sequence(accession: str, limiter: RateLimiter) -> str:
    cached = ANTIGENS / f"{accession}.fasta"
    if cached.exists():
        text = cached.read_text()
    else:
        response = http_get(UNIPROT_FASTA.format(acc=accession), limiter=limiter,
                            timeout=90)
        if response.status_code != 200 or not response.text.startswith(">"):
            raise RuntimeError(f"UniProt returned {response.status_code} for {accession}")
        text = response.text
        cached.write_text(text)
    return "".join(line.strip() for line in text.splitlines()[1:])


def main() -> int:
    cfg = targets_config()["targets"]
    limiter = RateLimiter(3)
    conn = connect()
    summary = {}

    for name, spec in cfg.items():
        accession = spec["uniprot"]
        try:
            full = fetch_sequence(accession, limiter)
        except Exception as exc:
            log(f"antigens: {name} ({accession}) FAILED: {exc}")
            manifest_write("antigens", {"id": name, "status": f"failed:{exc}"[:120]})
            continue
        start, end = spec["range"]
        if len(full) < end:
            log(f"antigens: {name} ({accession}) UniProt entry is {len(full)} aa but "
                f"the construct runs to {end}; the isoform may have changed")
            manifest_write("antigens", {"id": name, "status": "failed:short entry"})
            continue
        # UniProt numbering is 1-based and inclusive at both ends.
        construct = full[start - 1:end]
        if len(construct) != spec["length"]:
            log(f"antigens: {name} construct is {len(construct)} aa, "
                f"config/targets.yaml says {spec['length']}. Using the sliced "
                f"sequence and flagging the disagreement")
        digest = hashlib.md5(construct.encode()).hexdigest()
        stoich = int(spec["stoichiometry"])
        path = ANTIGENS / f"{name}.fasta"
        header = (f">{name}|{accession}|{spec['construct']}|{start}-{end}"
                  f"|x{stoich}|md5:{digest}")
        path.write_text(header + "\n" + "\n".join(
            construct[i:i + 60] for i in range(0, len(construct), 60)) + "\n")

        conn.execute("UPDATE targets SET antigen_md5 = ? WHERE target = ?",
                     (digest, name))
        summary[name] = {
            "uniprot": accession, "full_length": len(full),
            "range": [start, end], "chain_length": len(construct),
            "stoichiometry": stoich, "antigen_length": len(construct) * stoich,
            "complex_length_estimate": len(construct) * stoich + 125,
            "md5": digest, "path": str(path.relative_to(ANTIGENS.parent.parent)),
            "config_length": spec["length"],
            "length_matches_config": len(construct) == spec["length"],
        }
        manifest_write("antigens", {"id": name, "status": "ok", **summary[name]})
        log(f"antigens: {name} {accession} {start}-{end} -> {len(construct)} aa "
            f"x{stoich} = {len(construct) * stoich} aa antigen, md5 {digest[:8]}"
            + ("" if len(construct) == spec["length"]
               else f" (config says {spec['length']} aa)"))

    conn.commit()
    conn.close()
    write_json(ANTIGENS / "antigens.json", summary)
    order = sorted(summary, key=lambda k: summary[k]["complex_length_estimate"])
    log(f"antigens: {len(summary)} constructs cached. Cheapest first by complex "
        f"size: {', '.join(order)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
