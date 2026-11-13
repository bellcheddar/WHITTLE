"""Reference table and provenance. Spec 8.3.

`config/references.yaml` holds the DOIs and roles; this module resolves each
one against the Crossref API **at build time on the Studio** and bakes static
JSON. Crossref is never called per visitor.

Adapted from `debris/pipeline/verify_doi.py` (see REUSE.md). The title-matching
check is the part carried over: a DOI that resolves to a different paper is a
worse error than a DOI that does not resolve, because it looks right. What is
new here is that a reference may be citable by `repo` plus `commit` rather than
by DOI, which several of these tools need, and that a database row shows an
accession date rather than a DOI alone (spec 8.3).

A row that cannot be verified keeps its hand-written fallback title, is flagged
`doi_verified = 0`, and the About tab marks it unverified (spec 0.2 fallback).
Nothing is silently dropped and nothing unverified is presented as verified.

    pixi run refs
"""
from __future__ import annotations

import datetime as dt
import difflib
import re
import subprocess
import sys

from pipeline.common import (ANTIGENS, MODELS, RateLimiter, connect, git_sha,
                             http_get, load_yaml, log, manifest_write, run_id)

CROSSREF = "https://api.crossref.org/works/{doi}"
TITLE_MATCH_FLOOR = 0.60


def normalise(title: str) -> str:
    return re.sub(r"[^a-z0-9 ]", " ", (title or "").lower()).strip()


def resolve_doi(doi: str, fallback_title: str, limiter: RateLimiter) -> dict:
    """Crossref metadata for one DOI, with a title-match sanity check."""
    try:
        response = http_get(CROSSREF.format(doi=doi), limiter=limiter, timeout=45)
    except Exception as exc:
        return {"doi_verified": 0, "title": fallback_title,
                "verify_note": f"Crossref unreachable: {type(exc).__name__}"}
    if response.status_code != 200:
        return {"doi_verified": 0, "title": fallback_title,
                "verify_note": f"Crossref returned {response.status_code}"}
    try:
        message = response.json()["message"]
    except Exception:
        return {"doi_verified": 0, "title": fallback_title,
                "verify_note": "Crossref returned something unparseable"}

    title = (message.get("title") or [fallback_title])[0]
    authors = message.get("author") or []
    if authors:
        first = authors[0]
        name = f"{first.get('family', '')} {first.get('given', '')[:1]}".strip()
        author_text = name + (" et al." if len(authors) > 1 else "")
    else:
        author_text = ""
    issued = (message.get("issued") or {}).get("date-parts") or [[None]]
    year = issued[0][0] if issued and issued[0] else None
    container = (message.get("container-title") or [""])[0]

    # A DOI that resolves to a different paper looks right and is worse than one
    # that does not resolve at all.
    ratio = difflib.SequenceMatcher(
        None, normalise(title), normalise(fallback_title)).ratio()
    verified = 1 if (not fallback_title or ratio >= TITLE_MATCH_FLOOR) else 0
    note = ("" if verified else
            f"Crossref title does not match the expected one "
            f"(similarity {ratio:.2f}): Crossref says \"{title[:90]}\"")
    return {
        "doi_verified": verified,
        "title": title if verified else fallback_title,
        "crossref_title": title,
        "authors": author_text,
        "year": year,
        "container": container,
        "url": message.get("URL") or f"https://doi.org/{doi}",
        "title_similarity": round(ratio, 3),
        "verify_note": note,
    }


def resolve_commit(repo: str) -> str:
    """The head commit of a GitHub repository, for a reference with no DOI."""
    match = re.search(r"github\.com/([^/]+)/([^/.]+)", repo or "")
    if not match:
        return ""
    owner, name = match.groups()
    try:
        response = http_get(
            f"https://api.github.com/repos/{owner}/{name}/commits?per_page=1",
            timeout=30)
        if response.status_code == 200:
            payload = response.json()
            if payload:
                return str(payload[0].get("sha", ""))[:12]
    except Exception:
        pass
    return ""


# ---------------------------------------------------------------- provenance

def tool_versions() -> dict[str, str]:
    """Version of every tool that actually ran, read from the environment."""
    import importlib.metadata as md

    wanted = {
        "anarci": "ANARCI", "anarcii": "anarcii", "immunebuilder": "ImmuneBuilder",
        "tnp": "TNP", "ablang2": "AbLang2", "promb": "promb (OASis peptides)",
        "fair-esm": "ESM-2", "scikit-learn": "scikit-learn", "torch": "PyTorch",
        "openmm": "OpenMM", "freesasa": "FreeSASA", "biopython": "Biopython",
        "gemmi": "gemmi", "tmtools": "TM-align (tmtools)",
        "prodigy-prot": "PRODIGY", "numpy": "NumPy", "scipy": "SciPy",
    }
    out: dict[str, str] = {}
    for dist, label in wanted.items():
        try:
            out[label] = md.version(dist)
        except Exception:
            continue
    for binary, label in (("mmseqs", "MMseqs2"), ("hmmscan", "HMMER"),
                          ("mkdssp", "DSSP")):
        try:
            proc = subprocess.run([binary, "--version"], capture_output=True,
                                  text=True, timeout=20)
            text = (proc.stdout or proc.stderr).strip().splitlines()
            if text:
                out[label] = text[0][:60]
        except Exception:
            continue
    # boltz-community lives in the other environment.
    try:
        from pipeline.common import ROOT

        python = ROOT / ".pixi" / "envs" / "cofold" / "bin" / "python"
        if python.exists():
            proc = subprocess.run(
                [str(python), "-c",
                 "import importlib.metadata as m;print(m.version('boltz-community'))"],
                capture_output=True, text=True, timeout=60)
            if proc.returncode == 0 and proc.stdout.strip():
                out["boltz-community"] = proc.stdout.strip()
    except Exception:
        pass
    return out


def snapshot_dates() -> dict[str, str]:
    import json
    from pathlib import Path

    out: dict[str, str] = {}
    sabdab = Path("data/sabdab/snapshot.json")
    if sabdab.exists():
        try:
            out["sabdab_nano"] = json.loads(sabdab.read_text())["snapshot_date"]
        except Exception:
            pass
    for name, path in (("uniprot", ANTIGENS / "antigens.json"),
                       ("nbthermo", Path("data/nbthermo/tm_dataset.json"))):
        if path.exists():
            out[name] = dt.datetime.fromtimestamp(
                path.stat().st_mtime).strftime("%Y-%m-%d")
    return out


def model_hashes() -> dict[str, str]:
    import hashlib
    import json

    out: dict[str, str] = {}
    metrics = MODELS / "whittle_tm_v1" / "metrics.json"
    if metrics.exists():
        payload = json.loads(metrics.read_text())
        out["tm model version"] = payload.get("version", "")
        out["tm model ESM-2"] = payload.get("esm_model", "")
        out["tm model grouped R squared"] = str(
            (payload.get("cv") or {}).get("grouped", {}).get("r2", ""))
        out["tm model training set size"] = str(payload.get("n_training", ""))
    dataset = MODELS.parent / "data" / "nbthermo" / "tm_dataset.json"
    if dataset.exists():
        digest = hashlib.sha256(dataset.read_bytes()).hexdigest()[:16]
        out["tm training set hash"] = digest
    return out


def main() -> int:
    cfg = load_yaml("references.yaml")
    limiter = RateLimiter(3)
    conn = connect()
    rows = []
    order = 0
    verified = unverified = 0

    for kind in ("tools", "databases", "benchmarks"):
        for entry in cfg.get(kind) or []:
            order += 1
            doi = (entry.get("doi") or "").strip()
            fallback = entry.get("title") or ""
            resolved: dict = {"doi_verified": 0, "title": fallback}
            if doi:
                resolved = resolve_doi(doi, fallback, limiter)
            commit = (entry.get("commit") or "").strip()
            if commit in ("", "RESOLVE_AT_BUILD"):
                commit = resolve_commit(entry.get("repo", ""))
            snapshot = (entry.get("snapshot_date") or "").strip()
            if snapshot == "RESOLVE_AT_BUILD":
                snapshot = snapshot_dates().get(entry["key"], "")
            if resolved.get("doi_verified"):
                verified += 1
            else:
                unverified += 1
                if resolved.get("verify_note"):
                    log(f"references: {entry['key']} unverified: "
                        f"{resolved['verify_note']}")
            rows.append({
                "key": entry["key"],
                "kind": kind.rstrip("s"),
                "title": resolved.get("title") or fallback,
                "authors": resolved.get("authors", ""),
                "year": resolved.get("year"),
                "container": resolved.get("container", ""),
                "doi": doi or None,
                "doi_verified": resolved.get("doi_verified", 0),
                "url": resolved.get("url") or entry.get("repo") or "",
                "repo": entry.get("repo", ""),
                "commit_sha": commit,
                "licence": entry.get("licence", ""),
                "version_used": tool_versions().get(entry["key"], ""),
                "snapshot_date": snapshot,
                "role": entry.get("role", ""),
                "sort_order": order,
            })

    versions = tool_versions()
    for row in rows:
        for label, version in versions.items():
            if label.lower().startswith(row["key"].replace("_", "-")[:6]):
                row["version_used"] = row["version_used"] or version

    columns = list(rows[0])
    conn.execute('DELETE FROM "references"')
    conn.executemany(
        f'INSERT INTO "references" ({", ".join(columns)}) VALUES '
        f'({", ".join(":" + c for c in columns)})', rows)

    provenance = {"run id": run_id(), "git sha": git_sha(),
                  "built": dt.datetime.now().strftime("%Y-%m-%d %H:%M")}
    provenance.update({f"version: {k}": v for k, v in versions.items()})
    provenance.update({f"snapshot: {k}": v for k, v in snapshot_dates().items()})
    provenance.update(model_hashes())
    conn.execute("DELETE FROM provenance")
    conn.executemany("INSERT INTO provenance (key, value) VALUES (?,?)",
                     sorted(provenance.items()))
    conn.commit()

    # Spec 8.3 allows `repo` plus `commit` as the citable identifier instead of
    # a DOI, so a row with no DOI is not an unverified row: it is a differently
    # cited one. The check is that every row carrying a DOI verifies, and that
    # every row without one carries a repository and a commit instead.
    with_doi = [r for r in rows if r["doi"]]
    doi_failures = [r["key"] for r in with_doi if not r["doi_verified"]]
    repo_cited = [r for r in rows if not r["doi"]]
    no_identifier = [r["key"] for r in repo_cited
                     if not (r["repo"] and r["commit_sha"])]
    conn.execute("INSERT OR REPLACE INTO qc (key, value, passed) VALUES (?,?,?)",
                 ("references_verified",
                  f"{len(with_doi) - len(doi_failures)} of {len(with_doi)} DOIs "
                  f"verified against Crossref"
                  + (f"; unverified: {', '.join(doi_failures)}" if doi_failures else "")
                  + f"; {len(repo_cited)} rows cited by repository and commit"
                  + (f"; NO IDENTIFIER: {', '.join(no_identifier)}"
                     if no_identifier else ""),
                  1 if not doi_failures and not no_identifier else 0))
    conn.commit()
    conn.close()
    manifest_write("references", {"id": "references", "status": "ok",
                                  "rows": len(rows), "verified": verified})
    log(f"references: {len(rows)} rows. {len(with_doi) - len(doi_failures)} of "
        f"{len(with_doi)} DOIs verified against Crossref, "
        f"{len(repo_cited)} rows cited by repository and commit instead "
        f"(spec 8.3 allows both). {len(provenance)} provenance entries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
