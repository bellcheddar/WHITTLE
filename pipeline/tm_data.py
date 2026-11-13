"""Training labels for the thermostability model. Spec 6.3 item 1.

Two public sources, both cited to CC-BY articles, both stored only as derived
values in the repository (the raw downloads live under the gitignored
`data/nbthermo/`):

1. **The NanoMelt curated set**, Supplementary Dataset 1 of Ramon, Sormanni et
   al. (2025) mAbs, doi:10.1080/19420862.2024.2442750. 640 nanobody sequences
   with a measured apparent melting temperature and the experimental method,
   from the journal's CC-BY supplementary files. The NanoMelt code itself is
   CC-BY-NC-SA and is neither installed nor used (spec 4.2).
2. **NbThermo**, Valdés-Tresanco et al. (2023) Database, doi:10.1093/database/baad021.
   548 entries; 519 carry an AHo-numbered sequence. `Sequence.Raw` is populated
   for only 44 entries, so the loader reads `Sequence.Aho` and strips the gap
   characters: reading `Raw` instead silently yields 44 rows.

Apparent Tm depends on the method and the buffer, so the method is carried
through as a column and the cross-validation stratifies on it. Where a source
gives several measurements for one sequence the median is taken and the spread
is recorded, because that spread is the floor on any honest error bar: the
NanoMelt authors report a mean pairwise difference of 2.4 degrees C between
repeat measurements of the same nanobody.

    pixi run python -m pipeline.tm_data
"""
from __future__ import annotations

import json
import re
import statistics
import sys
from pathlib import Path

from pipeline.common import DATA, log, write_json

NBTHERMO = DATA / "nbthermo"
VALID = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

# NbThermo records a Tm per experimental method. The keys are its own.
NBTHERMO_METHODS = ["nanoDSF", "DSF (SYPRO)", "DSC", "Circular dichroism", "Other"]


def _clean(sequence: str | None) -> str:
    if not sequence:
        return ""
    seq = str(sequence).upper().replace("-", "").replace(".", "").strip()
    seq = "".join(c for c in seq if c.isalpha())
    return seq if VALID.match(seq) else ""


def load_nanomelt() -> list[dict]:
    """Supplementary Dataset 1: 640 curated apparent-Tm nanobodies."""
    import pandas as pd

    path = NBTHERMO / "KMAB_A_2442750_SM6729.xlsx"
    if not path.exists():
        log(f"tm_data: {path.name} is absent; the NanoMelt set is skipped")
        return []
    frame = pd.read_excel(path)
    frame.columns = [str(c).strip() for c in frame.columns]
    seq_col = next(c for c in frame.columns if c.lower().startswith("sequence"))
    tm_col = next(c for c in frame.columns if "melting" in c.lower())
    method_col = next((c for c in frame.columns if "method" in c.lower()), None)
    name_col = next((c for c in frame.columns if c.lower() == "name"), None)
    source_col = next((c for c in frame.columns if c.lower() == "source"), None)

    out: list[dict] = []
    for _, row in frame.iterrows():
        seq = _clean(row[seq_col])
        try:
            tm = float(row[tm_col])
        except (TypeError, ValueError):
            continue
        if not seq or not (0 < tm < 140):
            continue
        out.append({
            "sequence": seq,
            "tm": round(tm, 2),
            "method": str(row[method_col]).strip() if method_col else "unknown",
            "name": str(row[name_col]).strip() if name_col else "",
            "origin": str(row[source_col]).strip() if source_col else "",
            "source": "nanomelt_sd1",
        })
    log(f"tm_data: NanoMelt Supplementary Dataset 1 gave {len(out)} labelled rows")
    return out


def load_nanomelt_external() -> list[dict]:
    """Supplementary Dataset 2: 83 nanobodies held out from NanoMelt's training.

    Kept apart as an external test set: a model trained on Dataset 1 and scored
    on Dataset 2 is the only honest statement this build can make about
    generalisation to sequences nobody fitted on.
    """
    import pandas as pd

    path = NBTHERMO / "KMAB_A_2442750_SM6730.xlsx"
    if not path.exists():
        return []
    # The sheet carries a title row and two blank rows above the real header,
    # so the header row is located rather than assumed.
    probe = pd.read_excel(path, header=None)
    header_row = None
    for i in range(min(12, len(probe))):
        values = [str(v).strip().lower() for v in probe.iloc[i].tolist()]
        if "sequence" in values and any(v.startswith("measured tm") for v in values):
            header_row = i
            break
    if header_row is None:
        log("tm_data: the external sheet's header could not be located; skipped")
        return []
    frame = pd.read_excel(path, header=header_row)
    frame.columns = [str(c).strip() for c in frame.columns]
    seq_col = next((c for c in frame.columns if c.lower().startswith("sequence")), None)
    tm_col = next((c for c in frame.columns if c.lower().startswith("measured tm")), None)
    method_col = next((c for c in frame.columns if "techni" in c.lower()), None)
    if seq_col is None or tm_col is None:
        log("tm_data: the external sheet's columns could not be located; skipped")
        return []
    out: list[dict] = []
    for _, row in frame.iterrows():
        seq = _clean(row[seq_col])
        try:
            tm = float(row[tm_col])
        except (TypeError, ValueError):
            continue
        if not seq or not (0 < tm < 140):
            continue
        out.append({"sequence": seq, "tm": round(tm, 2),
                    "method": str(row[method_col]).strip() if method_col else "unknown",
                    "name": "", "origin": "", "source": "nanomelt_sd2_external"})
    log(f"tm_data: NanoMelt Supplementary Dataset 2 gave {len(out)} external rows")
    return out


def load_nbthermo() -> list[dict]:
    """NbThermo. Read Sequence.Aho and strip gaps, not Sequence.Raw."""
    path = NBTHERMO / "nbthermo_database.json"
    if not path.exists():
        log(f"tm_data: {path.name} is absent; NbThermo is skipped")
        return []
    entries = json.loads(path.read_text())
    out: list[dict] = []
    raw_only = 0
    for entry in entries:
        seq_block = entry.get("Sequence") or {}
        seq = _clean(((seq_block.get("Aho") or {}).get("value")))
        if not seq:
            seq = _clean(((seq_block.get("Raw") or {}).get("value")))
            if seq:
                raw_only += 1
        if not seq:
            continue
        tm_block = entry.get("Tm") or {}
        for method in NBTHERMO_METHODS:
            value = (tm_block.get(method) or {}).get("value")
            try:
                tm = float(value)
            except (TypeError, ValueError):
                continue
            if not (0 < tm < 140):
                continue
            out.append({
                "sequence": seq, "tm": round(tm, 2), "method": method,
                "name": str(((entry.get("Info") or {}).get("Name") or {}).get("value", "")),
                "origin": str(((entry.get("Origin") or {}).get("Source") or {}).get("value", "")),
                "doi": str(((entry.get("Info") or {}).get("Reference") or {}).get("DOI", "")),
                "source": "nbthermo",
            })
    log(f"tm_data: NbThermo gave {len(out)} measurement rows across "
        f"{len({r['sequence'] for r in out})} sequences "
        f"({raw_only} fell back to the Raw field)")
    return out


def collapse(rows: list[dict]) -> list[dict]:
    """One row per sequence: the median Tm, with the spread kept.

    The spread is not a nuisance. It is the measurement floor for the error
    bar, and a sequence measured three ways with a 9 degree range should not be
    presented as a crisp training label.
    """
    by_sequence: dict[str, list[dict]] = {}
    for row in rows:
        by_sequence.setdefault(row["sequence"], []).append(row)
    out: list[dict] = []
    for seq, group in by_sequence.items():
        values = [g["tm"] for g in group]
        spread = (max(values) - min(values)) if len(values) > 1 else 0.0
        methods = sorted({g["method"] for g in group})
        out.append({
            "sequence": seq,
            "tm": round(statistics.median(values), 2),
            "tm_spread": round(spread, 2),
            "n_measurements": len(values),
            "methods": ";".join(methods),
            "primary_method": methods[0],
            "name": next((g["name"] for g in group if g.get("name")), ""),
            "sources": ";".join(sorted({g["source"] for g in group})),
            "length": len(seq),
        })
    return sorted(out, key=lambda r: r["sequence"])


def build() -> dict:
    training = collapse(load_nanomelt() + load_nbthermo())
    external_raw = load_nanomelt_external()
    external = collapse(external_raw)
    # A sequence in both is not external any more: drop it from the held-out
    # set rather than quietly scoring the model on something it was fitted on.
    train_sequences = {r["sequence"] for r in training}
    overlap = [r for r in external if r["sequence"] in train_sequences]
    external = [r for r in external if r["sequence"] not in train_sequences]

    payload = {
        "training": training,
        "external": external,
        "n_training": len(training),
        "n_external": len(external),
        "n_external_dropped_as_overlapping": len(overlap),
        "tm_range": [min((r["tm"] for r in training), default=None),
                     max((r["tm"] for r in training), default=None)],
        # Most repeats are the same measurement reaching both sources, so the
        # median spread is 0 and says nothing. The spread over sequences that
        # genuinely disagree is the number that bounds the error bar.
        "median_spread_of_repeated_sequences": round(statistics.median(
            [r["tm_spread"] for r in training if r["n_measurements"] > 1] or [0.0]), 2),
        "mean_spread_where_measurements_disagree": round(statistics.mean(
            [r["tm_spread"] for r in training if r["tm_spread"] > 0] or [0.0]), 2),
        "max_spread": round(max(
            [r["tm_spread"] for r in training] or [0.0]), 2),
        "n_with_repeats": sum(1 for r in training if r["n_measurements"] > 1),
        "n_where_measurements_disagree": sum(1 for r in training if r["tm_spread"] > 0),
        "sources": {
            "nanomelt": {
                "doi": "10.1080/19420862.2024.2442750",
                "what": "Supplementary Dataset 1, 640 curated apparent-Tm nanobodies",
                "licence": "CC BY, as the article",
                "note": "The NanoMelt model itself is CC-BY-NC-SA and is not used.",
            },
            "nbthermo": {
                "doi": "10.1093/database/baad021",
                "what": "NbThermo, 548 entries; AHo sequences with gaps stripped",
                "licence": "the repository carries no licence file; cited to the "
                           "CC-BY article, derived values only",
            },
        },
    }
    write_json(NBTHERMO / "tm_dataset.json", payload)
    log(f"tm_data: {payload['n_training']} unique training sequences, "
        f"{payload['n_external']} external ({payload['n_external_dropped_as_overlapping']} "
        f"dropped as overlapping), Tm {payload['tm_range'][0]} to "
        f"{payload['tm_range'][1]} degrees C, "
        f"{payload['n_with_repeats']} sequences measured more than once, of which "
        f"{payload['n_where_measurements_disagree']} disagree between measurements "
        f"by a mean of {payload['mean_spread_where_measurements_disagree']} and at "
        f"most {payload['max_spread']} degrees C")
    return payload


def load() -> dict:
    path = NBTHERMO / "tm_dataset.json"
    if path.exists():
        return json.loads(path.read_text())
    return build()


if __name__ == "__main__":
    build()
    sys.exit(0)
