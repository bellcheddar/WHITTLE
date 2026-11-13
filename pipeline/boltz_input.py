"""Boltz input generation for Tier 2. Spec 6.4 and 5.2.

Adapted from BoltzMaker's `_chain_entry` / `_build_yaml_doc` shape (see
REUSE.md). The pieces carried over are the chain-entry structure, the
homo-oligomer `id: [A, B, C]` form that shares one sequence across several
chain ids, and the discipline of writing one YAML per complex with a
deterministic stem.

Two things here are WHITTLE's own and follow from spec 5.2:

**TNFa and TL1A are co-folded as homotrimers.** Three identical chains, by
decision, because the epitope commonly spans two subunits and a monomer run
would be wrong rather than merely cheaper. If a trimer will not fit in memory
the fallback is fewer seeds and then skipping the complex, never substituting a
monomer.

**The antigen MSA is cached once per target and reused across every candidate
for that target** (spec 6.4 item 1). That single decision removes most of the
MSA cost. The nanobody chain gets `msa: empty`: the VHH framework is
near-germline and an MSA adds little, which is why the spec says to skip it.
"""
from __future__ import annotations

from pathlib import Path

import yaml

from pipeline.common import ANTIGENS, MSA_CACHE, targets_config

# Chain ids. The antigen takes A (and B, C for a trimer); the nanobody is N, so
# a chain id never has to be reassigned when the stoichiometry changes.
ANTIGEN_CHAINS = ["A", "B", "C"]
NANOBODY_CHAIN = "N"


def antigen_sequence(target: str) -> str:
    path = ANTIGENS / f"{target}.fasta"
    if not path.exists():
        raise FileNotFoundError(
            f"no antigen construct for {target}; run pipeline.antigens first")
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


def antigen_msa(target: str) -> str | None:
    """The cached per-target antigen MSA, or None if it has not been built."""
    for suffix in (".a3m", ".csv"):
        path = MSA_CACHE / f"{target}{suffix}"
        if path.exists() and path.stat().st_size > 0:
            return str(path)
    return None


def chain_entry(chain_id, sequence: str, msa: str | None = None) -> dict:
    entry = {"id": chain_id, "sequence": sequence}
    if msa is not None:
        entry["msa"] = msa
    return {"protein": entry}


def build_doc(target: str, candidate_id: str, nanobody: str) -> dict:
    """One Boltz YAML document for one nanobody against one target."""
    spec = targets_config()["targets"][target]
    stoich = int(spec["stoichiometry"])
    if stoich > len(ANTIGEN_CHAINS):
        raise ValueError(f"{target} asks for {stoich} chains; only "
                         f"{len(ANTIGEN_CHAINS)} chain ids are reserved")
    antigen = antigen_sequence(target)
    msa = antigen_msa(target)

    sequences = []
    # A homotrimer is one entry with three chain ids, not three entries: Boltz
    # shares the MSA and the embedding across them, which is both correct and
    # much cheaper than three independent chains.
    ids = ANTIGEN_CHAINS[:stoich] if stoich > 1 else ANTIGEN_CHAINS[0]
    sequences.append(chain_entry(ids, antigen, msa if msa else "empty"))
    # The nanobody: no MSA. The VHH framework is near-germline, so an MSA adds
    # little and would cost an MSA-server round trip per candidate (spec 6.4).
    sequences.append(chain_entry(NANOBODY_CHAIN, nanobody, "empty"))
    return {"version": 1, "sequences": sequences}


def stem(target: str, candidate_id: str) -> str:
    return f"{target}__{candidate_id}"


def write_input(target: str, candidate_id: str, nanobody: str, out_dir: Path
                ) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{stem(target, candidate_id)}.yaml"
    doc = build_doc(target, candidate_id, nanobody)
    with path.open("w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, default_flow_style=False)
    return path


def complex_length(target: str, nanobody: str) -> int:
    spec = targets_config()["targets"][target]
    return len(antigen_sequence(target)) * int(spec["stoichiometry"]) + len(nanobody)


def describe(target: str) -> dict:
    spec = targets_config()["targets"][target]
    antigen = antigen_sequence(target)
    stoich = int(spec["stoichiometry"])
    return {
        "target": target,
        "construct": spec["construct"],
        "range": spec["range"],
        "chain_length": len(antigen),
        "stoichiometry": stoich,
        "antigen_length": len(antigen) * stoich,
        "complex_length_estimate": len(antigen) * stoich + 125,
        "msa_cached": antigen_msa(target) is not None,
        "trimer": stoich > 1,
        "ecd_scoped": spec["construct"] == "ecd_only",
    }
