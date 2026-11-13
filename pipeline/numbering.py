"""IMGT numbering and CDR definition. Spec 6.1 items 1-3, and the 0.2 fallback ladder.

Ladder, in the order spec 0.2 sets out:

1. **ANARCI** (`anarci` module, HMMER-backed). The house default.
2. **anarcii** (PyTorch, no HMMER). Same IMGT scheme, different engine.
3. **Regex against the VHH framework consensus.** Marks the row
   `numbering_method = "regex_fallback"` so nothing downstream reads it as an
   alignment result.

Whichever route runs, this module returns the same shape, so Tier 0 does not
know or care which one answered.

IMGT CDR boundaries (the definition this project uses throughout):
CDR1 27-38, CDR2 56-65, CDR3 105-117. The VHH hallmark tetrad sits at IMGT
42, 49, 50 and 52 (spec 6.1 item 3).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache

IMGT_CDR = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
# Framework-2 positions that distinguish a VHH from a VH. The camelid hallmark
# set is the one in the literature: 42F/Y, 49E, 50R, 52G/L/F. A VH-type residue
# at any of them (42V, 49G, 50L, 52W) is a solubility and stability liability
# because it leaves the former VL interface hydrophobic and exposed.
HALLMARK_POSITIONS = [42, 49, 50, 52]
HALLMARK_VHH = {42: set("FYHL"), 49: set("EDQ"), 50: set("RCLKQ"), 52: set("GLFADS")}
HALLMARK_VH = {42: set("V"), 49: set("GA"), 50: set("LMV"), 52: set("WY")}


@dataclass
class Numbering:
    method: str
    scheme: str = "imgt"
    chain_type: str = ""
    v_gene: str = ""
    v_identity: float | None = None
    # position -> single-letter residue, insertion codes collapsed onto the
    # integer position in order (IMGT 111A, 112A and so on inside CDR3)
    numbered: list[tuple[int, str, str]] = field(default_factory=list)
    regions: dict[str, str] = field(default_factory=dict)
    ok: bool = True
    detail: str = ""

    def at(self, position: int) -> str:
        """The residue at an IMGT position, or '' if the position is absent."""
        for pos, ins, aa in self.numbered:
            if pos == position and ins == " ":
                return aa
        return ""


# --------------------------------------------------------------------- route 1: ANARCI

@lru_cache(maxsize=1)
def _anarci():
    try:
        from anarci import run_anarci  # type: ignore

        return run_anarci
    except Exception:
        return None


def by_anarci(sequence: str) -> Numbering | None:
    run_anarci = _anarci()
    if run_anarci is None:
        return None
    try:
        # assign_germline is what puts the V gene in `details`; without it the
        # key is simply absent. allowed_species=None keeps alpaca and camel in,
        # which is the whole point for a VHH.
        _, numbered, details, _ = run_anarci(
            [("q", sequence)], scheme="imgt", allowed_species=None,
            assign_germline=True)
    except Exception as exc:
        return Numbering(method="anarci", ok=False, detail=f"{type(exc).__name__}: {exc}")
    if not numbered or not numbered[0]:
        return Numbering(method="anarci", ok=False, detail="no domain found")
    domain = numbered[0][0][0]
    info = (details[0] or [{}])[0] or {}
    rows: list[tuple[int, str, str]] = []
    for (pos, ins), aa in domain:
        if aa == "-":
            continue
        rows.append((int(pos), ins or " ", aa))
    v_gene, v_identity = "", None
    germlines = info.get("germlines") or {}
    entry = germlines.get("v_gene")
    if entry:
        # ANARCI returns (("species", "IGHV3-3"), identity).
        try:
            v_gene = str(entry[0][1])
            v_identity = float(entry[1])
        except (IndexError, TypeError, ValueError):
            v_gene, v_identity = str(entry), None
    return Numbering(
        method="anarci",
        chain_type=str(info.get("chain_type", "")),
        v_gene=v_gene,
        v_identity=v_identity,
        numbered=rows,
    )


# -------------------------------------------------------------------- route 2: anarcii

@lru_cache(maxsize=1)
def _anarcii():
    try:
        from anarcii import Anarcii  # type: ignore

        model = Anarcii(seq_type="antibody", batch_size=64, cpu=True, verbose=False,
                        mode="accuracy")
        return model
    except Exception:
        return None


def by_anarcii_batch(sequences: list[tuple[str, str]]) -> dict[str, Numbering]:
    """anarcii is a neural model: batching it is the only way it is not slow.

    Returns {name: Numbering}. A sequence it cannot number is absent from the
    result, and the caller drops to the regex route for it.
    """
    model = _anarcii()
    if model is None:
        return {}
    out: dict[str, Numbering] = {}
    try:
        model.number(dict(sequences))
        scheme_out = model.to_scheme("imgt")
    except Exception:
        return {}
    for name, result in scheme_out.items():
        try:
            numbering = result.get("numbering") if isinstance(result, dict) else None
            if not numbering:
                continue
            rows: list[tuple[int, str, str]] = []
            for entry in numbering:
                (pos, ins), aa = entry
                if aa in ("-", "X"):
                    continue
                rows.append((int(pos), ins or " ", aa))
            if not rows:
                continue
            out[str(name)] = Numbering(
                method="anarcii",
                chain_type=str(result.get("chain_type", "")) if isinstance(result, dict) else "",
                v_gene="",
                numbered=rows,
            )
        except Exception:
            continue
    return out


# --------------------------------------------------------------- route 3: regex fallback

# The VHH framework consensus, as anchors rather than as a whole alignment.
# FR1 ends at the conserved C at IMGT 23 preceded by ...SC(A/T)AS; FR2 opens at
# W(F/V)RQ; FR3 opens at the K/R-G-R/K-F/Y-T/A-I/L-S block; FR4 opens at WG.GT.
RX_FR1_CYS = re.compile(r"C[A-Z]{2}[AS]")           # Cys22 then the AAS/TAS motif
RX_FR2 = re.compile(r"W[FVYL]RQ")                    # Trp36 .. the FR2 opener
RX_FR3 = re.compile(r"R[FY][AGTSIV][IVT][SLT]RD")    # the FR3 opener after CDR2
RX_FR4 = re.compile(r"WG[QKRLE]G")                   # Trp103 then the FR4 opener
RX_CDR3_CYS = re.compile(r"[YF][YFHC]C")             # the Cys96 anchor before CDR3


def by_regex(sequence: str) -> Numbering:
    """Derive the three CDRs from framework anchors. Spec 0.2, last resort.

    This is deliberately conservative: when an anchor is missing the region is
    left empty rather than guessed, and `detail` says which anchor failed.
    """
    seq = sequence.strip().upper()
    missing: list[str] = []

    m_cys1 = RX_FR1_CYS.search(seq[:35])
    m_fr2 = RX_FR2.search(seq[25:60])
    m_fr3 = RX_FR3.search(seq[55:95])
    m_cys3 = None
    for match in RX_CDR3_CYS.finditer(seq[80:]):
        m_cys3 = match  # the last one before FR4 is the Cys96 anchor
    m_fr4 = RX_FR4.search(seq[85:])

    cdr1 = cdr2 = cdr3 = ""
    # Offsets are counted off the conserved anchors against an IMGT-numbered VHH,
    # not guessed: the W of the W(F/V)RQ block is IMGT 41, so IMGT 38 (the last
    # CDR1 residue) sits three before it, and IMGT 56 (the first CDR2 residue)
    # sits fifteen after it. The eight residues between IMGT 65 and the FR3
    # R-F-T-I-S-R-D opener are framework.
    w_abs = (25 + m_fr2.start()) if m_fr2 else None
    if m_cys1 and w_abs is not None:
        cdr1 = seq[m_cys1.end(): w_abs - 2]
    else:
        missing.append("cdr1 anchors")
    if w_abs is not None and m_fr3:
        start = w_abs + 15
        end = 55 + m_fr3.start() - 8
        cdr2 = seq[start:end] if 0 < start < end else ""
        if not cdr2:
            missing.append("cdr2 span")
    else:
        missing.append("cdr2 anchors")
    if m_cys3 and m_fr4:
        cdr3 = seq[80 + m_cys3.end(): 85 + m_fr4.start()]
    else:
        missing.append("cdr3 anchors")

    # Synthesise IMGT positions for the hallmark tetrad only: the tetrad sits in
    # FR2, a fixed offset from the W(F/V)RQ anchor (IMGT 41 is that W, so 42 is
    # the next residue and 49, 50 and 52 follow at +7, +8 and +10).
    rows: list[tuple[int, str, str]] = []
    if w_abs is not None:
        for offset, pos in ((1, 42), (8, 49), (9, 50), (11, 52)):
            idx = w_abs + offset
            if idx < len(seq):
                rows.append((pos, " ", seq[idx]))

    return Numbering(
        method="regex_fallback",
        chain_type="H",
        numbered=sorted(rows),
        regions={"cdr1": cdr1, "cdr2": cdr2, "cdr3": cdr3},
        ok=not missing,
        detail=("missing: " + ", ".join(missing)) if missing else "",
    )


# ------------------------------------------------------------------------- the ladder

def regions_from_numbering(n: Numbering) -> dict[str, str]:
    if n.regions:
        return n.regions
    out: dict[str, str] = {}
    for name, (lo, hi) in IMGT_CDR.items():
        out[name] = "".join(aa for pos, _ins, aa in n.numbered if lo <= pos <= hi)
    return out


def hallmark(n: Numbering) -> tuple[int, int, str]:
    """(intact, n_camelid, detail). Spec 6.1 item 3.

    `intact` is strict: all four of IMGT 42, 49, 50 and 52 carry a camelid
    residue. `n_camelid` is the graded 0-4 count, which is what the UI shows,
    because "three of four, with a VH-type Trp at 52" is the useful statement
    and a bare boolean throws that away. `detail` names each position and marks
    a VH-type residue with `!VH`, which is the solubility and stability
    liability the spec asks to see.
    """
    parts: list[str] = []
    n_camelid = 0
    for pos in HALLMARK_POSITIONS:
        aa = n.at(pos)
        if not aa:
            parts.append(f"{pos}:?")
            continue
        if aa in HALLMARK_VHH[pos]:
            parts.append(f"{pos}{aa}")
            n_camelid += 1
        elif aa in HALLMARK_VH[pos]:
            parts.append(f"{pos}{aa}!VH")
        else:
            parts.append(f"{pos}{aa}?")
    return (1 if n_camelid == len(HALLMARK_POSITIONS) else 0), n_camelid, " ".join(parts)


def number_many(sequences: list[tuple[str, str]]) -> dict[str, Numbering]:
    """Number a batch down the ladder. Returns {name: Numbering} for every input."""
    out: dict[str, Numbering] = {}
    run_anarci = _anarci()
    if run_anarci is not None:
        for name, seq in sequences:
            result = by_anarci(seq)
            if result is not None and result.ok and result.numbered:
                result.regions = regions_from_numbering(result)
                out[name] = result
    remaining = [(n, s) for n, s in sequences if n not in out]
    if remaining:
        batched = by_anarcii_batch(remaining)
        for name, result in batched.items():
            result.regions = regions_from_numbering(result)
            out[name] = result
    for name, seq in sequences:
        if name not in out:
            out[name] = by_regex(seq)
    return out


def available_methods() -> dict[str, bool]:
    """Which routes are installed. Deliberately does not instantiate anarcii:
    building its model loads torch weights, and the fallback is only reached if
    ANARCI has already failed."""
    try:
        import anarcii  # noqa: F401

        has_anarcii = True
    except Exception:
        has_anarcii = False
    return {"anarci": _anarci() is not None, "anarcii": has_anarcii,
            "regex_fallback": True}
