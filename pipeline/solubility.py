"""Sequence-based solubility: SWI and Protein-Sol. Spec 6.2 item 4.

Two published predictors, reimplemented here from their published
specifications so that nothing has to be called at serve time and nothing
depends on a web server.

**SWI** (Solubility-Weighted Index). Bhandari, Gardner and Lim (2020),
Bioinformatics 36, doi:10.1093/bioinformatics/btaa578. SWI is the mean of 20
published per-residue weights over the sequence, and the probability of
solubility is a logistic function of it with two published constants. Both the
weights and the constants are from the authors' own reference implementation at
github.com/Gardner-BinfLab/SoDoPE_paper_2020.

**Protein-Sol**. Hebditch, Carballo-Amador, Charonis, Curtis and Warwicker
(2017), Bioinformatics 33, doi:10.1093/bioinformatics/btx345. The paper does
not publish the model coefficients; the authors' downloadable reference package
does, in `seq_reference_data.txt`, together with a plain-English statement of
the algorithm. That package is GPL-3 (Copyright 2017 Jim Warwicker and Max
Hebditch), which spec 4.2 anticipates for this component. The code below is an
independent Python implementation written from that specification, not a
translation of their perl; the reference statistics it carries are the
published fitting constants, attributed here and in the reference table.

The algorithm, as the reference package states it: ten features are used, each
clamped to the interval between its average in the least-soluble and the
most-soluble tails of the Niwa E. coli cell-free dataset, linearly interpolated
onto the solubility values of those two tails, and combined as a weighted mean
with weights equal to the normalised absolute z-score difference per feature.

Both implementations are checked against the authors' own worked examples in
`pipeline/tests/test_solubility.py`. That test is the point: a solubility
column that silently disagrees with the tool it claims to be is worse than no
column at all.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# --------------------------------------------------------------------------- SWI
# Bhandari, Gardner and Lim (2020), doi:10.1093/bioinformatics/btaa578.
SWI_WEIGHTS = {
    "A": 0.8356471476582918, "C": 0.5208088354857734, "E": 0.9876987431418378,
    "D": 0.9079044671339564, "G": 0.7997168496420723, "F": 0.5849790194237692,
    "I": 0.6784124413866582, "H": 0.8947913996466419, "K": 0.9267104557513497,
    "M": 0.6296623675420369, "L": 0.6554221515081433, "N": 0.8597433107431216,
    "Q": 0.7894346483482080, "P": 0.8235328714705341, "S": 0.7440908318492778,
    "R": 0.7712466317693457, "T": 0.8096922697856334, "W": 0.6374678690957594,
    "V": 0.7357837119163659, "Y": 0.6112801822947587,
}
# Logistic fit published with the weights: prob = 1 / (1 + exp(-(A * swi + B)))
SWI_A = 81.0581
SWI_B = -62.7775


def swi(sequence: str) -> tuple[float, float]:
    """(SWI, probability of solubility). Selenocysteine is read as cysteine,
    as the reference implementation does; anything else unknown is skipped."""
    seq = sequence.strip().upper().replace("U", "C")
    values = [SWI_WEIGHTS[aa] for aa in seq if aa in SWI_WEIGHTS]
    if not values:
        return float("nan"), float("nan")
    index = sum(values) / len(values)
    probability = 1.0 / (1.0 + math.exp(-(SWI_A * index + SWI_B)))
    return round(index, 5), round(probability, 5)


# -------------------------------------------------------------------- Protein-Sol
# Fitting constants from the Protein-Sol reference package
# (protein-sol.manchester.ac.uk, GPL-3, Copyright 2017 Warwicker and Hebditch),
# derived from the Niwa et al. (2009) PNAS 106:4201 E. coli cell-free dataset.
#
# Solubility of the two tails, as a percentage.
PS_SOL_LOW = 5.208
PS_SOL_TOP = 113.241
PS_SOL_POP = 53.347

# feature -> (average in the LOW tail, average in the TOP tail, weight)
# The weights are the normalised absolute z-score differences and sum to 1.
PS_FEATURES: dict[str, tuple[float, float, float]] = {
    "K-R":        (-2.651,    1.042,   0.119),
    "length":     (437.096, 188.055,   0.166),
    "H":          (  2.923,    2.259,  0.053),
    "L":          ( 10.676,    9.202,  0.061),
    "V":          (  6.001,    7.112,  0.057),
    "D+E":        ( 11.395,   13.657,  0.095),
    "F+W+Y":      (  8.918,    6.447,  0.101),
    "abs_charge": (  0.0160,   0.0470, 0.120),
    "fold_index": (  0.1500,   0.0890, 0.077),
    "entropy":    (  4.113,    3.995,  0.150),
}

# Kyte and Doolittle (1982) JMB 157:105, normalised to 0-1 as FoldIndex uses it.
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


@dataclass
class ProteinSol:
    percent: float          # predicted solubility as a percentage
    scaled: float           # the 0-1 value the server plots
    population: float       # the dataset population average on the same scale
    features: dict[str, float]


def _composition(seq: str) -> dict[str, float]:
    """Percentage composition, which is what the reference statistics are in."""
    n = len(seq)
    return {aa: 100.0 * seq.count(aa) / n for aa in AMINO_ACIDS} if n else {}


def protein_sol_features(sequence: str) -> dict[str, float]:
    seq = "".join(c for c in sequence.strip().upper() if c in AMINO_ACIDS)
    n = len(seq)
    if not n:
        return {}
    comp = _composition(seq)
    # Net charge from K, R, D and E only, which is what FoldIndex and the
    # reference statistics use; per residue, and taken as an absolute value.
    net = (seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")) / n
    abs_charge = abs(net)
    # Kyte-Doolittle mean, normalised onto 0-1 across the scale's -4.5 to 4.5.
    kd = sum(KD_HYDROPATHY[aa] for aa in seq) / n
    kd_norm = (kd + 4.5) / 9.0
    # Uversky et al. (2000) Proteins 41:415, the FoldIndex form.
    fold_index = 2.785 * kd_norm - abs_charge - 1.151
    # Shannon entropy of the composition, in bits.
    entropy = 0.0
    for aa in AMINO_ACIDS:
        p = comp[aa] / 100.0
        if p > 0:
            entropy -= p * math.log2(p)
    return {
        "K-R": comp["K"] - comp["R"],
        "length": float(n),
        "H": comp["H"],
        "L": comp["L"],
        "V": comp["V"],
        "D+E": comp["D"] + comp["E"],
        "F+W+Y": comp["F"] + comp["W"] + comp["Y"],
        "abs_charge": abs_charge,
        "fold_index": fold_index,
        "entropy": entropy,
        # Carried for the record, not used by the ten-feature model.
        "kyte_doolittle_normalised": kd_norm,
        "net_charge_per_residue": net,
    }


def protein_sol(sequence: str) -> ProteinSol:
    features = protein_sol_features(sequence)
    if not features:
        return ProteinSol(float("nan"), float("nan"), float("nan"), {})
    total = 0.0
    for name, (low, top, weight) in PS_FEATURES.items():
        x = features[name]
        # Clamp into the interval the two tails span. The reference text writes
        # this as "if x < xL, x = xL or if x > xT, x = xT", which assumes
        # xT > xL; length, H, L, F+W+Y, fold index and entropy all run the
        # other way, so the clamp is written against the interval itself.
        lo, hi = (low, top) if low <= top else (top, low)
        x = min(max(x, lo), hi)
        fraction = (x - low) / (top - low) if top != low else 0.0
        total += weight * (PS_SOL_LOW + fraction * (PS_SOL_TOP - PS_SOL_LOW))
    span = PS_SOL_TOP - PS_SOL_LOW
    return ProteinSol(
        percent=round(total, 3),
        scaled=round((total - PS_SOL_LOW) / span, 4),
        population=round((PS_SOL_POP - PS_SOL_LOW) / span, 4),
        features={k: round(v, 4) for k, v in features.items()},
    )


def both(sequence: str) -> dict[str, float]:
    """What Tier 1 stores: SWI, its probability, and the Protein-Sol 0-1 value."""
    index, probability = swi(sequence)
    ps = protein_sol(sequence)
    return {
        "swi": index,
        "swi_prob_soluble": probability,
        "protein_sol": ps.scaled,
        "protein_sol_percent": ps.percent,
        "protein_sol_population": ps.population,
    }


if __name__ == "__main__":
    import sys

    for seq in sys.argv[1:]:
        print(seq[:30], both(seq))
