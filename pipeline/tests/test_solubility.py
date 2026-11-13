"""SWI and Protein-Sol checked against their authors' own worked examples.

The fixtures below are taken from the reference packages, not from this
implementation, so this test can fail. That is the point: a solubility column
that silently disagrees with the tool it claims to be is worse than no column.

Protein-Sol fixture: `seq_prediction.txt` shipped in the authors' downloadable
sequence-prediction package (protein-sol.manchester.ac.uk), for the two
sequences in its own `example.fasta`.
"""
from __future__ import annotations

import math

import pytest

from pipeline.solubility import protein_sol, protein_sol_features, swi

# E. coli homoserine kinase, UniProt P00547, the first entry in the authors'
# example.fasta.
P00547 = (
    "MVKVYAPASSANMSVGFDVLGAAVTPVDGALLGDVVTVEAAETFSLNNLGRFADKLPSEPRENIVYQCWERFCQE"
    "LGKQIPVAMTLEKNMPIGSGLGSSACSVVAALMAMNEHCGKPLNDTRLLALMGELEGRISGSIHYDNVAPCFLGG"
    "MQLMIEENDIISQQVPGFDEWLWVLAYPGIKVSTAEARAILPAQYRRQDCIAHGRHLAGFIHACYSRQPELAAKL"
    "MKDVIAEPYRERLLPGFRQARQAVAEIGAVASGISGSGPTLFALCDKPETAQRVADWLGKNYLQNQEGFVHICRL"
    "DTAGARVLEN"
)
# Human lysozyme C, UniProt P61626, the second entry.
P61626 = (
    "MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIF"
    "QINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV"
)

# SEQUENCE PREDICTIONS lines from the authors' seq_prediction.txt:
#   >P00547              41.535, 0.336, 0.446
#   >sp|P61626|LYSC...   64.448, 0.548, 0.446
EXPECTED = {
    P00547: {"percent": 41.535, "scaled": 0.336},
    P61626: {"percent": 64.448, "scaled": 0.548},
}


@pytest.mark.parametrize("sequence", list(EXPECTED))
def test_protein_sol_matches_reference_output(sequence: str) -> None:
    result = protein_sol(sequence)
    expected = EXPECTED[sequence]
    assert result.percent == pytest.approx(expected["percent"], abs=0.3), (
        f"predicted {result.percent} percent, the reference package says "
        f"{expected['percent']}")
    assert result.scaled == pytest.approx(expected["scaled"], abs=0.005)


def test_protein_sol_population_average_is_the_published_one() -> None:
    # The 0.45 that the server plots as the population average.
    assert protein_sol(P00547).population == pytest.approx(0.446, abs=0.002)


def test_protein_sol_features_are_internally_consistent() -> None:
    f = protein_sol_features(P00547)
    # FoldIndex is Uversky's expression in the normalised hydropathy and the
    # absolute net charge per residue. If this drifts, the fold-index feature is
    # no longer the one the reference statistics were fitted against.
    expected = 2.785 * f["kyte_doolittle_normalised"] - f["abs_charge"] - 1.151
    assert f["fold_index"] == pytest.approx(expected, abs=1e-6)
    # Shannon entropy of a 20-letter composition cannot exceed log2(20).
    assert 0 < f["entropy"] <= math.log2(20)
    assert f["length"] == len(P00547)


def test_swi_is_the_mean_of_the_published_weights() -> None:
    index, probability = swi(P00547)
    # Every published weight lies in (0.5, 1.0), so their mean must too.
    assert 0.5 < index < 1.0
    assert 0.0 <= probability <= 1.0
    # A poly-lysine sequence is at the top of the weight table and a poly-
    # cysteine one at the bottom: the index has to order them that way round.
    assert swi("K" * 50)[0] > index > swi("C" * 50)[0]


def test_swi_handles_selenocysteine_as_cysteine() -> None:
    assert swi("UUUU")[0] == swi("CCCC")[0]


def test_empty_sequence_does_not_raise() -> None:
    assert math.isnan(swi("")[0])
    assert math.isnan(protein_sol("").percent)
