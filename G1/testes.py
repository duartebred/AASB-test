import pytest

from resolucao import cols_conservadas, prot2re, maior_inversao
from resolucao import perc_conservacao


def test_cols_conservadas_ignore_gap_true():
    alin = ['ATGA-AC',
            'AA-AT-C',
            'ATCA-AG',
            'ATGATAG',
            'ATCTCCG']
    assert cols_conservadas(alin) == [0, 1, 3, 5]


def test_cols_conservadas_ignore_gap_false():
    alin = ['ATGA-AC',
            'AA-AT-C',
            'ATCA-AG',
            'ATGATAG',
            'ATCTCCG']
    assert cols_conservadas(alin, ignore_gap=False) == [0, 1, 3]


@pytest.mark.parametrize("perc, expected", [
    (3, 80),
    (2, 80),
    (4, 75),
    (1, 100),
    (5, 67)
])
def test_perc_conservacao_ignore_gap_true(perc, expected):
    alin = ['ATGA-AC',
            'AA-AT-C',
            'ATCA-AG',
            'ATGATAG',
            'ATCTCCG']

    # 1 0.8 0.5 0.8 0.66 0.75 0.6
    # ordenar
    # 0.5 0.6 0.67 0.75 0.8 0.8 1
    result = perc_conservacao(alin, perc)
    assert result == expected


@pytest.mark.parametrize("perc, expected", [
    (3, 80),
    (2, 80),
    (4, 60),
    (1, 100),
    (5, 60),
    (6, 40)
])
def test_perc_conservacao_ignore_gap_false(perc, expected):
    alin = ['ATGA-AC',
            'AA-AT-C',
            'ATCA-AG',
            'ATGATAG',
            'ATCTCCG']
    # 1 0.8 0.4 0.8 0.4 0.6 0.6
    # ordenar
    # 0.4 0.4 0.6 0.6 0.8 0.8 1
    result = perc_conservacao(alin, perc, ignore_gap=False)
    assert result == expected


@pytest.mark.parametrize("seq, expected", [
    ("ATAGTGATCGAT", (6, [1], [2])),
    ("ABABA", (4, [0], [1])),
    ("ABCDEFG", (0, [], [])),
    ("ABCDCBA", (6, [0], [1])),
    ("ABCDFFDCBA", (9, [0], [1])),
    ("ABCDDCBA", (7, [0], [1])),
    ("ABCDCBABCD", (6, [0], [1])),
    ("ABCDCBADCBA", (6, [0], [1])),
    ("ABCDEDCBA", (8, [0], [1])),
    ("AAABBBAABBBAA", (11, [1], [2])),
    ("ABCDEFEDCBA", (10, [0], [1])),
    ("ABBBBA", (5, [0], [1])),
    ("AB", (0, [], [])),
    ("ABCDEFGIJSFEDKLMNFED", (3, [3], [10, 17]))

])
def test_maior_inversao(seq, expected):
    assert maior_inversao(seq) == expected


@pytest.mark.parametrize("proteina, esperado", [
    ('A', ['GC[ACGT]']),
    ('L', ['(CT[ACGT]|TT[AG])', '(TT[AG]|CT[ACGT])']),
    ("MSRL", ["ATG(AG[CT]|TC[ACGT])(CG[ACGT]|AG[AG])(TT[AG]|CT[ACGT])",
              "ATG(AG[CT]|TC[ACGT])(CG[ACGT]|AG[AG])(CT[ACGT]|TT[AG])",
              "ATG(AG[CT]|TC[ACGT])(AG[AG]|CG[ACGT])(TT[AG]|CT[ACGT])",
              "ATG(AG[CT]|TC[ACGT])(AG[AG]|CG[ACGT])(CT[ACGT]|TT[AG])",

              "ATG(TC[ACGT]|AG[CT])(CG[ACGT]|AG[AG])(TT[AG]|CT[ACGT])",
              "ATG(TC[ACGT]|AG[CT])(CG[ACGT]|AG[AG])(CT[ACGT]|TT[AG])",
              "ATG(TC[ACGT]|AG[CT])(AG[AG]|CG[ACGT])(TT[AG]|CT[ACGT])",
              "ATG(TC[ACGT]|AG[CT])(AG[AG]|CG[ACGT])(CT[ACGT]|TT[AG])"]),
])
def test_prot2re(proteina, esperado):
    assert prot2re(proteina) in esperado

