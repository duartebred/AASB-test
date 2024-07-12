import math

import pytest

from resolucao import Tree


@pytest.fixture
def tree_example():
    s1 = Tree('ACGT')
    s2 = Tree('ACT')
    s3 = Tree('ACG')
    b1 = Tree(s1, s3)
    return Tree(b1, s2)


@pytest.fixture
def tree_example2():
    s1 = Tree('ACGT')
    s2 = Tree('ACT')
    s3 = Tree('ACG')
    s4 = Tree('ACGTA')
    s5 = Tree('AAAA')
    b1 = Tree(s1, s3)
    b2 = Tree(s2, s4)
    b3 = Tree(s5, b1)
    return Tree(b3, b2)


@pytest.mark.parametrize("seq1, seq2, expected_dist", [
    ('ACG', 'ACGT', 1),
    ('ACT', 'ACG', 2),
    ('ACGT', 'ACGT', 0),
    ('ACG', 'AAAA', math.inf),
    ('ACG', 'ACGTA', math.inf),
])
def test_dist(tree_example, seq1, seq2, expected_dist):
    assert tree_example.dist(seq1, seq2) == expected_dist


@pytest.mark.parametrize("seq1, seq2, expected_dist", [
    ('ACG', 'ACGT', 1),
    ('ACT', 'ACG', 4),
    ('ACGT', 'ACGT', 0),
    ('ACG', 'AAAA', 2),
    ('ACG', 'ACGTA', 4),
])
def test_dist(tree_example2, seq1, seq2, expected_dist):
    print(tree_example2)
    assert tree_example2.dist(seq1, seq2) == expected_dist


@pytest.mark.parametrize("seq, expected_depth", [
    ('ACGT', 2),
    ('ACT', 1),
    ('AAAA', math.inf),
    ('ACG', 2),
    ('ACGTA', math.inf),
])
def test_depth(tree_example, seq, expected_depth):
    assert tree_example.depth(seq) == expected_depth


