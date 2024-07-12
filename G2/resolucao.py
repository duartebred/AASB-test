import math


class Tree:
    def __init__(self, *args):
        def is_seq(seq):
            return type(seq) is str and all(b in "ACGT" for b in seq)

        if len(args) == 1:
            assert is_seq(args[0])
            self.leaf = args[0]
        else:
            assert len(args) == 2
            assert all(type(s) is Tree for s in args)
            self.left, self.right = args

    def dist(self, seq1, seq2):
        """
        Devolve a distância na árvore entre seq1 e seq2 ou infinito
        caso ambas as sequências não se encontrem na árvore
        """
        if not hasattr(self, 'leaf'):
            if seq1 == seq2:
                return 0

            left_dist = self.left.dist(seq1, seq2)
            right_dist = self.right.dist(seq1, seq2)
            cross_dist = self.left.depth(seq1) + self.right.depth(seq2) + 1 if self.left.contains(
                seq1) and self.right.contains(seq2) else math.inf
            cross_dist = min(cross_dist, self.right.depth(seq1) + self.left.depth(seq2) + 1 if self.right.contains(
                seq1) and self.left.contains(seq2) else math.inf)
            return min(left_dist, right_dist, cross_dist)
        return math.inf

    def depth(self, seq):
        """
        Devolve a profundidade na árvore onde está a sequência seq
        ou infinito caso não esteja na árvore
        """
        if hasattr(self, 'leaf'):
            return 0 if self.leaf == seq else math.inf
        else:
            return 1 + min(self.left.depth(seq), self.right.depth(seq))

    def contains(self, seq):
        """Verifica se a árvore contém a sequência seq"""
        if hasattr(self, 'leaf'):
            return self.leaf == seq
        return self.left.contains(seq) or self.right.contains(seq)

    def __str__(self, indent=''):
        if hasattr(self, 'leaf'):
            return indent + self.leaf
        return f"""{indent}[
{indent} {self.left.__str__(indent=indent + ' ')}
{indent} {self.right.__str__(indent=indent + ' ')}
{indent}]"""

    def __repr__(self):
        return str(self)


