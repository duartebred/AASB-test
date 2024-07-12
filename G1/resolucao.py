import re


def cols_conservadas(alinhamento, perc=0.7, ignore_gap=True):
    colunas_conservadas = [] 
    for col in range(len(alinhamento[0])):
        contador = {}
        total_validos = 0

        for seq in alinhamento:
            simbolo = seq[col]
            if ignore_gap and simbolo == '-':
                continue

            total_validos += 1

            if simbolo == '-':
                continue

            if simbolo in contador:
                contador[simbolo] += 1
            else:
                contador[simbolo] = 1

        # se os simbolos forem todos traços na mesma coluna
        if total_validos == 0:
            continue

        proporcao_maxima = max(contador.values()) / total_validos

        if proporcao_maxima >= perc:
            colunas_conservadas.append(col)

    return colunas_conservadas


def perc_conservacao(alinhamento, num, ignore_gap=True):
    min_perc_per_column = []

    for col in range(len(alinhamento[0])):
        contador = {}
        total_validos = 0

        for seq in alinhamento:
            simbolo = seq[col]
            if ignore_gap and simbolo == '-':
                continue

            total_validos += 1

            if simbolo == '-':
                continue

            if simbolo in contador:
                contador[simbolo] += 1
            else:
                contador[simbolo] = 1

        if total_validos == 0:
            min_perc_per_column.append(0)
            continue

        proporcao_maxima = max(contador.values()) / total_validos
        min_perc_per_column.append(proporcao_maxima)

    min_perc_per_column.sort()
    return round(min_perc_per_column[-num] * 100, 0)


def maior_inversao(seq):
    max_len = 0
    seq_indices = []
    inv_indices = []

    # Percorre todas as sub-sequências possíveis
    for i in range(len(seq)):
        for j in range(i + 1, len(seq) + 1):
            subseq = seq[i:j]
            reversed_subseq = subseq[::-1]

            # Verifica se a sequência invertida está na sequência original
            if reversed_subseq in seq and seq.find(reversed_subseq) != i:
                # Calcula o tamanho da sequência (ou sub-sequência)
                subseq_len = len(subseq)

                # Verifica se é maior do que a maior sequência encontrada até agora
                if subseq_len > max_len:
                    max_len = subseq_len

                    seq_matches = list(re.finditer(subseq, seq))
                    inv_matches = list(re.finditer(reversed_subseq, seq))
                    seq_indices = [match.start() for match in seq_matches]
                    inv_indices = [match.start() for match in inv_matches]

    return max_len, seq_indices, inv_indices


def prot2re(proteina):
    # Tabela de códons de DNA
    DNA_table = {
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],
        'C': ['TGC', 'TGT'],
        'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'F': ['TTC', 'TTT'],
        'G': ['GGA', 'GGC', 'GGG', 'GGT'],
        'H': ['CAC', 'CAT'],
        'I': ['ATA', 'ATC', 'ATT'],
        'K': ['AAA', 'AAG'],
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCA', 'CCC', 'CCG', 'CCT'],
        'Q': ['CAA', 'CAG'],
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
        '_': ['TAA', 'TAG', 'TGA']
    }

    expr = ''
    for amino in proteina:
        codon_list = DNA_table[amino]
        if len(codon_list) == 1:
            expr += codon_list[0]
        else:
            prefixos = set(codon[:2] for codon in codon_list)
            if len(prefixos) > 1:
                expr += '('

            for i, prefixo in enumerate(prefixos):
                expr += prefixo
                codons = [codon[-1] for codon in codon_list if codon.startswith(prefixo)]
                if len(codons) == 1:
                    expr += codons[0]
                else:
                    expr += '[' + ''.join(codons) + ']'
                if i < len(prefixos) - 1:
                    expr += '|'

            if len(prefixos) > 1:
                expr += ')'

    return expr
