import sys

def is_start_codon(codon):
    return codon in ('ATG', 'GTG', 'TTG')
# include ambiguous bases
rc_n2n = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'R': 'Y',
        'Y': 'R',
        'S': 'S',
        'W': 'W',
        'K': 'M',
        'M': 'K',
        'B': 'V',
        'V': 'B',
        'D': 'H',
        'H': 'D',
        }
def rc(nucleotides):
    return ''.join(reversed([rc_n2n[n] for n in nucleotides]))

codon2aa = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAU': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACU': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGU': 'S',
    'AUA': 'I', 'AUC': 'I', 'AUG': 'M', 'AUU': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAU': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCU': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGU': 'R',
    'CUA': 'L', 'CUC': 'L', 'CUG': 'L', 'CUU': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAU': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCU': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGU': 'G',
    'GUA': 'V', 'GUC': 'V', 'GUG': 'V', 'GUU': 'V',
    'UAA': '*', 'UAC': 'Y', 'UAG': '*', 'UAU': 'Y',
    'UCA': 'S', 'UCC': 'S', 'UCG': 'S', 'UCU': 'S',
    'UGA': 'W', 'UGC': 'C', 'UGG': 'W', 'UGU': 'C',
    'UUA': 'L', 'UUC': 'F', 'UUG': 'L', 'UUU': 'F',
}


# There is a particular file from ProGenomes that contained a null character
# (\0) in the middle of the sequence. This confused the `prodigal_sm` algorithm
# and it returned a mistranslated sequence.
# The sequence is: GMSC10.100AA.409_088_494
# Thus, we have to hardcode the translation for this particular sequence.
SPECIAL_CASE_NUCLEOTIDES = 'ATGGCTAAGGGGCAATCTTTACAAGATCCGTTCCTGAACGCATTGCGTCGGGAACGTGTTCCAGTTTCTATTTATTTGGTGAATGGTATTAAGCTGCAAGGTCAAATCGAGTCCTTTGATCAGTTCGTG'
SPECIAL_CASE_AA = 'MRRLTREHCDELISIPMAGSVSSLNVSVATGICLFEAVRQRT'

for k, v in list(codon2aa.items()):
    if 'U' in k:
        codon2aa[k.replace('U', 'T')] = v
def translate(nucleotides):
    """Translate nucleotides to amino acids."""
    if nucleotides == SPECIAL_CASE_NUCLEOTIDES:
        return SPECIAL_CASE_AA
    if not is_start_codon(nucleotides[:3]):
        nucleotides = rc(nucleotides)
        if not is_start_codon(nucleotides[:3]):
            raise ValueError(f"Invalid start codon: {nucleotides[:3]}")

    rs = [codon2aa.get(nucleotides[i:i+3], 'X') for i in range(0, len(nucleotides), 3)]
    if rs[-1] != '*' and rs[-1] != 'W':
        raise ValueError(f"Unterminated sequence: {nucleotides}")
    if rs[0] != 'M':
        rs[0] = 'M'
    return ''.join(rs[:-1])

