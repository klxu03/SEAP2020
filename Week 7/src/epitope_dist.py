"""
The function that returns the epitope distance between two strands
"""

def get_epitope_distance(seq, ref_seq, matrix_dict, weights):
    num = 0
    denom = 0
    for i in range(0, len(seq) - 1):
        seq_letter = seq[i]
        ref_letter = ref_seq[i]
        if (seq_letter == '-'):
            seq_letter = '*'
        if (ref_letter == '-'):
            ref_letter = '*'

        print('seq_letter', seq_letter)
        print('ref_letter', ref_letter)
        num += weights[i] * (matrix_dict[seq_letter, ref_letter])
        denom += weights[i]
        print('i', i, 'new num', num, 'new denom', denom)
        # print('seq[i]', seq[i], 'ref_seq[i]', ref_seq[i], 'weights[i]', weights[i])

    return num/denom
