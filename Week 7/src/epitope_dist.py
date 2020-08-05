"""
The function that returns the epitope distance between two strands
"""

# The blosum parameter is a boolean value that is true if I am using the blosum matrix_dict, and false if I am not using the blosum matrix
def get_epitope_distance(seq, ref_seq, matrix_dict, weights, blosum):
    num = 0
    denom = 0
    for i in range(0, len(seq) - 1):
        seq_letter = seq[i]
        ref_letter = ref_seq[i]

        if (blosum):
            if (seq_letter == '-'):
                seq_letter = '*'
            if (ref_letter == '-'):
                ref_letter = '*'
            
            num += weights[i] * matrix_dict[seq_letter, ref_letter]
        else:
            min_val_in_pmbec = -0.13
            if (seq_letter == '-' or ref_letter == '-'):
                num += weights[i] * min_val_in_pmbec 
            else:
                num += weights[i] * matrix_dict[seq_letter, ref_letter]

        denom += weights[i]
        # print('i', i, 'new num', num, 'new denom', denom)
        # print('seq[i]', seq[i], 'ref_seq[i]', ref_seq[i], 'weights[i]', weights[i])

    return num/denom
