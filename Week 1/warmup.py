""" Warmup Tasks  """
import numpy as np
from pathlib import Path

""" 
A function to get the potential N-linked glycolysation sites (PNGS) for a given Amino Acid sequence
https://www.hiv.lanl.gov/content/sequence/GLYCOSITE/glycosite.html
"""
def PNGS(seq):
    counter = 0
    for x in range(len(seq) - 3):
        if (seq[x] == 'N' and seq[x+1] != 'P' and (seq[x+2] == 'S' or seq[x+2] == 'T') and seq[x+3] != 'P'):
            counter += 1
    return counter

# Should return 2
print('Number of PNGS: ' + str(PNGS('PPNNSNNSNPSN')))

""" I/O the 136 HIV panel fasta sequence  """
base_path = Path(__file__).parent
file_path = (base_path / './files/seap2020/136_panel_with_4lts.fa').resolve()
# file_path = (base_path / './files/seap2020/test.fa').resolve()
# (x, y), where x is the sequence index and panel_dict[x] is the actual amino acid sequence 
panel_dict = {} 
# Bootstrapped coordinate compression for panel_dict sequence name and its corresponding sequence index
# panel_compression[sequence name] == sequence index
panel_compression = {}

panel_compression_counter = 0
with open(file_path, 'r') as f:
    prev = ''
    for x, line in enumerate(f):
        if (x % 2 == 0):
            prev = line[1:]
        else:
            panel_dict[panel_compression_counter] = line
            panel_compression[line] = panel_compression_counter
            panel_compression_counter += 1
sequence_length = len(panel_dict[0])   

f.close()

"""
# Printing out the panel_dict
for x in panel_dict:
    print('Strand Name: ', x, ' Sequence: ', panel_dict[x])
"""

"""
Determining the consensus sequence in O(n * k)
where n is the number of characters in a sequence, and k is the number of sequences
"""
# print('Each sequence\'s length: ', sequence_length)
consensus_sequence = ''
for x in range(sequence_length - 1):
    letter_counter = np.zeros((27), dtype=int)
    for y in panel_dict:
        if (panel_dict[y][x] != '-'):
            letter_counter[ord(panel_dict[y][x]) - 65] += 1
        else:
            letter_counter[26] += 1
    # best_letter[0][0] gives the first index of the indices that have the highest value
    # best_letter[0] gives all the indices where that element's value is the same as the max
    # For example, if array is [1, 5, 5], best_letter[0][0] gives 1, whereas best_letter[0] gives [1 2]
    best_letter = np.where(letter_counter == np.amax(letter_counter))
    if (best_letter[0][0] != 26):
        consensus_sequence += chr(best_letter[0][0] + 65)
    else:
        consensus_sequence += '-'

print('Consensus Sequence:', consensus_sequence)
