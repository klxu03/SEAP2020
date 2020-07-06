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
