import numpy as np
from pathlib import Path

class Panel:

    def __init__(self, filepath):
        # Importing in the file path (have to add a . because it is in /~/src instead of /~
        self.filepath = '.' +  filepath

        """ I/O the 136 HIV panel fasta sequence """
        base_path = Path(__file__).parent
        file_path = (base_path / self.filepath).resolve()
        # (x, y), where x is the sequence index and panel_dict[x] is the actual amino acid sequence
        panel_dict = {}
        # Bootstrapped coordinate compression for panel_dict sequence name and its corresponding sequence index
        # panel_compression[sequence name] == sequence index
        panel_compression = {}

        """ Inputting and reading the file """
        panel_compression_counter = 0
        with open(file_path, 'r') as f:
            prev = ''
            for x, line in enumerate(f):
                if (x % 2 == 0):
                    prev = line[1:]
                else:
                    panel_dict[panel_compression_counter] = line
                    prev = prev.split('\n')[0] # IO adds an extra \n to the end of the prev name
                    # print('prev:', prev, 'panel_compression_counter:', panel_compression_counter)
                    panel_compression[prev] = panel_compression_counter
                    panel_compression_counter += 1
        sequence_length = len(panel_dict[0])   
        self.sequence_length = sequence_length

        f.close()

        """
        # Printing out the panel_dict
        for x in panel_dict:
            print('Strand Number: ', x, ' Sequence: ', panel_dict[x])
        """

        """
        Determining the consensus sequence in O(n * k)
        where n is the number of characters in a sequence, and k is the number of sequences
        """
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

        """ Making everything a self so that other methods can access the variables  """
        self.consensus_sequence = consensus_sequence
        # (x, y), where x is the sequence index and panel_dict[x] is the actual amino acid sequence
        self.panel_dict = panel_dict
        # panel_compression[sequence name] == sequence index
        self.panel_compression = panel_compression

    """ Returning the consensus sequence  """
    def get_consensus_sequence(self):
        return self.consensus_sequence

    """ Get a HIV sequence from the panel  """
    def get_seq(self, seq_header):
        return self.panel_dict[self.panel_compression[seq_header]]

    """ Getting the length of each sequence """
    def get_seq_length(self):
        return self.sequence_length

