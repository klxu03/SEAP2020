import numpy as np
from pathlib import Path
import pandas as pd

class PMBEC:

    def __init__(self, filepath):
        # Importing in the file path (have to add a . because it is in /~/src instead of /~
        self.filepath = '.' +  filepath

        """ I/O the 136 HIV panel fasta sequence """
        base_path = Path(__file__).parent
        file_path = (base_path / self.filepath).resolve()

        letters = "A C D E F G H I K L M N P Q R S T V W Y".split()
        # letter_compression[index] == letter
        # letter_compression[0] == 'A'
        letter_compression = {} 
        for num, letter in enumerate(letters, start=0):
            letter_compression[num] = letter        

        # Matrix is the Panda Dataframe Showing Blosum
        numpy_matrix = np.zeros((24, 24))
        
        # blosum_dict[x_letter, y_letter] == Value
        # For example, blosum_dict['A', 'R'] will be -1 
        blosum_dict = {}

        with open(file_path, 'r') as f:
            for x, line in enumerate(f):
                input = line.split()
                # print('input', input)
                y = 0
                if x != 0 and x < 25:
                    x_letter = letter_compression[x - 1]
                    for character in input:
                        if y != 0:
                            y_letter = letter_compression[y - 1]
                            el = float(character)
                            numpy_matrix[x - 1][y - 1] = el
                            blosum_dict[x_letter, y_letter] = el
                        y += 1

        """ panda_matrix = pd.DataFrame(data=numpy_matrix[0:, 0:],
                index = [i for i in letters],
                columns = [i for i in letters]) """

        # print(panda_matrix)
        # Setting everything to self so that other functions can access these variables
        self.blosum_dict = blosum_dict
        # self.panda_matrix = panda_matrix

    def get_value(self, letter1, letter2):
        return self.blosum_dict[letter1, letter2]

    def get_panda_matrix(self):
        return self.panda_matrix

    def get_blosum_dict(self):
        return self.blosum_dict
