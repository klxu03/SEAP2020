import numpy as np
from pathlib import Path

class IC50:

    def __init__(self, filepath, number_of_sequences):
        # Importing in the file path (have to add a . because it is in /~/src instead of /~
        self.filepath = '.' +  filepath

        """ I/O the IC50 concentration """
        base_path = Path(__file__).parent
        file_path = (base_path / self.filepath).resolve()

        """ ic50_arr[0] is the weight of the 0th index (1st) strand """
        ic50_arr = np.empty(number_of_sequences)

        with open(file_path, 'r') as f:
            for x, line in enumerate(f):
                # print('(line.split())[1]', (line.split())[1])
                ic50_arr[x] = (line.split())[1]

        self.ic50_arr = ic50_arr

    def get_ic50_arr(self):
        return self.ic50_arr
