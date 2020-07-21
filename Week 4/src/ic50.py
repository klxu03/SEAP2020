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
        ic50_arr = np.empty((number_of_sequences, 2))
        # ic50_panel_compression[0] is going to output the sequence name of the 1st IO amino acid
        ic50_panel_compression = np.empty((number_of_sequences), dtype="object")

        with open(file_path, 'r') as f:
            for x, line in enumerate(f):
                # Input into the panel_compression
                ic50_arr[x][0] = x
                ic50_panel_compression[x] = (line.split())[0]

                # Inputing into the ic50_arr
                ic50_arr[x][1] = (line.split())[1]

        self.ic50_panel_compression = ic50_panel_compression
        self.ic50_arr = ic50_arr

    def get_ic50_arr(self):
        return self.ic50_arr[:,1]

    def get_lowest_ic50_sequences(self, amount):
        ic50_arr_sort = self.ic50_arr[self.ic50_arr[:,1].argsort()]
        sequences = np.empty(amount, dtype="object")
        for i in range(amount):
            # print('ic50_arr_sort[i][0]', int(ic50_arr_sort[i][0]))
            sequences[i] = self.ic50_panel_compression[int(ic50_arr_sort[i][0])]

        # print('self.ic50_panel_compression[7]', self.ic50_panel_compression[7])
        # print(self.ic50_panel_compression[int(ic50_arr_sort[0][0])])
        # print('ic50_arr_sort', ic50_arr_sort)
        return sequences
