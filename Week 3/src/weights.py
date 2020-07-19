from pathlib import Path

class Weights:

    def __init__(self, filepath):
        # Importing in the file path (have to add a . because it is in /~/src instead of /~
        self.filepath = '.' +  filepath

        """ I/O the weights file """
        base_path = Path(__file__).parent
        file_path = (base_path / self.filepath).resolve()

        """ site_dict['G.44'] = 0.000 (weight)  """
        site_dict = {}

        """ site_dict_panel_compression[0 (0th in IO)] = 'G.44' """
        site_dict_panel_compression = {}

        """ AA_at_an_order[0 (0th in IO)] = 'V' (1st Amino Acid) """
        AA_at_an_order = {}

        """ Inputting and reading the file """
        AA_counter = 0
        with open(file_path, 'r') as f:
            for x, line in enumerate(f):
                input = line.split()
                # print('Site', input[0], 'Amino Acid', input[1], 'Weight', input[2])
                site_dict[input[0]] = input[2]
                site_dict_panel_compression[AA_counter] = input[0]
                AA_at_an_order[AA_counter] = input[1]
                AA_counter += 1

        # print('AA_counter', AA_counter)
        self.site_dict = site_dict
        self.site_dict_panel_compression = site_dict_panel_compression
        self.AA_at_an_order = AA_at_an_order

    """ Use case: print('weights.get_weight_by_site("G.44")', weights.get_weight_by_site('G.44')) """
    def get_weight_by_site(self, site):
        return self.site_dict[site]

    """ Use case: print('weights.get_weight_by_order(53)', weights.get_weight_by_order(53)) """
    def get_weight_by_order(self, number):
        return self.site_dict[self.site_dict_panel_compression[number]]

    def get_aa_at_an_order(self):
        return self.AA_at_an_order
