""" Package Imports  """
import numpy as np

""" Python File Imports"""
from src.pngs import PNGS
from src.panel import Panel
from src.blosum import BLOSUM
from src.weights import Weights
from src.epitope_dist import get_epitope_distance
from src.ic50 import IC50

""" Relative Python Paths """
rel_panel_path = './files/seap2020/136_panel_with_4lts.fa'
rel_weight_path = './files/seap2020/vrc01_wts.4lts.txt'
rel_blosum_path = './files/seap2020/BLOSUM62.txt'
rel_ic50_path = './files/seap2020/vrc01_ic50.txt'

""" Instantiating Each Class """
panel = Panel(rel_panel_path)
blosum = BLOSUM(rel_blosum_path)
weights = Weights(rel_weight_path)
weight_array_modified = np.zeros(panel.get_seq_length())
ic50 = IC50(rel_ic50_path, (panel.get_number_of_seq() - 2))

class Main:
    
    def get_consensus_sequence():
        return panel.get_consensus_sequence()

    def get_blosum_dict():
        return blosum.get_blosum_dict()

    def get_ic50_weights():
        return ic50.get_ic50_arr()

""" Turn weights -> weight_array that is usable for epitope distance """
def init_weight_array_modified():
    VRC_seq = panel.get_seq_from_name('#4lst_G(VRC01)')
    AA_order = weights.get_aa_at_an_order()
    counter = 0

    for i in range(0, panel.get_seq_length() - 1):
        # print('i', i)
        if (VRC_seq[i] != '-'):
            weight_array_modified[i] = weights.get_weight_by_order(counter)
            # print('counter', counter)
            # print('VRC_seq[i]', VRC_seq[i], 'AA_order[counter]', AA_order[counter])
            counter += 1


def epitope_distance(consensus_sequence, blosum_dict, ic50_weights):
    init_weight_array_modified()
    panel_ic50 = np.empty((2, (panel.get_number_of_seq() - 2)))  
    for i in range(1, panel.get_number_of_seq() - 1):
        panel_ic50[0][i - 1] = get_epitope_distance(panel.get_seq(i), consensus_sequence, blosum_dict, weight_array_modified)
        panel_ic50[1][i - 1] = ic50_weights[i - 1]
        # print(panel_ic50[0][i - 1], panel_ic50[1][i - 1])

    return panel_ic50

""" Testing out the Functions """

"""
print(weight_array_modified)

print('Number of PNGS: ' + str(PNGS('PPNNSNNSNPSN')))

print('Consensus sequence:', panel.get_consensus_sequence())
print('Test:', panel.get_seq('HXB2.DG'))

print('blosum.get_value("A", "R")', blosum.get_value('A', 'R'))
print('Panda BLOSUM Matrix')
print(blosum.get_panda_matrix())

print('4lst_G(VRC01) Sequence', panel.get_seq('#4lst_G(VRC01)'))

print('weights.get_weight_by_site("G.44")', weights.get_weight_by_site('G.44'))
print('weights.get_weight_by_order(53)', weights.get_weight_by_order(53))
"""
