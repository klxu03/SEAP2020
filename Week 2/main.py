""" Package Imports  """
import numpy as np
# from pathlib import Path

""" Python File Imports"""
from src.pngs import PNGS
from src.panel import Panel
from src.blosum import BLOSUM

""" Relative Python Paths """
rel_panel_path = './files/seap2020/136_panel_with_4lts.fa'
rel_weight_path = './files/seap2020/vrc01_wts.4lts.txt'
rel_blosum_path = './files/seap2020/BLOSUM62.txt'

print('Number of PNGS: ' + str(PNGS('PPNNSNNSNPSN')))

panel = Panel(rel_panel_path)
print('Consensus sequence:', panel.get_consensus_sequence())
print('Test:', panel.get_seq('HXB2.DG'))

blosum = BLOSUM(rel_blosum_path)
print('blosum.get_value("A", "R")', blosum.get_value('A', 'R'))
print('Panda BLOSUM Matrix')
print(blosum.get_panda_matrix())
