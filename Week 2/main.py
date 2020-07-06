from src.pngs import PNGS

print('Number of PNGS: ' + str(PNGS('PPNNSNNSNPSN')))

from src.panel import Panel
relPath = './files/seap2020/136_panel_with_4lts.fa'
panel = Panel(relPath)
print('Consensus sequence:', panel.get_consensus_sequence())
panel.get_seq('HXB2.DG')
# print('Test:', panel.get_seq('HXB2.DG'))
