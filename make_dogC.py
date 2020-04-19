from model_maker import Catcher, pyrosetta
import json

pymol = pyrosetta.PyMOLMover()
dogC = Catcher(lyx=9, asx=121, glh=70, asx_type='ASN', cut_resi=105, other_res=['WAT'],
               params_folder='params',
               iso_constraint_file='constraints/iso.dogC.cst',
               trans_constraint_file='constraints/ASA-LYX.dogC.cst')
## Starting pose
print('Starting pose')
pose = dogC.load_pose_from_file('data/RrgA.altered.pdb')
pymol.pymol_name('init')
pymol.apply(pose)
dogC.relax_with_ED(pose, 'data/2ww8.ccp4')
pymol.apply(pose)
logbook = {}
s = dogC.get_score_panel(pose, save_variants=True, filename='models/00_initial')
s['description'] = 'PDB:2WW8 734-860 energy minimised against CCP4 map'
logbook['native'] = s
json.dump(logbook, open('scores.json', 'w'))
# G109T
print('G109T')
G109T = dogC.make_mutant(pose, 'G109T')
s = dogC.get_score_panel(G109T, save_variants=True, filename='models/01a_G109T')
s['description'] = 'PDB:2WW8 734-860 G109T'
logbook['G109T'] = s
json.dump(logbook, open('scores.json', 'w'))
pymol.pymol_name('G109T')
pymol.apply(G109T)
# N115G
print('N115G')
N115G = dogC.make_mutant(pose, 'N115G')
s = dogC.get_score_panel(N115G, save_variants=True, filename='models/01b_N115G')
s['description'] = 'PDB:2WW8 734-860 N115G'
logbook['N115G'] = s
json.dump(logbook, open('scores.json', 'w'))
pymol.pymol_name('N115G')
pymol.apply(N115G)
# G109T N115G
print('G109T N115G')
base = dogC.make_mutant(G109T, 'N115G')
s = dogC.get_score_panel(base, save_variants=True, filename='models/02_dogC')
s['description'] = 'PDB:2WW8 734-860 G109T N115G "DogC"'
logbook['dogC'] = s
json.dump(logbook, open('scores.json', 'w'))
pymol.pymol_name('DogC')
pymol.apply(base)
