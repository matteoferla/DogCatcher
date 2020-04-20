from model_maker import Catcher, pyrosetta
import json

pymol = pyrosetta.PyMOLMover()
dogC = Catcher(lyx=9, asx=121, glh=70, asx_type='ASN', cut_resi=105, other_res=['WAT'],
               params_folder='params',
               iso_constraint_file='constraints/iso.dogC.cst',
               trans_constraint_file='constraints/ASA-LYX.dogC.cst')

## Starting pose
print('Starting pose')
#pose = dogC.load_pose_from_file('data/RrgA.altered.pdb')
pose = dogC.load_pose_from_file('../RrgA.relaxed.pdb')
pymol.pymol_name('init')
pymol.apply(pose)
# dogC.relax_with_ED(pose, 'data/2ww8.ccp4')
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
#A75P
print('A75P')
A75P = dogC.make_mutant(base, 'A75P')
dogC.relax_loop(A75P, 73, 80)
s = dogC.get_score_panel(A75P, save_variants=True, filename='models/03_A75P')
s['description'] = 'A75P'
logbook['A75P'] = s
json.dump(logbook, open('scores.json', 'w'))

pair_A = dogC.make_double_mutant(A75P, ['N11D', 'N13T'])
s = dogC.get_score_panel(pair_A, save_variants=True, filename='models/04a_N11D_N13T')
s['description'] = 'N11D N13T A75P'
logbook['N11D N13T'] = s

pair_B = dogC.make_double_mutant(A75P, ['D4E', 'K59T'])
s = dogC.get_score_panel(pair_B, save_variants=True, filename='models/04b_D4E_K59T')
s['description'] = 'D4E K59T A75P'
logbook['D4E K59T'] = s

pair_C = dogC.make_double_mutant(A75P, ['A87E', 'I101A'])
s = dogC.get_score_panel(pair_C, save_variants=True, filename='models/04c_A87E_I101A')
s['description'] = 'A75P A87E I101A'
logbook['A87E I101A'] = s

quad = dogC.make_double_mutant(pair_A, ['D4E', 'K59T'])
s = dogC.get_score_panel(quad, save_variants=True, filename='models/05_D4E_N11D_N13T_K59T')
s['description'] = 'D4E N11D N13T K59T'
logbook['D4E N11D N13T K59T A75P'] = s

for letter, resi in (('d', 'A38P'), ('e','Y45G'), ('f','N47D'), ('g','N92D'), ('h','A87E')):
    x = dogC.make_mutant(A75P, resi)
    s = dogC.get_score_panel(x, save_variants=True, filename=f'models/04{letter}_{resi}')
    s['description'] = f'A75P {resi}'
    logbook[resi] = s
    json.dump(logbook, open('scores.json', 'w'))

pair_D = dogC.make_double_mutant(A75P, ['A87E', 'I101A'])
s = dogC.get_score_panel(pair_D, save_variants=True, filename='models/04i_N47D_N92D')
s['description'] = 'N47D A75P N92D'
logbook['N47D A75P N92D'] = s

aqua = dogC.make_double_mutant(quad, ['N92D', 'N47D'])
s = dogC.get_score_panel(aqua, save_variants=True, filename='models/06_N47D_N92D')
s['description'] = 'D4E N11D N13T N47D A75P K59T N92D'
logbook['D4E N11D N13T N47D A75P K59T N92D'] = s

F69I = dogC.make_mutant(aqua, 'F69I')
s = dogC.get_score_panel(F69I, save_variants=True, filename='models/07a_F69I')
s['description'] = '+ F69I'
logbook['F69I'] = s
json.dump(logbook, open('scores.json', 'w'))

Q89R = dogC.make_mutant(aqua, 'Q89R')
s = dogC.get_score_panel(Q89R, save_variants=True, filename='models/07b_Q89R')
s['description'] = '+ Q89R'
logbook['Q89R'] = s
json.dump(logbook, open('scores.json', 'w'))

A87S = dogC.make_mutant(aqua, 'A87S')
s = dogC.get_score_panel(A87S, save_variants=True, filename='models/07c_A87S')
s['description'] = '+ A87S'
logbook['A87S'] = s
json.dump(logbook, open('scores.json', 'w'))

phage = dogC.make_double_mutant(aqua, ['Q89R', 'A87S', 'F69I'])
s = dogC.get_score_panel(phage, save_variants=True, filename='models/08_F69I_A87S_Q89R')
s['description'] = '+ F69I A87S Q89R'
logbook['F69I A87S Q89R'] = s
json.dump(logbook, open('scores.json', 'w'))