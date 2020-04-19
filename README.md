# DogCatcher
Manuscript in preparation...

# Catcher

``model_maker.Catcher`` is class to create catcher variants as see in ``make_dogC.py`` using PyRosetta.
Class initialisation just sets the parameters.

    dogC = Catcher(lyx=9, asx=121, glh=70, asx_type='ASN', cut_resi=105, other_res=['WAT'],
                   params_folder='params',
                   iso_constraint_file='constraints/iso.dogC.cst',
                   trans_constraint_file='constraints/ASA-LYX.dogC.cst')

A pose can be loaded via 

    pose = dogC.load_pose_from_file('data/RrgA.altered.pdb')
    
which can be relaxed with a CCP4 map:

    dogC.relax_with_ED(pose, 'data/2ww8.ccp4')
    
or its isopeptide only:

    dogC.relax_isopeptide(pose)
    
or around a residue:

    dogC.relax_around_mover(pose, 123, 'A', cycles=5, distance=5, cartesian=False)
    
Then variants can be made (return a new pose):

    mutant = dogC.make_mutant('A123D')
    catcher = dogC.make_catcher_only(pose)
    tag = dogC.make_tag_only(pose)
    unreacted = dogC.make_unreacted(pose)
    trans = dogC.make_transition(pose)
    
Alternatively, get all the scores:

    s = dogC.get_score_panel(pose, save_variants=True, filename='xxxxx')
