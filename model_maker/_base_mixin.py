import pyrosetta

pyrosetta.init(extra_options='-no_optH false -mute all -load_PDB_components false')

import os, re

from typing import Optional, List, Union


class _Base:

    def __init__(self, lyx: int, asx: int, glh: int,
                 asx_type: str,  # ASP || ASN?
                 cut_resi: int,
                 other_res: Optional[List[Union[int, str]]] = None,
                 params_folder: str = 'params',
                 iso_constraint_file: Optional[str] = None,
                 trans_constraint_file: Optional[str] = None,
                 chain: str = 'A'):
        self.params_folder = params_folder
        self.lyx = lyx
        self.asx = asx
        self.glh = glh
        self.asx_type = asx_type
        self.cut_resi = cut_resi
        self.iso_constraint_file = iso_constraint_file
        self.trans_constraint_file = trans_constraint_file
        if other_res is None:
            self.other_res = []
        else:
            self.other_res = other_res
        self.chain = chain

    def load_pose_from_file(self, filename: str, constraint_file: Optional[str] = None) -> pyrosetta.Pose:
        """
        Loads a pose from filename with the params in the params_folder

        :param filename:
        :param constraint_file:
        :return:
        """
        pose = pyrosetta.Pose()
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([os.path.join(self.params_folder, file) for file in os.listdir(self.params_folder)
                             if os.path.splitext(file)[1] == '.params'])
        pyrosetta.generate_nonstandard_residue_set(pose, params_paths)
        pyrosetta.rosetta.core.import_pose.pose_from_file(pose, filename)
        if constraint_file is None and self.iso_constraint_file is not None:
            constraint_file = self.iso_constraint_file
        if constraint_file:
            setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
            setup.constraint_file(constraint_file)
            setup.apply(pose)
        return pose

    def _print_resi_in_vector(self, vector):
        ResidueVector = pyrosetta.rosetta.core.select.residue_selector.ResidueVector
        print(ResidueVector(vector))

    def _print_scores_per_res(self, pose):
        data = pose.energies().residue_total_energies_array()
        for i in range(data.size):
            print(pose.residue(i + 1).name1(), i + 1, data[i]['total_score'])


    def get_effective_pKa(self, pose):
        raise Exception('This causes a seg fault.')
        dummy = pose.clone()
        prot_rosetta_db = '/Users/matteo/miniconda3/lib/python3.7/site-packages/pyrosetta/database/chemical/residue_type_sets/fa_standard/residue_types/protonation_states/'
        lys_pos = dummy.pdb_info().pdb2pose(chain='A', res=self.lyx)
        asn_pos = dummy.pdb_info().pdb2pose(chain='A', res=self.asx)
        glu_pos = dummy.pdb_info().pdb2pose(chain='A', res=self.glh)
        # RESCON: 305 LIG n-conn= 1 n-poly= 0 n-nonpoly= 1 conn# 1 22 145 3
        dummy.conformation().sever_chemical_bond(seqpos1=lys_pos, res1_resconn_index=3, seqpos2=asn_pos,
                                                 res2_resconn_index=3)
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=lys_pos, new_res='ALA').apply(dummy)
        MutateResidue(target=asn_pos, new_res='ALA').apply(dummy)
        MutateResidue(target=glu_pos, new_res='ALA').apply(dummy)
        dummy.delete_residue_slow(pose.total_residue())
        params_paths = pyrosetta.rosetta.utility.vector1_string()
        params_paths.extend([os.path.join(prot_rosetta_db, resn + '.params') for resn in ('LYS_D',)])
        pyrosetta.generate_nonstandard_residue_set(dummy, params_paths)
        return pyrosetta.rosetta.protocols.simple_moves.ReportEffectivePKA().apply(dummy)
