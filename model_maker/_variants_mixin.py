import re

from typing import Optional, List, Union

from ._relax_mixin import _Relax, pyrosetta


class _Make(_Relax):

    def make_mutant(self, pose: pyrosetta.Pose, mutation, chain='A',
                    constraint_file: Optional[str] = None) -> pyrosetta.Pose:
        """
        Make a point mutant (``A23D``).

        :param pose: pose
        :param mutation:
        :param chain:
        :param constraint_file:
        :return:
        """
        mutant = pose.clone()
        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_file(constraint_file)
        setup.apply(mutant)
        name3 = {'A': 'ALA',
                 'C': 'CYS',
                 'D': 'ASP',
                 'E': 'GLU',
                 'F': 'PHE',
                 'G': 'GLY',
                 'H': 'HIS',
                 'I': 'ILE',
                 'L': 'LEU',
                 'K': 'LYS',
                 'M': 'MET',
                 'N': 'ASN',
                 'P': 'PRO',
                 'Q': 'GLN',
                 'R': 'ARG',
                 'S': 'SER',
                 'T': 'THR',
                 'V': 'VAL',
                 'W': 'TRP',
                 'Y': 'TYR'}
        pose2pdb = pose.pdb_info().pdb2pose
        rex = re.match('(\w)(\d+)(\w)', mutation)
        r = pose2pdb(res=int(rex.group(2)), chain=chain)
        rn = pose.residue(r).name1()
        assert rn == rex.group(1), f'residue {r}(pose)/{rex.group(2)}(pdb) is a {rn}, not a {rex.group()}'
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=r, new_res=name3[rex.group(3)]).apply(mutant)
        self.relax_around_mover(mutant, int(rex.group(2)), 'A', distance=7, cycles=15)
        return mutant

    def make_quick_unreacted(self, pose) -> pyrosetta.Pose:
        split_pose = pyrosetta.Pose()
        split_pose.assign(pose)
        lys_pos = split_pose.pdb_info().pdb2pose(chain=self.chain, res=self.lyx)
        asn_pos = split_pose.pdb_info().pdb2pose(chain=self.chain, res=self.asx)
        # RESCON: 305 LIG n-conn= 1 n-poly= 0 n-nonpoly= 1 conn# 1 22 145 3
        split_pose.conformation().sever_chemical_bond(seqpos1=lys_pos,
                                                      res1_resconn_index=3,
                                                      seqpos2=asn_pos,
                                                      res2_resconn_index=3)
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=lys_pos, new_res='ALA').apply(split_pose)
        MutateResidue(target=lys_pos, new_res='LYS').apply(split_pose)
        MutateResidue(target=asn_pos, new_res='ALA').apply(split_pose)
        MutateResidue(target=asn_pos, new_res=self.asx_type).apply(split_pose)
        self.relax_isopeptide(split_pose, distance=5, cycles=3)
        return split_pose

    def make_unreacted(self, pose) -> pyrosetta.Pose:
        split_pose = self.make_quick_unreacted(pose)
        setup = pyrosetta.rosetta.protocols.constraint_movers.ClearConstraintsMover()
        setup.apply(split_pose)
        self.relax_isopeptide(split_pose, cartesian=True, distance=1, cycles=15)
        self.relax_isopeptide(split_pose)
        return split_pose

    def make_unbound(self, pose) -> pyrosetta.Pose:
        print('This fails to break the bond.')
        split_pose = self.make_quick_unreacted(pose)
        get_resn = lambda r: split_pose.residue(r).name3()
        last_pos = split_pose.pdb_info().pdb2pose(chain='A', res=self.cut_resi - 1)
        tag_pos = split_pose.pdb_info().pdb2pose(chain='A', res=self.cut_resi)
        # RESCON: 305 LIG n-conn= 1 n-poly= 0 n-nonpoly= 1 conn# 1 22 145 3
        split_pose.conformation().sever_chemical_bond(seqpos1=last_pos,
                                                      res1_resconn_index=2,
                                                      seqpos2=tag_pos,
                                                      res2_resconn_index=1)
        ## Move away
        xyz = pyrosetta.rosetta.numeric.xyzVector_double_t()
        xyz.x = 500.0
        xyz.y = 0.0
        xyz.z = 0.0
        for r in range(tag_pos, split_pose.total_residue()):  # what about that water?
            for a in range(1, split_pose.residue(r).natoms() + 1):
                split_pose.residue(r).set_xyz(a, split_pose.residue(r).xyz(a) + xyz)
        ## termini
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        MutateResidue(target=last_pos, new_res=f'{get_resn(last_pos)}:CtermProteinFull').apply(split_pose)
        MutateResidue(target=tag_pos, new_res=f'{get_resn(tag_pos)}:NtermProteinFull').apply(split_pose)
        ft = pyrosetta.FoldTree()
        ft.add_edge(1, self.cut_resi - 1, -1)  # chain
        ft.add_edge(self.cut_resi - 1, self.cut_resi, 1)  # jump
        ft.add_edge(self.cut_resi, self.total_residue(), -1)  # chain
        ft.add_edge(127, 128, 2)  # jump
        # ft.add_edge(127, 128, -1) # water
        print(str(ft), ft.check_fold_tree())
        split_pose.fold_tree(ft)
        relax = self.relax_iso(split_pose)
        relax.apply(split_pose)
        return split_pose

    def make_catcher_only(self, pose):
        catcher = self.make_quick_unreacted(pose)
        setup = pyrosetta.rosetta.protocols.constraint_movers.ClearConstraintsMover()
        setup.apply(catcher)
        pyrosetta.rosetta.protocols.grafting.delete_region(catcher, self.cut_resi, pose.total_residue())
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 5)
        relax.apply(catcher)
        return catcher

    def make_tag_only(self, pose) -> pyrosetta.Pose:
        tag = self.make_quick_unreacted(pose)
        pyrosetta.rosetta.protocols.grafting.delete_region(tag, 1, self.cut_resi - 1)
        setup = pyrosetta.rosetta.protocols.constraint_movers.ClearConstraintsMover()
        setup.apply(tag)
        scorefxn = pyrosetta.get_fa_scorefxn()
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 15)
        relax.apply(tag)
        return tag

    def make_transition(self, pose, constraint_file: Optional[str] = None) -> pyrosetta.Pose:
        trans = pose.clone()
        MutateResidue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue
        lys_pos = trans.pdb_info().pdb2pose(chain='A', res=self.lyx)
        asn_pos = trans.pdb_info().pdb2pose(chain='A', res=self.asx)
        MutateResidue(target=lys_pos, new_res=f'LYX').apply(trans)
        if self.asx_type == 'ASN':
            MutateResidue(target=asn_pos, new_res=f'ASA').apply(trans)
        elif self.asx_type == 'ASP':
            MutateResidue(target=asn_pos, new_res=f'ASL').apply(trans)
        else:
            raise ValueError(f'What is a {self.asx_type}')
        setup = pyrosetta.rosetta.protocols.constraint_movers.ClearConstraintsMover()
        setup.apply(trans)
        if constraint_file:
            constraint_file = self.trans_constraint_file
        setup = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
        setup.constraint_file(constraint_file)
        setup.apply(trans)
        self.relax_isopeptide(pose, cartesian=True, distance=1, cycles=15)
        self.relax_isopeptide(pose, distance=7, cycles=3)
        return trans
