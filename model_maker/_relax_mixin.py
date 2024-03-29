from ._base_mixin import pyrosetta, _Base

from typing import Optional, List, Union


class _Relax(_Base):

    def _cst_score(self, scorefxn):
        stm = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
        for contype_name in ("atom_pair_constraint", "angle_constraint", "dihedral_constraint"):
            contype = stm.score_type_from_name(contype_name)
            scorefxn.set_weight(contype, 5)

    def relax_with_ED(self, pose, ccp4_file: str) -> None:
        """
        Relaxes ``pose`` based on the ccp4 electron density map provided.

        :param pose:
        :param ccp4_file: download map from ePDB
        :return: Relaxes pose in place
        """
        scorefxnED = pyrosetta.get_fa_scorefxn()
        ED = pyrosetta.rosetta.core.scoring.electron_density.getDensityMap(ccp4_file)
        sdsm = pyrosetta.rosetta.protocols.electron_density.SetupForDensityScoringMover()
        sdsm.apply(pose)
        ## Set ED constraint
        elec_dens_fast = pyrosetta.rosetta.core.scoring.ScoreType.elec_dens_fast
        scorefxnED.set_weight(elec_dens_fast, 30)
        ## Set generic constraints
        self._cst_score(scorefxnED)
        ## Relax
        for w in (30, 20, 10):
            scorefxnED.set_weight(elec_dens_fast, w)
            relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxnED, 5)
            relax.apply(pose)

    def relax_around_mover(self,
                           pose: pyrosetta.Pose,
                           resi: int, chain: str,
                           scorefxn=None, cycles=5, distance=5, cartesian=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain.

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param scorefxn:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :param cartesian:
        :return:
        """
        if scorefxn is None:
            scorefxn = pyrosetta.get_fa_scorefxn()
            self._cst_score(scorefxn)
        movemap = pyrosetta.MoveMap()
        ####
        resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(resi_sele, distance=distance, include_focus_in_subset=True)
        n = neigh_sele.apply(pose)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def relax_pair_mover(self,
                           pose: pyrosetta.Pose,
                           resi_pair: List[int], chain: str,
                           scorefxn=None, cycles=5, distance=5, cartesian=False) -> None:
        """
        Relaxes pose ``distance`` around resi:chain.

        :param resi: PDB residue number.
        :param chain:
        :param pose:
        :param scorefxn:
        :param cycles: of relax (3 quick, 15 thorough)
        :param distance:
        :param cartesian:
        :return:
        """
        if scorefxn is None:
            scorefxn = pyrosetta.get_fa_scorefxn()
            self._cst_score(scorefxn)
        movemap = pyrosetta.MoveMap()
        ####
        or_sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
        for resi in resi_pair:
            resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
            resi_sele.set_index(pose.pdb_info().pdb2pose(chain=chain, res=resi))
            or_sele.add_residue_selector(resi_sele)
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        neigh_sele = NeighborhoodResidueSelector(or_sele, distance=distance, include_focus_in_subset=True)
        n = neigh_sele.apply(pose)
        movemap.set_bb(allow_bb=n)
        movemap.set_chi(allow_chi=n)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def relax_isopeptide(self, pose, cartesian=False, distance=7, cycles=5) -> None:
        or_sele = pyrosetta.rosetta.core.select.residue_selector.OrResidueSelector()
        for res in [self.lyx, self.asx, self.glh, *self.other_res]:
            if isinstance(res, int):
                resi_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
                resi_sele.set_index(pose.pdb_info().pdb2pose(chain=self.chain, res=res))
                or_sele.add_residue_selector(resi_sele)
            elif isinstance(res, str):
                resn_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueNameSelector()
                resn_sele.set_residue_name3(res)
                or_sele.add_residue_selector(resn_sele)
            else:
                raise ValueError
        NeighborhoodResidueSelector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector
        ns = NeighborhoodResidueSelector(or_sele, distance=distance, include_focus_in_subset=True)
        nv = ns.apply(pose)
        movemap = pyrosetta.MoveMap()
        movemap.set_bb(allow_bb=nv)
        movemap.set_chi(allow_chi=nv)
        movemap.set_jump(1, True)
        if cartesian:
            scorefxn = pyrosetta.create_score_function('ref2015_cart')
        else:
            scorefxn = pyrosetta.create_score_function('ref2015')
        self._cst_score(scorefxn)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(cartesian)
        relax.apply(pose)

    def relax_loop(self, pose, start, stop):
        span_sele = pyrosetta.rosetta.core.select.residue_selector.ResidueSpanSelector(start, stop)
        movemap = pyrosetta.MoveMap()
        sv = span_sele.apply(pose)
        movemap.set_bb(allow_bb=sv)
        movemap.set_chi(allow_chi=sv)
        movemap.set_jump(1, True)
        scorefxn = pyrosetta.create_score_function('ref2015_cart')
        self._cst_score(scorefxn)
        omega = pyrosetta.rosetta.core.scoring.ScoreType.omega
        fa_rep = pyrosetta.rosetta.core.scoring.ScoreType.fa_rep
        scorefxn.set_weight(omega, 0.1)
        scorefxn.set_weight(fa_rep, 0.1)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 20)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(True)
        relax.apply(pose)
        scorefxn = pyrosetta.get_fa_scorefxn()
        self._cst_score(scorefxn)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, 20)
        relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
        relax.set_movemap(movemap)
        relax.cartesian(False)
        relax.apply(pose)

    def relax_simple(self, pose, cycles=15):
        scorefxn = pyrosetta.get_fa_scorefxn()
        self._cst_score(scorefxn)
        relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn, cycles)
        relax.apply(pose)
