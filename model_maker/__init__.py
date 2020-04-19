########################################################################################################################

__doc__ = \
    """
    ...
    """
__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2019 A.D."
__license__ = "MIT"
__version__ = "0.4"
__citation__ = ""

########################################################################################################################
from typing import Optional, List, Union

from ._base_mixin import _Base, pyrosetta  ### __init__ and basic io
from ._relax_mixin import _Relax  ### relax methods
from ._variants_mixin import _Make  ### making variants

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint


class Catcher(_Make):  # _Make <- _Relax <- _Base

    def get_score_panel(self, pose, save_variants=True, filename='pose'):
        scorefxn = pyrosetta.get_fa_scorefxn()
        catcher = self.make_catcher_only(pose)
        tag = self.make_tag_only(pose)
        unreacted = self.make_unreacted(pose)
        trans = self.make_transition(pose)
        d = {'isopeptide': scorefxn(pose),
             'Lazaridis-Karplus solvatation term': pose.energies().total_energies_array()['fa_sol'][0],
             'catcher': scorefxn(catcher),
             'tag': scorefxn(tag),
             'unreacted': scorefxn(unreacted),
             'trans': scorefxn(trans),
             'sequence': pose.sequence(),
             'scorefxn': pyrosetta.get_fa_scorefxn().get_name(),
             'pI': IsoelectricPoint(pose.sequence()).pi(),
             }
        if save_variants:
            pose.dump_pdb(f'{filename}.isopeptide.pdb')
            catcher.dump_pdb(f'{filename}.catcher.pdb')
            tag.dump_pdb(f'{filename}.tag.pdb')
            unreacted.dump_pdb(f'{filename}.unreacted.pdb')
            trans.dump_pdb(f'{filename}.transition.pdb')
        return d
