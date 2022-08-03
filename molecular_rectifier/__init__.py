from __future__ import annotations
########################################################################################################################

__doc__ = \
    """
Fix issue in auto-merging.
    """

__author__ = "Matteo Ferla. [Github](https://github.com/matteoferla)"
__email__ = "matteo.ferla@gmail.com"
__date__ = "2021 A.D."
__license__ = "MIT"
__version__ = "0.1.8"
__citation__ = "See Fragmenstein"

########################################################################################################################

from ._base import _RectifierBase  # provides the __init__ and shared methods
from ._ring import _RectifierRing  # fixes rings
from ._odd import _RectifierOdd  # fixes specific oddities
from ._valence import _RectifierValence  # fixes valence
from ._overclose import _RectifierOverclose  # fixes extreme overlaps
from rdkit import Chem
from rdkit.Chem import AllChem

class Rectifier(_RectifierOverclose, _RectifierRing, _RectifierOdd, _RectifierValence):
    """
    Fixes the nastiness.

    Do note that that the Chem.Mol does not get modified in place.
    ``.rwmol`` does and is a Chem.RWMol. The ``.mol`` is a Chem.Mol.

    The steps can be found in ``.modifications``.

    The .log log is not given a handler.

    New atoms with have the bool prop ``_Novel``.

    Does not link distant atoms. For that see joining methods in Monster.

    >>> Rectifier(mol).fix().mol
    """

    def fix(self, catchErrors:bool=False, iteration:int=0) -> Rectifier:
        """
        `CatchErrors` is the command in Sanitize, but it is not used there.
        """
        self.log.debug('============== Rectifier fix started =============')
        self.log.debug(self.mol_summary)
        self.absorb_overclose(distance_cutoff=0.8)  # from _RectifierOverclose
        self.log.debug(self.mol_summary)
        self.log.debug('------------ overclose safeguard done ------------')
        try:
            if self.has_issues():
                self._fix_aromatic_rings()
                Chem.SetAromaticity(self.rwmol)
                flags = AllChem.SANITIZE_ALL ^ AllChem.SANITIZE_KEKULIZE
                self._update_cache(sanitize=True, flags=flags)
                self.log.debug(self.mol_summary)
                self.log.debug('------------ Aromatic correction done ------------')
        except Exception as error:
            self.log.info(f'{error.__class__.__name__}: {error}')
            self._update_cache(sanitize=False)
        # self._preemptive_protonate() this is a weird test
        self.fix_rings()  # from _RectifierRing
        self.log.debug(self.mol_summary)
        self.log.debug('------------ ring fixes done ------------')
        self.prevent_oddities()  # from _RectifierOdd
        self.ununspecified_bonds()  # from _RectifierValence
        self.triage_rings()  # from _RectifierValence
        self.log.debug(self.mol_summary)
        self.log.debug('------------ oddity fixes done ------------')
        Chem.Cleanup(self.rwmol)
        self._adjust_Hs()
        self.log.debug(self.mol_summary)
        self.log.debug('------------ hydrogens done ------------')
        if self.has_issues():
            self.fix_issues()  # from _RectifierValence
        self.log.debug(self.mol_summary)
        self.log.debug('------------ valence drama done ------------')
        self._update_cache(sanitize=True)
        # round trip to RWMol
        if self.has_issues():
            self.log.debug('Some issues remain.')
            self.fix_aromatic_radicals()  # from _RectifierRing
        self.rwmol = Chem.RWMol(self.mol)  # noqa it is a word
        self.no_radicals()  # just in case
        self.log.debug(self.mol_summary)
        self.log.debug('------------ check against radicals done ------------')
        iteration += 1
        if not self.has_issues():
            pass
        if iteration < 3:
            self.log.debug('Some issues remain... Trying again!')
            self.fix(catchErrors=catchErrors, iteration=iteration)
        elif iteration == 3:
            self.strip_hydrogens()
            problems = Chem.DetectChemistryProblems(self.rwmol)
            for p in problems:
                if p.GetType() == 'KekulizeException':
                    for i in p.GetAtomIndices():
                        self.log.debug(f'KekulizeException downgrade_substituents last resort: {i}')
                        self.downgrade_substituents(self.rwmol.GetAtomWithIdx(i))
            self.log.debug('Some issues remain. Stripping hydrogens...')
            self.fix(catchErrors=catchErrors, iteration=iteration)
        elif iteration == 4:
            problems = Chem.DetectChemistryProblems(self.rwmol)
            for p in problems:
                if p.GetType() == 'KekulizeException':
                    for i in p.GetAtomIndices():
                        self.log.debug(f'KekulizeException downgrade_ring last resort: {i}')
                        self.downgrade_ring(self.rwmol.GetAtomWithIdx(i), hard=True)
            self.fix(catchErrors=catchErrors, iteration=iteration)
        else:
            pass # burn...
        # check:
        CaughtError = Exception if catchErrors else ()  # noqa Uppercase as its a class
        try:
            self._update_cache(sanitize=False)  # catchErrors is true...
            Chem.SanitizeMol(self.rwmol, catchErrors=False)
        except CaughtError as error:
            self.log.info(f'{error.__class__.__name__}: {error}')
        return self
