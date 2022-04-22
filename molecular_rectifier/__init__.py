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
__version__ = "0.1"
__citation__ = "See Fragmenstein"

########################################################################################################################

from ._base import _RectifierBase  # provides the __init__ and shared methods
from ._ring import _RectifierRing  # fixes rings
from ._odd import _RectifierOdd  # fixes specific oddities
from ._valence import _RectifierValence  # fixes valence
from ._overclose import _RectifierOverclose  # fixes extreme overlaps
from rdkit import Chem


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

    def fix(self) -> Rectifier:
        self.absorb_overclose()  # from _RectifierOverclose
        self.fix_rings()  # from _RectifierRing
        self.prevent_oddities()  # from _RectifierOdd
        self.ununspecified_bonds()  # from _RectifierValence
        self.triage_rings()  # from _RectifierValence
        Chem.Cleanup(self.rwmol)
        self.fix_issues()  # from _RectifierValence
        Chem.SanitizeMol(self.rwmol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)
        return self
