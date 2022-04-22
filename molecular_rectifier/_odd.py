########################################################################################################################

__doc__ = \
    """
    This add the oddity fixing functionality. Generally odd corner cases.
    Formerly part of collapse_ring.py
    """

########################################################################################################################

from rdkit import Chem

from ._base import _RectifierBase


class _RectifierOdd(_RectifierBase):

    def prevent_oddities(self):
        self._prevent_allene()

    # ====== Private methods ===========================================================================================

    def _prevent_allene(self):
        for atom in self.rwmol.GetAtoms():
            if atom.GetAtomicNum() < 14:
                n = []
                for bond in atom.GetBonds():
                    if bond.GetBondType().name in ('DOUBLE', 'TRIPLE'):
                        n.append(bond)
                    else:
                        pass
                if len(n) > 2:
                    # this is a mess!
                    self.log.info(f'Allene issue: {n} double bonds on {atom.GetSymbol()} atom {atom.GetIdx()}!')
                    for bond in n:
                        bond.SetBondType(Chem.BondType().SINGLE)
                elif len(n) == 2:
                    # downgrade the higher bonded one!
                    others = [a for bond in n for a in (bond.GetBeginAtom(), bond.GetEndAtom()) if
                              a.GetIdx() != atom.GetIdx()]
                    others = sorted(others, key=lambda atom: sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()]))
                    self.log.info(f'Allene removed between {atom.GetIdx()} and {[a.GetIdx() for a in others]}')
                    self.rwmol.GetBondBetweenAtoms(atom.GetIdx(), others[-1].GetIdx()).SetBondType(Chem.BondType.SINGLE)
                else:
                    pass
            else:
                continue
        self.modifications.append(self.mol)

    def _prevent_overclose(self):
        pass

