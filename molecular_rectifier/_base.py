########################################################################################################################

__doc__ = \
    """
    This has the methods common to ring and odd.
    Formerly part of collapse_ring.py
    """

########################################################################################################################

import logging
from typing import Tuple

from rdkit import Chem


class _RectifierBase:
    log = logging.getLogger(__name__)

    def __init__(self, mol: Chem.Mol,
                 atoms_in_bridge_cutoff: int = 2,
                 valence_correction: str = 'element',
                 debug: bool = False):
        """
        Instantiates but does not call ``fix`` (or its specific methods).

        :param mol: Does not get edited. ``.mol`` does (but is a``Chem.RWMol``, so use ``mol.GetMol()``)
        # how many bridge atoms can be deleted? (0 = preserves norbornane, 1 = preserves monsterantane)

        :param valence_correction:
        :param debug:
        """
        self.debug = bool(debug)
        if valence_correction not in ('charge', 'element'):
            raise ValueError(f'valence_correction "{valence_correction} id not charge/element')
        self.valence_correction = str(valence_correction)
        self.atoms_in_bridge_cutoff = int(atoms_in_bridge_cutoff)
        self.rwmol = Chem.RWMol(mol)
        self.modifications = []  # keeping track of steps
        self._valence_mode = 'max'
        self._iterations_done = 0
        self._subiterations_done = 0

    @property
    def mol(self) -> Chem.Mol:
        return self.rwmol.GetMol()

    @property
    def name(self) -> str:
        if self.rwmol.HasProp("_Name"):
            return self.rwmol.GetProp("_Name")
        else:
            return 'compound'

    def _get_ring_info(self, mode='atom') -> Tuple[Tuple[int]]:
        """
        you cannot get ring info on an unsanitized mol. Ironically I need ring info for sanitization

        :param mode: bond|atom
        :return: same as mol.GetRingInfo().AtomRings() or .BondRings()
        """
        mol2 = Chem.Mol(self.mol)
        for bond in mol2.GetBonds():
            bond.SetBondType(Chem.BondType.UNSPECIFIED)
        for atom in mol2.GetAtoms():
            atom.SetIsAromatic(False)
            atom.SetAtomicNum(0)
        Chem.SanitizeMol(mol2)
        if mode == 'atom':
            return mol2.GetRingInfo().AtomRings()
        elif mode == 'bond':
            return mol2.GetRingInfo().BondRings()
        else:
            raise ValueError(f'Unknown mode {mode}')


    def _get_atom_valence(self, atom: Chem.Atom):
        """
        Cannot get the normal way as it cannot be sanitised.

        :param atom:
        :return:
        """
        valence = 0
        for bond in atom.GetBonds():
            valence += bond.GetBondTypeAsDouble()
        return valence - atom.GetFormalCharge()
