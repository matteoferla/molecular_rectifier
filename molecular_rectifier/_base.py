########################################################################################################################

__doc__ = \
    """
    This has the methods common to ring and odd.
    Formerly part of collapse_ring.py
    """

########################################################################################################################

import logging
from typing import Tuple, Union
from collections import Counter
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
        if not isinstance(mol, Chem.Mol):
            # anything that is not Chem.Mol or Chem.RWMol
            # isinstance(Chem.RWMol(), Chem.Mol) == True
            raise TypeError('not a molecule')
        # an unsanitised mol has no `.updatePropertyCache`
        # mol.updatePropertyCache(strict=False)
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

    def _downgrade_aromatic_bond(self, bond: Union[Chem.Bond, int]):
        """
        This is a bad way.
        """
        if isinstance(bond, int):
            bond: Chem.Bond = self.rwmol.GetBondWithIdx(bond)
        bond.SetBondType(Chem.BondType.DOUBLE)
        for atom in (bond.GetBeginAtom(), bond.GetBeginAtom()):
            atom.SetIsAromatic(False)

    def downgrade_ring(self, atom: Chem.Atom):
        ## very crappy way of doing this
        self.log.info(f'downgrading whole ring due to {atom.GetSymbol()} atom i={atom.GetIdx()}')
        atom.SetIsAromatic(False)
        ringinfo = self._get_ring_info(mode='atom')  # == sane_mol.GetRingInfo().AtomRings()
        get_atomrings = lambda ai: [ring for ring in ringinfo if ai in ring]
        atomrings = get_atomrings(atom.GetIdx())
        for atomring in atomrings:
            rnieghs = self._get_ring_neighbors(atomring) # list of pairs of indices that are neighbors in the ring
            for n1, n2 in rnieghs:
                for ai in (n1, n2):
                    atom: Chem.Atom = self.rwmol.GetAtomWithIdx(ai)
                    atom.SetIsAromatic(False)
                self.rwmol.GetBondBetweenAtoms(n1, n2).SetBondType(Chem.BondType.SINGLE)
        for atomring in atomrings:
            rnieghs = self._get_ring_neighbors(atomring)
            for n1, n2 in rnieghs:
                if self._get_valence_difference(self.rwmol.GetAtomWithIdx(n1)) <= -2 and \
                        self._get_valence_difference(self.rwmol.GetAtomWithIdx(n2)) <= -2:
                    self.rwmol.GetBondBetweenAtoms(n1, n2).SetBondType(Chem.BondType.DOUBLE)

    @property
    def mol_summary(self) -> str:
        parts = [f'Fragments: {len(Chem.GetMolFrags(self.rwmol))}',
                 'Bonds: ' +
                 str(Counter([bond.GetBondType().name for bond in self.rwmol.GetBonds()]).most_common()),
                 'Atoms: '+
                 str(Counter([atom.GetSymbol() for atom in self.rwmol.GetAtoms()]).most_common()),
                 ]
        return '; '.join(parts)

