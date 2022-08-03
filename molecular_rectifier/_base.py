########################################################################################################################

__doc__ = \
    """
    This has the methods common to ring and odd.
    Formerly part of collapse_ring.py
    """

########################################################################################################################

import logging
from typing import Tuple, Union, Sequence, List, Set
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

    def downgrade_ring(self, atom: Chem.Atom, hard=False):
        ## very crappy way of doing this
        self.log.info(f'downgrading whole ring due to {atom.GetSymbol()} atom i={atom.GetIdx()}')
        atom.SetIsAromatic(False)
        ringinfo = self._get_ring_info(mode='atom')  # == sane_mol.GetRingInfo().AtomRings()
        get_atomrings = lambda ai: [ring for ring in ringinfo if ai in ring]
        atomrings = get_atomrings(atom.GetIdx())
        for atomring in atomrings:
            rnieghs = self._get_ring_neighbors(atomring)  # list of pairs of indices that are neighbors in the ring
            for n1, n2 in rnieghs:
                for ai in (n1, n2):
                    atom: Chem.Atom = self.rwmol.GetAtomWithIdx(ai)
                    atom.SetIsAromatic(False)
                self.rwmol.GetBondBetweenAtoms(n1, n2).SetBondType(Chem.BondType.SINGLE)
        if hard:
            return
        # add double bonds strategically
        for atomring in atomrings:
            rnieghs = self._get_ring_neighbors(atomring)
            for n1, n2 in rnieghs:
                if self._get_valence_difference(self.rwmol.GetAtomWithIdx(n1)) <= -2 and \
                        self._get_valence_difference(self.rwmol.GetAtomWithIdx(n2)) <= -2:
                    self.rwmol.GetBondBetweenAtoms(n1, n2).SetBondType(Chem.BondType.DOUBLE)

    def downgrade_substituents(self, atom: Chem.Atom) -> bool:
        """
        Prequel to ``downgrade_ring``.
        Downgrades double bonded substituents even if they are 4n+2 obeying...
        """
        change = False
        ai: int = atom.GetIdx()
        bondrings = [br for br in self._get_ring_info(mode='bond') if ai in br]
        if len(bondrings) == 0:
            return False  # it is not in a ring
        for bond_idxs in bondrings:
            bonds: List[Chem.Bond] = list(map(self.rwmol.GetBondWithIdx, bond_idxs))
            # expand bonds to substituents
            exbond_idxs: Set[int] = set()
            for bond in bonds:
                for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                    exbond_idxs.update(map(Chem.Bond.GetIdx, atom.GetBonds()))
            subbonds: List[Chem.Bond] = list(map(self.rwmol.GetBondWithIdx, set(exbond_idxs).difference(bond_idxs)))
            for bond in subbonds:
                if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                    bond.SetBondType(Chem.BondType.SINGLE)
                    change = True
        return change

    def strip_hydrogens(self):
        """
        This is not ideal... last resort only
        """
        for i in list(range(self.rwmol.GetNumAtoms() - 1, 0 - 1, -1)):
            if self.rwmol.GetAtomWithIdx(i).GetSymbol() == 'H':
                self.rwmol.RemoveAtom(i)
        for a in self.rwmol.GetAtoms():
            a.SetNoImplicit(False)

    @property
    def mol_summary(self) -> str:
        parts = [f'Fragments: {len(Chem.GetMolFrags(self.rwmol))}',
                 'Bonds: ' +
                 str(Counter([bond.GetBondType().name for bond in self.rwmol.GetBonds()]).most_common()),
                 'Atoms: '+
                 str(Counter([atom.GetSymbol() for atom in self.rwmol.GetAtoms()]).most_common()),
                 ]
        return '; '.join(parts)

    def _update_cache(self, sanitize=False, flags:Chem.SanitizeFlags=Chem.SanitizeFlags.SANITIZE_ALL):
        if hasattr(self.rwmol, 'UpdatePropertyCache'):
            self.rwmol.UpdatePropertyCache(strict=False)
        if sanitize:
            Chem.SanitizeMol(self.rwmol, sanitizeOps=flags, catchErrors=True)

    def has_issues(self) -> bool:
        """
        ``Chem.DetectChemistryProblems(mol)`` is like running a test of ``Chem.SanitizeMol(mol)``,
        but it does not count radicals as a problem... in fact SanitizeMol will introduce radicals!
        This can be mimicked by ``Chem.AssignRadicals(mol)``.
        Therefore this method returns true if there is an issue according to ``Chem.DetectChemistryProblems(mol)``
        but also if sanitization would result in radicals.
        """
        if Chem.DetectChemistryProblems(self.rwmol):
            return True
        mol:Chem.Mol = self.mol  # this is a copy
        Chem.AssignRadicals(mol)
        return bool(sum(map(Chem.Atom.GetNumRadicalElectrons, mol.GetAtoms())))
