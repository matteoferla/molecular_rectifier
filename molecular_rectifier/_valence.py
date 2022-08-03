########################################################################################################################

__doc__ = \
    """
    This add the oddity fixing functionality. Generally odd corner cases.
    Formerly part of collapse_ring.py
    """



########################################################################################################################
from rdkit import Chem
from rdkit.Chem import AllChem
import itertools
from typing import List, Tuple, Union
from ._base import _RectifierBase


# when hits are combined they can result in odd valence and other issues.

class _RectifierValence(_RectifierBase):
    """
        Checks whether the valence is right and corrects by either shifting the element or adding a charge.
        With the following exceptions:

        * Texas carbon -> Sulfur
        * Hydrogen -> Fluoride shifted downwards
        * Carbon in aromatic ring -> single bond ring
    """
    # ========= Three main steps =======================================================================================

    def ununspecified_bonds(self) -> None:
        """
        This get run during instantiation. It is the first one run.
        BondType.UNSPECIFIED ==> BondType.SINGLE

        :return: None
        """
        for bond in self.rwmol.GetBonds():
            if bond.GetBondType().name == 'UNSPECIFIED':
                self.log.debug(f'Fixing unspecified bond {bond.GetIdx()}')
                bond.SetBondType(Chem.BondType.SINGLE)
                # debug:
                self.modifications.append(self.mol)
        self._update_cache()

    def triage_rings(self) -> None:
        """
        This get run during instantiation. It is the second one run.
        It is also run again by `fix_issues` the third step after a weird valence change.
        Deals with rings.

        :return: None
        """
        # upgrade to aromatic if aromatic.
        rings = self._get_ring_info()
        for ring in rings:
            if self._is_aromatic_ring(ring, rings):
                for i, ni in self._get_ring_neighbors(ring):
                    self.rwmol.GetBondBetweenAtoms(i, ni).SetBondType(Chem.BondType.AROMATIC)
        # downgrade to single if non-ring aromatic
        for i, atom in enumerate(self.rwmol.GetAtoms()):
            if any([i in r for r in rings]):
                continue  # ring
            else:  # non-ring
                for bond in atom.GetBonds():
                    if bond.GetBondType().name == 'AROMATIC':
                        self.log.debug(f'downgrading bond {i}')
                        bond.SetBondType(Chem.BondType.SINGLE)
        # aromatics
        for ring in sorted(rings, key=self._is_aromatic_ring):
            if self._is_aromatic_ring(ring, rings):  # the nonaromatic rings will be done first.
                for i in ring:
                    self.rwmol.GetAtomWithIdx(i).SetIsAromatic(True)
            else:
                for i in ring:
                    self.rwmol.GetAtomWithIdx(i).SetIsAromatic(False)
        self._update_cache()
        self.modifications.append(self.mol)  # may not have changed.

    def fix_issues(self, _previous=None) -> None:
        """
        This get run during instantiation. It is the third and final one run before sanitization.
        Deals with a variety of problems.
        It calls itself until no problems according to `Chem.DetectChemistryProblems` exits.
        It is a bit shoddy and any oddity likely steps from here. TODO

        :return: None
        """
        self.modifications.append(self.mol)  # may not have changed.
        problems = Chem.DetectChemistryProblems(self.rwmol)
        if self._iterations_done > 100:
            self.log.error(f'Iterations maxed out!')
            return None
        elif self._subiterations_done > 5:
            self.log.error(f'Unfixable')
            return None
        elif len(problems) == 0:
            return None
        else:
            self.log.debug(f'(Iteration: {self._iterations_done}) N problems {len(problems)}')
            p = problems[0]
            self.log.debug(f'(Iteration: {self._iterations_done}) Issue {p.GetType()}: {p.Message()}')
            if p.Message() == _previous:
                self.triage_rings()
            ############################################################
            if p.GetType() == 'KekulizeException':
                if p.Message() != _previous:
                    N = self._get_nitrogens(p.GetAtomIndices())
                    if len(N) > 0 and self._nitrogen_protonate(N, p.Message()):
                        pass  # been fixed.
                    else:
                        # triage rings should have altered any not ring atoms that are aromatic.
                        # self._get_ring_info()
                        # so it is likely a hetatom thing.
                        self.log.info(f'Ring triages seems to have failed. Is it a valence thing?')
                        valence_issues = [self._has_correct_valence(i) for i in p.GetAtomIndices()]
                        if not all(valence_issues):
                            for i in p.GetAtomIndices():
                                self.fix_valence(i)
                        else:
                            self.log.warning(f'Attempting default valency (not max)')
                            self._valence_mode = 'default'
                            for i in p.GetAtomIndices():
                                self.fix_valence(i)
                            self._valence_mode = 'max'
                else:
                    for i in p.GetAtomIndices():
                        if self.downgrade_substituents(self.rwmol.GetAtomWithIdx(i)):
                            break
                        self.downgrade_ring(self.rwmol.GetAtomWithIdx(i))
                    self.triage_rings()
            ############################################################
            elif p.GetType() == 'AtomKekulizeException' and 'non-ring atom' in p.Message():
                atom = self.rwmol.GetAtomWithIdx(p.GetAtomIdx())
                atom.SetIsAromatic(False)
                self.log.debug(f'Atom {p.GetAtomIdx()} set to non-aromatic.')
                for bond in atom.GetBonds():
                    bond.SetBondType(Chem.BondType.SINGLE)
            elif p.GetType() == 'AtomKekulizeException' and 'Aromatic bonds on non aromatic atom' in p.Message():
                atom = self.rwmol.GetAtomWithIdx(p.GetAtomIdx())
                self.log.debug(f'Atom {p.GetAtomIdx()} set to aromatic.')
                atom.SetIsAromatic(True)
            ############################################################
            elif p.GetType() == 'AtomValenceException':
                i = p.GetAtomIdx()
                self.fix_valence(i)
            else:
                self.log.error('???', p.GetType(), p.Message())
            self._iterations_done += 1
            if _previous != p.Message():
                self.log.debug(f'{self._iterations_done} appears successful.')
                self._subiterations_done = 0
            else:
                self._subiterations_done += 1
                self.log.debug(f'{self._iterations_done} appears unsuccessful.')
            return self.fix_issues(_previous=p.Message())

    # ========= Methods that circumvent the nonsanitization ============================================================

    def _get_valence_difference(self, atom: Chem.Atom) -> int:
        pt = Chem.GetPeriodicTable()
        valence = self._get_atom_valence(atom)
        if self._valence_mode == 'max':
            maxv = max(pt.GetValenceList(atom.GetAtomicNum()))
            return valence - maxv
        else:
            d = pt.GetDefaultValence(atom.GetAtomicNum())
            return valence - d

    def _has_correct_valence(self, atom: Union[Chem.Atom, int]):
        if isinstance(atom, Chem.Atom):
            return self._get_valence_difference(atom) <= 0
        elif isinstance(atom, int):
            atom = self.rwmol.GetAtomWithIdx(atom)
            return self._get_valence_difference(atom) <= 0

    # ========= rings ==================================================================================================

    def _get_atoms_at_fusion(self, ring, rings):
        fused = []
        if rings is not None:
            for other in rings:
                if other == ring:
                    pass
                elif not set(ring).isdisjoint(other):
                    fused.extend(set(ring).intersection(other))
        return fused

    def _is_aromatic_ring(self, ring: Tuple[int], rings=None) -> bool:
        """
        :param ring: GetRingInfo().AtomRings() entry
        :return:
        """
        fused = self._get_atoms_at_fusion(ring, rings)
        for i in list(set(ring).difference(fused)):
            atom_i = self.rwmol.GetAtomWithIdx(i)
            for n in atom_i.GetNeighbors():
                ni = n.GetIdx()
                if ni in ring:
                    if self.rwmol.GetBondBetweenAtoms(i, ni).GetBondType().name == 'AROMATIC':
                        return True
        else:
            return False

    def _get_ring_neighbors(self, ring: Tuple[int]) -> List[Tuple[int, int]]:
        """
        :param ring: GetRingInfo().AtomRings() entry
        :return: list of pairs of indices that are neighbors in the ring
        """
        rns = []
        for i in ring:
            atom = self.rwmol.GetAtomWithIdx(i)
            for n in atom.GetNeighbors():
                ni = n.GetIdx()
                if ni in ring:
                    rns.append((i, ni))
        return rns

    def _get_aroma(self, atom, this_bond):
        # determine if the bond of the atom is aromatic
        return [b for b in atom.GetBonds() if b.GetIdx() != this_bond and b.GetBondType().name == 'AROMATIC']

    def _get_other(self, bond, these_atoms):
        others = [a for a in (bond.GetBeginAtom(), bond.GetEndAtom()) if a.GetIdx() not in these_atoms]
        if others:
            other = others[0]
            other.SetIsAromatic(False)
            return other

        # if len(self._get_rings(atom.GetIdx())) == 1:
        #     for bond in atom.GetBonds():
        #         bond.SetBondType(Chem.BondType.SINGLE)
        #         other = self._get_other(bond, [atom.GetIdx()])
        #         aro = self._get_aroma(other, bond.GetIdx())
        #         if aro:
        #             aro[0].SetBondType(Chem.BondType.DOUBLE)
        #             doubleother = self._get_other(aro[0], [atom.GetIdx(), other.GetIdx()])
        #             for b in doubleother.GetBonds():
        #                 if b.GetBondType() == Chem.BondType.AROMATIC:
        #                     b.SetBondType(Chem.BondType.SINGLE)
        #                 neigh = self._get_other(b, [doubleother.GetIdx()])
        #                 if neigh:
        #                     neigh.SetIsAromatic(False)

    # ========= Sanitization based fixes ===============================================================================

    def _nitrogen_protonate(self, nitrogens, previous):
        """

        :param nitrogens: list of Nitrogens
        :param previous:
        :return:
        """
        def reset():
            for n in nitrogens:
                self.rwmol.GetAtomWithIdx(n).SetNumExplicitHs(0)

        reset()
        p = Chem.DetectChemistryProblems(self.rwmol)
        if len(p) == 0 or p[0].Message() != previous:
            return True
        for i in range(1, len(nitrogens)):
            for c in itertools.combinations(nitrogens, i):
                reset()
                for n in c:
                    self.rwmol.GetAtomWithIdx(n).SetNumExplicitHs(1)
                p = Chem.DetectChemistryProblems(self.rwmol)
                if len(p) == 0 or p[0].Message() != previous:
                    return True
        return False

        self.log.debug(f'KekulizeException likely caused by nitrogen')

    # ========= other helpers ==========================================================================================

    def _get_nitrogens(self, indices):
        """
        Called when ``KekulizeException`` happends during ``.fix_issues``

        :param indices:
        :return:
        """
        return [i for i in indices if self.rwmol.GetAtomWithIdx(i).GetSymbol() == 'N']

    # ========= shift/charge ===========================================================================================

    def _adjust_for_fix_valence(self, atom):
        df = self._get_valence_difference(atom)
        ori = atom.GetSymbol()
        if self.valence_correction == 'charge':
            atom.SetFormalCharge(df)
        elif self.valence_correction == 'element':
            # ## correct row
            n = atom.GetAtomicNum()
            if n == 1:
                atom.SetAtomicNum(8)
            elif n > 10:
                n = (n % 8) - 2 + 8
                atom.SetAtomicNum(n)
            # ## correct column
            if len(atom.GetNeighbors()) > 4:
                self._break_bonds(atom)
            # elif len(atom.GetNeighbors()) > 4 and n <= 16: # S...
            #     atom.SetAtomicNum(16)
            elif n - df < 6:  # C -> B no!
                for bond in atom.GetBonds():
                    bond.SetBondType(Chem.BondType.SINGLE)
            else:  # N, O, F etc.
                atom.SetAtomicNum(int(n - df))
            self.log.info(f'Shifting atom from {ori} to {atom.GetSymbol()}')
        else:
            raise ValueError(f'self.valence_correction can only be "element"/"charge" not {self.valence_correction}.')

    def _break_bonds(self, atom):
        """
        Extreme last ditch. Breaks off all non-ring bonds to atom.
        Will likely trigger emergency_joining.

        :param atom:
        :return:
        """
        self.log.warning(f'In molecule ({self.name}) breaking bond to atom {atom.GetIdx()}')
        ring_indices = [a for ring in self._get_ring_info() for a in ring]
        for neigh in atom.GetNeighbors():
            if neigh.GetIdx() in ring_indices or neigh.GetSymbol() == '*':
                continue
            else:
                self.rwmol.RemoveBond(atom.GetIdx(), neigh.GetIdx())


    def fix_valence(self, i: int):
        atom = self.rwmol.GetAtomWithIdx(i)
        atom.SetFormalCharge(0)
        atom.SetNumExplicitHs(0)
        self.log.debug(f'{i} {atom.GetSymbol()}: {len(atom.GetNeighbors())} bonds {self._get_atom_valence(atom)}')
        if self._has_correct_valence(atom):
            self.log.debug('\tValence seems correct')
            return None
        elif atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) > 4:
            self.log.debug('\ttexas carbon --> S')
            atom.SetAtomicNum(16)
        elif atom.GetSymbol() == 'C' and atom.GetIsAromatic() and len(atom.GetNeighbors()) == 4:
            self.log.debug('\tDowngrading ring')
            self.downgrade_ring(atom)
        elif atom.GetSymbol() == 'C' and not atom.GetIsAromatic():
            for bond in atom.GetBonds():
                bond.SetBondType(Chem.BondType.SINGLE)
        else:
            self._adjust_for_fix_valence(atom)
        # did it work?
        if self._has_correct_valence(atom):
            return self.rwmol
        else:
            return self.fix_valence(i)

    def _adjust_Hs(self):
        for atom in self.rwmol.GetAtoms():
            atom.SetNumRadicalElectrons(0)
        self._update_cache(sanitize=False)
        idxs = [i for ring in self._get_ring_info('atom') for i in ring]
        for i in idxs:
            # atom: Chem.Atom = self.rwmol.GetAtomWithIdx(i)
            # self.fix_valence(i)
            pass
        mol = AllChem.AddHs(self.rwmol, addCoords=bool(self.rwmol.GetNumConformers()))
        self.rwmol = Chem.RWMol(mol)

    def _preemptive_protonate(self):
        """
        This is a debug test in essence as it's a bad idea
        """
        for atom in self.rwmol.GetAtoms():
            if atom.GetSymbol() == 'N' and atom.GetIsAromatic() and self._get_atom_valence(atom) < 4:
                atom.SetFormalCharge(1)
                atom.SetNumExplicitHs(1)
                atom.SetNumRadicalElectrons(0)
