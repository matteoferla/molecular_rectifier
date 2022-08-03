########################################################################################################################

__doc__ = \
    """
    This add the ring fixing functionality. Fusing etc.
    Formerly part of collapse_ring.py
    """

########################################################################################################################

import itertools
from collections import Counter
from rdkit.Chem import rdqueries
from typing import Optional, Tuple, Dict, Union, List, Set

from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D

from ._base import _RectifierBase


# ======================================================================================================================


class _RectifierRing(_RectifierBase):

    def fix_rings(self):
        self._prevent_bonded_to_bridgeheads()
        self._fix_aromatic_rings()
        self.fix_aromatic_radicals()  # unlikely to not have an effect
        self._prevent_conjoined_ring()
        self._prevent_weird_rings()

    def fix_aromatic_radicals(self):
        """
        This happens on a super corner case and I do not know why.
        """
        self._update_cache(sanitize=True, flags=Chem.SanitizeFlags.SANITIZE_FINDRADICALS)
        ringbonds = self._get_ring_info(mode='bond')  # == sane_mol.GetRingInfo().BondRings()
        for ring_i, bonds_ring in enumerate(ringbonds): # self.rwmol.GetBonds():
            for bond_i in bonds_ring:
                bond = self.rwmol.GetBondWithIdx(bond_i)
                if not bond.GetEndAtom().GetNumRadicalElectrons():
                    continue
                if not bond.GetBeginAtom().GetNumRadicalElectrons():
                    continue
                self.log.info('Radical in aromatic structure present')
                # both sides are radicals...
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    bond.SetBondType(Chem.BondType.DOUBLE)
                    for atom in (bond.GetBeginAtom(), bond.GetEndAtom()):
                        atom.SetNumRadicalElectrons(0)
                        self._update_cache(sanitize=False)
                    # self._downgrade_aromatic_bond(bond_j) not needed.
        # add hydrogen?
        # rdqueries does not work on unsanitary molecules!
        #is_radical = rdqueries.NumRadicalElectronsEqualsQueryAtom(1)
        # for atom in self.rwmol.GetAtoms(): #: Chem.Atom
        #     if atom.GetSymbol() != 'C' and atom.GetNumRadicalElectrons():
        #         atom.SetFormalCharge(+1)
        #         atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
        #     atom.SetNumRadicalElectrons(0)

    def no_radicals(self):
        for atom in self.rwmol.GetAtoms():
            atom.SetNumRadicalElectrons(0)
        self._update_cache(sanitize=False)

    def _prevent_bonded_to_bridgeheads(self):
        """
        A 3-aromoatic-bond carbon cannot be bonded to anything else (except if a spiro)
        It removes the offending bond but does nothig more...
        """
        # list of dict shared, ring_atom_A_idxs, ring_atom_B_idxs
        bridge_infos = self._get_bridges()
        for info in bridge_infos:
            for idx in info['shared']:
                atom: Chem.Atom = self.rwmol.GetAtomWithIdx(idx)
                # middle atom:
                if len(atom.GetNeighbors()) <= 2:
                    continue
                # heavy-heavy atom:
                elif atom.GetAtomicNum() > 8:
                    continue
                # aromatic atom:
                elif atom.GetIsAromatic() or \
                   any([bond.GetBondType() == Chem.BondType.AROMATIC for bond in atom.GetBonds()]):
                    neigh_idxs: Set[int] = {neigh.GetIdx() for neigh in atom.GetNeighbors()}
                    other_idxs = neigh_idxs.difference(info['ring_atom_A_idxs'])\
                                           .difference(info['ring_atom_B_idxs'])
                    for idx in other_idxs:
                        self.log.info('non ring atom bonded to bridgehead')
                        self.rwmol.RemoveBond(atom.GetIdx(), idx)

    def make_all_rings_aromatic(self):
        """
        This is an extreme debug method not called by fix.
        """
        for i, bonds in enumerate(self.rwmol.GetRingInfo().BondRings()):
            if any([self.rwmol.GetBondWithIdx(j).GetBondType() == Chem.BondType.AROMATIC for j in bonds]):
                for j in bonds:
                    self.rwmol.GetBondWithIdx(j).SetBondType(Chem.BondType.AROMATIC)

    def _prevent_conjoined_ring(self) -> None:
        """
        This kills bridging bonds with not atoms in the bridge within rings.
        So it is bridged, fused and spiro safe.
        It removes only one bond, so andamantane/norbornane are safe.
        """
        c = Counter([i for ring in self._get_ring_info() for i in ring])
        nested = [k for k in c if c[k] >= 3]
        pairs = [(idx_a, idx_b) for idx_a, idx_b in itertools.combinations(nested, r=2) if
                 self.rwmol.GetBondBetweenAtoms(idx_a, idx_b) is not None]
        rank = sorted(pairs, key=lambda x: c[x[0]] + c[x[1]], reverse=True)
        if len(rank) > 0:
            idx_a, idx_b = rank[0]
            self.rwmol.RemoveBond(idx_a, idx_b)  # SetBoolProp('_IsRingBond') is not important
            self.log.info(f'Zero-atom bridged ring issue: bond between {idx_a}-{idx_b} removed')
            # re-run:
            self._prevent_conjoined_ring()
        self.modifications.append(self.mol)  # not self.rwmol as I want a snapshot

    def _get_bridges(self) -> List[Dict[str, Union[Tuple[int]]]]:
        """
        The output includes not only the bridgeheads but the middle ones too
        list of dict shared, ring_atom_A_idxs, ring_atom_A_idxs
        """
        bridge_info = []
        ringatoms = self._get_ring_info(mode='atom')  # GetRingInfo().AtomRings()
        for ring_A, ring_B in itertools.combinations(ringatoms, r=2):
            shared = set(ring_A).intersection(set(ring_B))
            if len(shared) == 0:
                continue
            bridge_info.append(dict(shared=shared,
                                    ring_atom_A_idxs=ring_A,
                                    ring_atom_B_idxs=ring_B)
                                   )
        return bridge_info


    def _prevent_weird_rings(self, iteration=0):
        bridge_info: List[Dict[str, Tuple[int]]] = self._get_bridges()
        for info in bridge_info:
            ring_A: Tuple[int] = info['ring_atom_A_idxs']
            ring_B: Tuple[int] = info['ring_atom_B_idxs']
            shared: Tuple[int] = info['shared']
            if len(shared) == 0:
                # no longer accessible
                self.log.debug(f'This molecule ({self.name}) has some separate rings')
                pass  # separate rings
            elif len(shared) < self.atoms_in_bridge_cutoff and \
                    self.atoms_in_bridge_cutoff >= 2 \
                    and len(ring_A) == len(ring_B):
                # adamantene/norbornane/tropinone kind of thing
                self.log.warning(f'This molecule ({self.name}) has a bridge: leaving/spiro')
                pass  # ideally check if planar...
            elif len(shared) == 1:
                self.log.debug(f'This molecule ({self.name}) has a spiro bicycle')
                pass  # spiro ring.
            elif len(shared) == 2:
                self.log.debug(f'This molecule ({self.name}) has a fused ring')
                if self.rwmol.GetBondBetweenAtoms(*shared) is not None:
                    pass  # indole/naphtalene
                    small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
                    if len(small) == 4:
                        self.log.warning(f'This molecule ({self.name}) ' +
                                             'has a benzo-azetine–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1CC2')
                        # benzo-azetine is likely an error: add and extra atom
                        a, b = set(small).difference(big)
                        self._place_between(a, b)
                    elif len(small) == 3:
                        self.log.warning(f'This molecule ({self.name}) '+
                                             'has a benzo-cyclopropane–kind-of-thing: expanding to indole')
                        # Chem.MolFromSmiles('C12CCCCC1C2')
                        # benzo-cyclopronane is actually impossible at this stage.
                        a = list(set(small).difference(big))[0]
                        for b in shared:
                            self._place_between(a, b)
                    else:
                        pass  # indole and nathalene
                elif (len(ring_A), len(ring_B)) == (6, 6):
                    raise Exception('This is utterly impossible')
                else:
                    self.log.warning(f'mysterious ring system {len(ring_A)} + {len(ring_B)}')
                    pass  # ????
            elif len(shared) < self.atoms_in_bridge_cutoff:
                # adamantene/norbornane/tropinone kind of thing
                self.log.warning(f'This molecule ({self.name}) has a bridge: leaving/spiro')
                pass  # ideally check if planar...
            else:
                self.log.warning(f'This molecule ({self.name}) has a bridge that will be removed')
                self._prevent_bridge_ring(ring_A)
                # start from scratch.
                if iteration < 3:
                    self._prevent_weird_rings(iteration=iteration+1)
                else:
                    self.log.warning(f'Too many trials to remove bridge.')
        self.modifications.append(self.mol)  # not self.rwmol as I want a snapshot

    # ===== Dependencies of prevent weird rings =================================================================================

    def _place_between(self, a: int, b: int, aromatic: Optional[bool] = None, atomic_number: int = 6) -> None:
        """
        Places an C atom, possibly of type aromatic, between atom of index a, and of b.

        :param a: index of atom A
        :param b: index of atom B
        :param aromatic: bool of aromaticity (False = Single, None = copy, True = aromatic)
        :param atomic_number: Carbon is 6.
        :return:
        """
        oribond = self.rwmol.GetBondBetweenAtoms(a, b)
        if oribond is None:
            self.log.critical(f'FAIL. There should be a bond btween {a} and {b}')
            return None  # fail
        elif aromatic is True:
            bt = Chem.BondType.AROMATIC
        elif aromatic is False:
            bt = Chem.BondType.SINGLE
        else:
            bt = oribond.GetBondType()
        idx = self.rwmol.AddAtom(Chem.Atom(atomic_number))
        neoatom = self.rwmol.GetAtomWithIdx(idx)
        atom_a = self.rwmol.GetAtomWithIdx(a)
        atom_b = self.rwmol.GetAtomWithIdx(b)
        if aromatic:
            neoatom.SetIsAromatic(True)
            atom_a.SetIsAromatic(True)
            atom_b.SetIsAromatic(True)
        # prevent constraints
        neoatom.SetBoolProp('_Novel', True)
        atom_a.SetBoolProp('_Novel', True)
        atom_b.SetBoolProp('_Novel', True)
        # fix position
        conf = self.rwmol.GetConformer()
        pos_A = conf.GetAtomPosition(a)
        pos_B = conf.GetAtomPosition(b)
        x = pos_A.x / 2 + pos_B.x / 2
        y = pos_A.y / 2 + pos_B.y / 2
        z = pos_A.z / 2 + pos_B.z / 2
        conf.SetAtomPosition(idx, Point3D(x, y, z))
        # fix bonds
        self.rwmol.RemoveBond(a, b)
        self.rwmol.AddBond(a, idx, bt)
        self.rwmol.AddBond(b, idx, bt)

    def _prevent_bridge_ring(self, examplar: Tuple[int]) -> None:
        # examplar is ring
        ringatoms = self._get_ring_info()  # GetRingInfo().AtomRings()
        ringatoms = [ring for ring in ringatoms if set(ring).intersection(examplar)]
        ring_idx = list(range(len(ringatoms)))
        shared_count = {}
        for ra, rb in itertools.combinations(ring_idx, r=2):
            shared_count[(ra, rb)] = len(set(ringatoms[ra]).intersection(set(ringatoms[rb])))
        if len(shared_count) == 0:
            return None
        ra, rb = list(shared_count.keys())[0]
        shared = list(set(ringatoms[ra]).intersection(ringatoms[rb]))
        has_bond = lambda a, b: self.rwmol.GetBondBetweenAtoms(a, b) is not None
        pairs = [(a, b) for a, b in itertools.combinations(shared, r=2) if has_bond(a, b)]
        c = Counter([i for pair in pairs for i in pair])
        ring_A, ring_B = ringatoms[ra], ringatoms[rb]
        small, big = sorted([ring_A, ring_B], key=lambda ring: len(ring))
        inners = [i for i in c if c[i] > 1]
        x = list(set(shared).difference(inners))
        if len(x) != 2:
            self.log.critical(
                f'This is impossible. {ringatoms} share {shared} with {inners} in the inside and {x} on the edge?')
            return None
        a, b = x
        if len(big) > 6:
            self.log.warning(f'Removing {len(inners)} bridging atoms and replacing with fused ring')
            # bond the vertices
            bt = Chem.BondType.SINGLE  # ???
            if self.rwmol.GetBondBetweenAtoms(a, b) is None:
                self.rwmol.AddBond(a, b, bt)
            else:
                self.log.warning('This is really odd! Why is there a bond already??')
            # remove the middle atoms.
            for i in sorted(inners, reverse=True):
                self.rwmol.RemoveAtom(i)
        else:
            self.log.warning(f'Shriking the smaller ring to change from bridged to fused.')
            # get the neighbour in the small atom to a vertex.
            neighs = [neigh for neigh in self.rwmol.GetAtomWithIdx(a).GetNeighbors() if
                      neigh.GetIdx() not in shared and neigh.GetIdx() in small]
            neigh = sorted(neighs, key=lambda atom: atom.GetSymbol() != 'C')[0]
            bt = self.rwmol.GetBondBetweenAtoms(a, neigh.GetIdx()).GetBondType()
            self.rwmol.RemoveBond(a, neigh.GetIdx())
            new_neigh = [neigh for neigh in self.rwmol.GetAtomWithIdx(a).GetNeighbors() if neigh.GetIdx() in shared][0]
            self.rwmol.AddBond(neigh.GetIdx(), new_neigh.GetIdx(), bt)
            neigh.SetBoolProp('_Novel', True)
            new_neigh.SetBoolProp('_Novel', True)
            self.rwmol.GetAtomWithIdx(a).SetBoolProp('_Novel', True)

    def _fix_aromatic_rings(self):
        """
        Does not test for 4n+2. That is a headache due to charge of non-carbon elements...
        """
        ringbonds:Tuple[Tuple[int]] = self._get_ring_info(mode='bond')
        ringatoms: Tuple[Tuple[int]] = self._get_ring_info(mode='atom')
        get_bond = lambda i: self.rwmol.GetBondWithIdx(i)
        for i, bonds in enumerate(ringbonds):
            # find bonds that are not shared with other rings
            uniques = [bond_i for bond_i in bonds if sum([bond_i in ring for ring in ringbonds]) == 1]
            # ignore rings that are not aromatic
            if not any([get_bond(bond_i).GetBondType() == Chem.BondType.AROMATIC for bond_i in uniques]):
                continue
            # for special case for ring with 5 aromatic carbons
            composition: Set[str] = {get_bond(bond_i).GetBeginAtom().GetSymbol() for bond_i in bonds}.union(
                                            {get_bond(bond_i).GetEndAtom().GetSymbol() for bond_i in bonds}
                                        )
            # for special case of spiro rings
            if self.check_for_spiro(i, ringatoms):
                self._fix_nonaromatic(uniques)
            elif len(bonds) == 5 and composition == {'C'}:
                self._fix_penta_aromatic(bonds)
            elif len(bonds) != 6 and len(bonds) != 18 and composition == {'C'}:
                self._fix_nonaromatic(uniques)
            else:   # all bonds are aromatic
                self._set_atomatic(bonds)
            self._update_cache()

    def _set_atomatic(self, bonds):
        get_bond = lambda i: self.rwmol.GetBondWithIdx(i)
        for bond_i in bonds:
            bond = get_bond(bond_i)
            bond.SetBondType(Chem.BondType.AROMATIC)
            bond.GetBeginAtom().SetIsAromatic(True)
            bond.GetEndAtom().SetIsAromatic(True)

    def _fix_penta_aromatic(self, bonds):
        """
        The ring has an invalid composition for an aromatic ring.
        """
        self.log.debug(f'Fixing penta aromatic carbon ring "cyclopent-ine" to pyrrole')
        get_bond = lambda i: self.rwmol.GetBondWithIdx(i)
        atoms = [get_bond(bond_i).GetBeginAtom() for bond_i in bonds]
        # find the central atom
        for atom in atoms:  #: Chem.Atom
            if len(atom.GetNeighbors()) == 2:
                # downgrade to nitrogen
                atom.SetNumRadicalElectrons(0)
                atom.SetAtomicNum(7)
                atom.SetNumExplicitHs(1)
                return self._set_atomatic(bonds)

    def _fix_nonaromatic(self, bonds):
        self.log.debug(f'Downgrading non-aromatic ring to single bond')
        get_bond = lambda i: self.rwmol.GetBondWithIdx(i)
        for bond in map(get_bond, bonds):
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.GetBeginAtom().SetIsAromatic(False)
            bond.GetEndAtom().SetIsAromatic(False)

    def check_for_spiro(self, ring_idx: int, ringatoms: Tuple[Tuple[int]]) -> bool:
        ringatom_idxs = ringatoms[ring_idx]
        for atom_idx in ringatom_idxs:
            if len(self.rwmol.GetAtomWithIdx(atom_idx).GetNeighbors()) > 3:
                return True
        else:
            return False
