from ._base import _RectifierBase
import numpy as np
from rdkit import Chem

from typing import List, Union
from collections import Counter

# ======================================================================================================================


class _RectifierOverclose(_RectifierBase):
    """
    Fragmenstein monster is better equipped to deal with overclose atoms
    This simply absorbs all which are egrigiously too close
    It is an extreme fallback
    """
    def absorb_overclose(self, distance_cutoff: float):
        """
        Never say never... This really should have an effect as Fragmenstein monster is build for this.
        But stuff always goes wrong
        """
        if self.rwmol.GetNumConformers() == 0:
            return
        # make matrix
        distance_matrix: np.ndarray = Chem.Get3DDistanceMatrix(self.rwmol)
        # nan fill the self values
        np.fill_diagonal(distance_matrix, np.nan)
        # nan fill the upper triangle
        distance_matrix[np.triu_indices(distance_matrix.shape[0],)] = np.nan
        # get overclose indices as a list of two indices
        mergers: List[List[int]] = np.transpose( (distance_matrix < distance_cutoff).nonzero() ).tolist()  # noqa
        if len(mergers) == 0:
            return
        # sort the mergers by most common atoms
        commonality = dict(Counter(np.array(mergers).flatten()).most_common())
        mergers = sorted(mergers, key=lambda idxs: commonality[idxs[0]]+commonality[idxs[1]])
        # sort each pair by most common (and then bonded) first
        sort_indices = lambda i: self._get_atom_valence(self.rwmol.GetAtomWithIdx(i)) + 10 * commonality[i]
        mergers = [sorted(idxs, key=sort_indices, reverse=True) for idxs in mergers]
        self.log.critical(f'Absorption required for the atom pairs: {mergers}')
        # mark for deletion
        for _, i in mergers:
            self.rwmol.GetAtomWithIdx(i).SetBoolProp('DELETE', True)
        # copy bonds
        for keep_i, del_i in mergers:
            for neigh in self.rwmol.GetAtomWithIdx(del_i).GetNeighbors():  #: Chem.Atom
                neigh_i: int = neigh.GetIdx()
                if neigh_i == keep_i:
                    continue  # I expect keep is del's neigh
                del_bondtype: Chem.BondType =  self.rwmol.GetBondBetweenAtoms(neigh_i, del_i).GetBondType()
                previous_bond: Union[None, Chem.Bond] = self.rwmol.GetBondBetweenAtoms(keep_i, neigh_i)
                if not previous_bond:
                    self.rwmol.AddBond(keep_i, neigh_i, del_bondtype)
        # delete atoms
        while True:
            for atom in self.rwmol.GetAtoms():
                if atom.HasProp('DELETE'):
                    self.rwmol.RemoveAtom(atom.GetIdx())
                    break
            else:
                break
        # store copy (self.mol is property <= .rwmol
        self.modifications.append(self.mol)

