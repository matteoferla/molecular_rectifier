## Rules

### Minimum bond length

The method `.absorb_overclose` (from [_overclose.py](molecular_rectifier/_overclose.py)) has a default length of 0.8.

The atom to be absorbed is determined by the line

```python
sort_indices = lambda i: self._get_atom_valence(self.rwmol.GetAtomWithIdx(i)) + 10 * commonality[i]
```

where the ad-hoc 10x-weighted penalty `commonality` is the number of times the atom violates bond lengths.


NB. This happens before the method call `has_issues`.

### Rings

The method `._fix_aromatic_rings` (from (_ring.py)[molecular_rectifier/_ring.py]) does a bunch of checks,
if issues are present.
If the ring is a spiro ring, or a carbocyclic ring without 6 or 18 rings, the ring is reduced.
With a special case of 5-atom rings
(`_fix_penta_aromatic`) wherein a carbon is made into a nitrogen (ie. a pyrrole).

The method `.fix_rings` (from [_ring.py](molecular_rectifier/_ring.py)) is called regardlessly of issues
and enforces several rules that prevent unusual structures, not all illegal:

* The method `._prevent_bonded_to_bridgeheads` prevents a carbon forming 3-aromatic-bond _and_ a bond to something else.
* The method `.fix_aromatic_radicals` fixes radicals due to aromaticity.
* The method `._prevent_conjoined_ring` prevents 3- and 4- membered rings fused to another ring.
* The method `._prevent_weird_rings` prevents benzo-cyclopropane or benzo-azetine kind of molecules (legal, but weird).

The method `triage_rings` (from [_valence.py](molecular_rectifier/_valence.py)) fixes mixed bond issues 

### Allene

An allene is normal... but weird (`_prevent_allene` from [_odd.py](molecular_rectifier/_odd.py) correct it).

### Unspecified bonds

An unspecified bond (not racemic, but with unknown bond order) is corrected to a single bond
(`.ununspecified_bonds` see [_valence.py](molecular_rectifier/_valence.py)).

### Radical protons

The method `_adjust_Hs` (from [_valence.py](molecular_rectifier/_valence.py)) fixes radical protons.

### Valence

This is a nightmare. The method `.fix_issues` (from [_valence.py](molecular_rectifier/_valence.py))
iterates around the violating atoms and tries to find a workable protonation state.