# Molecular Rectifier
Given an RDKit molecule that does not sanitise, correct it until it does, regardless of the severity of the change.

![horror](horror.png)

## Install

Requires RDKit.

    pip3 install molecular-rectifier

## Beyond RDKit Sanitisation

The command `rdkit.Chem.SanitizeMol` fixes minor issues with the molecule.
However, more drastic changes such as valence correction and removal of weird bonds is not done,
hence the molecular rectifier!

    from molecular_rectifier import Rectifier
    
    recto = Rectifier(wrong_mol)
    recto.fix()
    fixed_mol = recto.mol
    # this works:
    Chem.SanitizeMol(fixed_mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL)

A lot of the corrections are based on arbitrary choices. These choices are manually encoded.
I am sure that sooner or later some ML module will do this, but for now the Rectifier does the job.

The attribute `valence_correction` controls how the valence is fixed up
—either by shifting the element (`'element'`) or by adding a charge (`'charge'`)

* protonating ring nitrogens if required
* forcing all atoms/bonds in a ring to be aromatic or not (unless part of another ring)
* Texas carbon -> Sulfur
* Hydrogen -> Fluoride shifted downwards
* Carbon in aromatic ring -> single bond ring
* kills bridging bonds with no atoms in the bridge within rings.
* correcting `BondType.UNSPECIFIED`
* preventing 3- and 4- membered rings fused to another ring
* preventing allene

## Rationale

This is used by [Fragmenstein](https://github.com/matteoferla/Fragmenstein) in the automatic merging mode (`combine`).
Originally part of it, but moved apart because it may be useful for other uses, outside of Fragmenstein.
Or a better module is found to fix the molecules.