import unittest, os
from molecular_rectifier import Rectifier
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem

class RectifierTester(unittest.TestCase):
    # name: [before, after]
    chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
               'benzoazetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
               # 'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'), # bridged hexane
               'allene': ('C=C=C', 'C=CC'),
               'benzocyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
               'benzoindolium': ('c12cccc(nc3)c1c3ccc2', 'C1=Nc2cccc3cccc1c23'),
               # 'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2'),
               'tetralin': ('c1ccc2c(c1)CCCC2', 'c1ccc2c(c1)CCCC2'),
               'mixed_ring': ('c1cccc2c1CCCC2', 'c1ccc2c(c1)CCCC2'),
               'mixe_ring2': ('C1CCCc2c1cccc2', 'c1ccc2c(c1)CCCC2'),
               'acenaphthylene': ('c1cc2c3c1cccc3ccc2', 'c1cc2c3c1cccc3ccc2'),
               }

    def test_cyclopentine(self):
        """
        A pyrrole is altered to have a carbon as opposed to a nitrogen
        """
        # aromatic cyclopent-ine -> cyclopentadiene
        name = 'cyclopentine'
        mol = Chem.MolFromSmiles('[nH]1cccc1')
        mol.SetProp('_Name', name)
        AllChem.EmbedMolecule(mol)
        mod = Chem.RWMol(mol)
        mod.GetAtomWithIdx(0).SetAtomicNum(6)
        mol = mod.GetMol()
        recto = Rectifier(mol, atoms_in_bridge_cutoff=3).fix()
        gotten = Chem.MolToSmiles(AllChem.RemoveHs(recto.mol))
        after = 'C1=CCC=C1'
        self.assertEqual(after, gotten, f'{name} failed {gotten} (expected {after})')

    def test_bad_ring(self):
        name = 'bad ring'
        after = 'c1ccc2c(c1)CCCC2'
        mol = Chem.MolFromSmiles(after)
        mol.SetProp('_Name', name)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
        AllChem.EmbedMolecule(mol)
        before = Chem.MolToSmiles(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(AllChem.RemoveHs(recto.mol))
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_bad_ring2(self):
        name = 'bad ring2'
        before = 'c1ccc2c(c1)CCCC2'
        after = 'c1ccc2ccccc2c1'
        mol = Chem.MolFromSmiles(before)
        mol.SetProp('_Name', name)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
        mol.GetBondBetweenAtoms(6, 7).SetBondType(Chem.BondType.AROMATIC)
        before = Chem.MolToSmiles(mol)
        AllChem.EmbedMolecule(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(AllChem.RemoveHs(recto.mol))
        self.assertEqual(after, gotten,  f'{name} failed {gotten} (expected {after})')

    def test_emergency_overclose(self):
        mol: Chem.Mol = Chem.MolFromSmiles('CC(C)(C)')
        conf = Chem.Conformer()
        for i in range(4):
            conf.SetAtomPosition(i, Geometry.Point3D(i/4, 0, 0))
        mol.AddConformer(conf)
        recto = Rectifier(mol).fix()
        self.assertEqual(AllChem.RemoveHs(recto.mol).GetNumAtoms(), 1)

from functools import partial


def test_factory(name):
    def test(self):
        before, after = self.chemdex[name]
        mol = Chem.MolFromSmiles(before)
        mol.SetProp('_Name', name)
        AllChem.EmbedMolecule(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(AllChem.RemoveHs(recto.mol))
        self.assertEqual(after, gotten, f'{name} failed {gotten} (expected {after}) from {before}')
    return test

for name in RectifierTester.chemdex:
    setattr(RectifierTester,
            f'test_{name}',
            test_factory(name))