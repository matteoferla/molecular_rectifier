import unittest, os
from molecular_rectifier import Rectifier
from rdkit import Chem
from rdkit.Chem import AllChem

class RectifierTester(unittest.TestCase):
    def test_rectifier(self):
        # name: [before, after]
        chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
                   'benzo-azetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
                   # 'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'), # bridged hexane
                   'allene': ('C=C=C', 'C=CC'),
                   'benzo-cyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
                   # 'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2'),
                   'mixed ring': ('c1cccc2c1CCCC2', 'c1ccc2c(c1)CCCC2'),
                   'mixed ring': ('C1CCCc2c1cccc2', 'c1ccc2c(c1)CCCC2'),
                   }

        for name in chemdex:
            before, after = chemdex[name]
            mol = Chem.MolFromSmiles(before)
            mol.SetProp('_Name', name)
            AllChem.EmbedMolecule(mol)
            recto = Rectifier(mol).fix()
            gotten = Chem.MolToSmiles(recto.mol)
            self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after}) from {before}')

    def test_cyclopentine(self):
        # aromatic cyclopent-ine -> cyclopentadiene
        name = 'cyclopentine'
        mol = Chem.MolFromSmiles('[nH]1cccc1')
        mol.SetProp('_Name', name)
        AllChem.EmbedMolecule(mol)
        mod = Chem.RWMol(mol)
        mod.GetAtomWithIdx(0).SetAtomicNum(6)
        mol = mod.GetMol()
        recto = Rectifier(mol, atoms_in_bridge_cutoff=3).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        after = 'C1=CCC=C1'
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')

    def test_bad_ring(self):
        name = 'bad ring'
        after = 'c1ccc2c(c1)CCCC2'
        mol = Chem.MolFromSmiles(after)
        mol.SetProp('_Name', name)
        mol.GetBondBetweenAtoms(0, 1).SetBondType(Chem.BondType.SINGLE)
        before = Chem.MolToSmiles(mol)
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(recto.mol)
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
        recto = Rectifier(mol).fix()
        gotten = Chem.MolToSmiles(recto.mol)
        self.assertEqual(gotten, after, f'{name} failed {gotten} (expected {after})')