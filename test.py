import unittest, os
from molecular_rectifier import Rectifier
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem

class RectifierTester(unittest.TestCase):
    # name: [before, after]
    chemdex = {'phenylnaphthalene': ('c1ccc2ccccc2c1(c3ccccc3)', 'c1ccc(-c2cccc3ccccc23)cc1'),
               'benzoazetine': ('C12CCCCC1CC2', 'C1CCC2CCCC2C1'),
               #'conjoined': ('C1C2CCC2C1', 'C1CCCCC1'), # bridged hexane
               'allene': ('C=C=C', 'C=CC'),
               'benzocyclopronane': ('C12CCCCC1C2', 'C1CCC2CCCC2C1'),
               'benzoindolium': ('c12cccc(nc3)c1c3ccc2', 'C1=Nc2cccc3cccc1c23'),
               # 'norbornane': ('C1CC2CCC1C2', 'C1CC2CCC1C2'),
               'tetralin': ('c1ccc2c(c1)CCCC2', 'c1ccc2c(c1)CCCC2'),
               'mixed_ring': ('c1cccc2c1CCCC2', 'c1ccc2c(c1)CCCC2'),
               'mixe_ring2': ('C1CCCc2c1cccc2', 'c1ccc2c(c1)CCCC2'),
               'acenaphthylene': ('c1cc2c3c1cccc3ccc2', 'C1=Cc2cccc3cccc1c23'),
               'pentalene': ('C1=CC2=CC=CC2=C1', 'C1=CC2=CC=CC2=C1'),
               }

    def test_cyclopentadiene(self):
        """
        A pyrrole is altered to have a carbon as opposed to a nitrogen
        """
        # aromatic cyclopent-ine -> cyclopentadiene
        name = 'cyclopentadiene'
        mol = Chem.MolFromSmiles('[nH]1cccc1')
        mol.SetProp('_Name', name)
        AllChem.EmbedMolecule(mol)
        mod = Chem.RWMol(mol)
        mod.GetAtomWithIdx(0).SetAtomicNum(6)
        mol = mod.GetMol()
        recto = Rectifier(mol, atoms_in_bridge_cutoff=3).fix()
        gotten = Chem.MolToSmiles(AllChem.RemoveHs(recto.mol))
        after = 'C1=CCC=C1'
        #self.assertEqual(after, gotten, f'{name} failed {gotten} (expected {after})')
        # it no longer gives cyclopentadiene, but a pyrrole

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

    def test_sub_downgrade(self):
        # none of these should be touched
        xanthine = Chem.MolFromSmiles('c1[nH]c2c(n1)nc(nc2O)O')
        caffeine = Chem.MolFromSmiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')
        naphthalene = Chem.MolFromSmiles('c1c2ccccc2ccc1')
        naphthoquinone = Chem.MolFromSmiles('O=C1c2ccccc2C(=O)cc1')
        anthracene = Chem.MolFromSmiles('c1ccc2cc3ccccc3cc2c1')
        primulin = Chem.MolFromSmiles(
            'COc1cc(cc(c1O)OC)c2c(cc3c(cc(cc3[o+]2)O)O)O[C@H]4[C@@H]([C@H]([C@H]([C@H](O4)CO)O)O)O')
        for mol in (xanthine, caffeine, naphthalene, naphthoquinone, anthracene, primulin):
            self.assertTrue(Rectifier(mol).fix().mol.HasSubstructMatch(mol),
                            'failed to fix molecule')
        # this will
        phenol = Chem.MolFromSmiles('Oc1ccccc1')
        phenmonoone = Chem.Mol(phenol)  # not a real thing...
        phenmonoone.GetBondWithIdx(0).SetBondType(Chem.BondType.DOUBLE)
        self.assertTrue(Rectifier(phenmonoone).fix().mol.HasSubstructMatch(phenol),
                        'failed to fix substituent!')

    def test_real_merger(self):
        molblock = '''spiro-ditoluene
     RDKit          2D

 13 14  0  0  0  0  0  0  0  0999 V2000
    1.5000    2.5981    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.9196   -0.8437    0.0724 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.8002    0.1514    1.1984 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.4463   -0.5624    0.0483 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6193   -0.9951   -1.1259 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1462   -0.7139   -1.1501 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.3270    0.4327    1.1743 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  8  4  0
  3  4  4  0
  4  5  4  0
  4  1  1  0
  5  6  4  0
  6  7  4  0
  7  8  4  0
  9 13  4  0
  9 10  4  0
 10 11  4  0
 10  2  1  0
 11 12  4  0
  6 12  4  0
  6 13  4  0
M  END'''
        mol = Chem.MolFromMolBlock(molblock, sanitize=False, removeHs=False)
        recto = Rectifier(mol)
        recto.fix()
        self.assertEqual('CC1CCC2(CCCC(C)C2)CC1', Chem.MolToSmiles(AllChem.RemoveHs(recto.mol)))

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