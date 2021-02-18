from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import sys

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    raise ModuleNotFoundError('This 3.6+ script **requires** rdkit which cannot be pip installed.' +
                              ' To install try either ' +
                              'conda install -c conda-forge rdkit or ' +
                              'sudo apt-get/brew install python3-rdkit or visit rdkit documentation.')

setup(
    name='molecular_rectifier',
    version='0.1.1',
    packages=find_packages(),
    url='https://github.com/matteoferla/molecular_rectifier',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description='Given an RDKit molecule that does not sanitise, correct it until it does, '+
                'regardless of the severity of the change.'
)
