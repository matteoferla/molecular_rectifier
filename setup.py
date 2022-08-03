from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import sys, os

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ---------- Non pip modules  ------------------------------------------------------------------------------------------

if not util.find_spec('rdkit'):
    warn('This 3.6+ script **requires** rdkit which cannot be pip installed.' +
                              ' To install try either ' +
                              'conda install -c conda-forge rdkit or ' +
                              'sudo apt-get/brew install python3-rdkit or visit rdkit documentation.')

# ------------ description ---------------------------------------------------------------------------------------------

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    __doc__ = f.read()
descr = __doc__.split('\n')[1] # non-title line.

# ------------ Setup ---------------------------------------------------------------------------------------------------

setup(
    name='molecular_rectifier',
    version='0.1.8',
    python_requires='>3.6',
    packages=find_packages(),
    url='https://github.com/matteoferla/molecular_rectifier',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description=descr,
    long_description=__doc__,
    long_description_content_type='text/markdown',
)
