from setuptools import setup, find_packages
from warnings import warn
from importlib import util
import sys, os

if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print(sys.version_info)
    raise SystemError('Module written for Python 3.6+.')

# ------------ description ---------------------------------------------------------------------------------------------

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    __doc__ = f.read()
descr = __doc__.split('\n')[1] # non-title line.

# ------------ Setup ---------------------------------------------------------------------------------------------------

setup(
    name='molecular_rectifier',
    version='0.1.10.2',
    python_requires='>3.6',
    packages=find_packages(),
    install_requires=['rdkit'],
    url='https://github.com/matteoferla/molecular_rectifier',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@gmail.com',
    description=descr,
    long_description=__doc__,
    long_description_content_type='text/markdown',
)
