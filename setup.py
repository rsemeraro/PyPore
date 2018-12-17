import os
from setuptools import setup

reqs = ['mpi4py==3.0.0', 'numpy==1.15.0', 'h5py==2.8.0', 'h5py_cache==1.0', 'plotly==2.7.0', 'python_dateutil==2.7.3', 'ntpath', 'biopython']

if os.name != 'nt':
    reqs = reqs + ['pysam==0.13']

setup(name='PyPore',
      version='0.1',
      description='A python scripts suite for the analysis and alignment of Nanopore data.',
      url='https://github.com/rsemeraro/PyPore',
      author='Roberto Semeraro',
      author_email='robe.semeraro@gmail.com',
      license='LICENSE.txt',
      scripts=['pypore'],
      packages=['lib'],
      install_requires=reqs
)