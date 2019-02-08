import os
from setuptools import setup

try:
    import ntpath
except:
    print("Seems you don't have ntpath installed!")  

reqs = ['numpy==1.15.0', 'h5py>=2.8.0', 'plotly==2.7.0', 'python_dateutil==2.7.3', 'biopython']

if os.name != 'nt':
    reqs = reqs + ['pysam==0.13']

setup(name='PyPore',
      version='0.1',
      description='A python scripts suite for the analysis and alignment of Nanopore data.',
      url='https://github.com/rsemeraro/PyPore',
      requires=['python (>=2.7, <3.0)'],
      author='Roberto Semeraro',
      author_email='robe.semeraro@gmail.com',
      license='LICENSE.txt',
      install_requires=reqs,
      package_data={'lib': ['*py'], 'lib.minimap2': ['*'], 'lib.SAMtools': ['*']},
      packages=['lib', 'lib.minimap2', 'lib.SAMtools'],
      include_package_data=True,
      entry_points={
          'console_scripts': [
              'pypore=lib.pypore:main'
          ],
       },         
)
