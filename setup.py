from setuptools import setup
from os import path

location = path.abspath(path.dirname(__file__))

with open(path.join(location, 'README.md')) as fh:
    longdesc = fh.read()

setup(name='seq-qc',
      version='0.16.0',
      description='utilities for performing various preprocessing steps on '
          'sequencing reads',
      long_description=longdesc,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      keywords='bioinformatics sequence preprocessing quality control',
      url='https://github.com/Brazelton-Lab/seq_qc/',
      author='Christopher Thornton',
      author_email='christopher.thornton@utah.edu',
      license='GPLv2',
      packages=['seq_qc',],
      include_package_data=True,
      zip_safe=False,
      install_requires=['screed',],
      entry_points={
          'console_scripts': [
              'qtrim = seq_qc.qtrim:main',
              'filter_replicates = seq_qc.filter_replicates:main',
              'interleave_pairs = seq_qc.interleave_pairs:main'
          ]
      }
      )
