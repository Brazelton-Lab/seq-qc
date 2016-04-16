from setuptools import setup, Extension

bernoulli = Extension('bernoulli', sources=['seq_qc/bernoullimodule.c'])

setup(name='seq-qc',
      version='1.1.3',
      packages=['seq_qc',],
      description='utilities for performing various preprocessing steps on '
          'sequencing reads',
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
      download_url = 'https://github.com/fBrazelton-Lab/seq_qc/tarball/v1.1.0',
      author='Christopher Thornton',
      author_email='christopher.thornton@utah.edu',
      license='GPLv2',
      include_package_data=True,
      zip_safe=False,
      install_requires=['screed',],
      ext_modules=[bernoulli,],
      entry_points={
          'console_scripts': [
              'qtrim = seq_qc.qtrim:main',
              'filter_replicates = seq_qc.filter_replicates:main',
              'error_filter = seq_qc.error_filter:main'
          ]
      }
      )
