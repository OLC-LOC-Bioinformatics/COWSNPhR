#!/usr/bin/env python3
"""
Set up the package
"""

__author__ = 'adamkoziol'

# Standard inputs
from distutils.util import convert_path
import os

# Third party inputs
from setuptools import setup, find_packages

# Find the version
version = {}
with open(convert_path(os.path.join('cowsnphr_src', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)


setup(
    name='cowsnphr',
    version=version['__version__'],
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'cowsnphr.py = cowsnphr_src.cowsnphr:main',
            'cowsnphr = cowsnphr_src.cowsnphr:main',
        ],
    },
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@inspection.gc.ca',
    description='Single Nucleotide Variant Calling Pipeline',
    url='https://github.com/OLC-LOC-Bioinformatics/COWSNPhR',
    long_description=open('README.md').read(),
)
