#!/usr/bin/env python
from setuptools import setup, find_packages
import os
__author__ = 'stuber, adamkoziol'

setup(
    name="cowsnphr",
    version="0.0.30",
    packages=find_packages(),
    include_package_data=True,
    scripts=[os.path.join('cowsnphr', 'cowsnphr.py')],
    license='MIT',
    author='Adam Koziol',
    author_email='adam.koziol@canada.ca',
    description='Single Nucleotide Variant Calling Pipeline',
    url='https://github.com/OLC-LOC-Bioinformatics/COWSNPhR',
    long_description=open('README.md').read(),
)
