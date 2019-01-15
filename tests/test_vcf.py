#!/usr/bin/env python3
from vsnp.vcf import VCF
import pytest
import os
__author__ = 'adamkoziol'


def test_invalid_path():
    with pytest.raises(AssertionError):
        assert VCF(path='not_a_real_path')


def test_valid_path():
    assert VCF(path=os.getcwd())
