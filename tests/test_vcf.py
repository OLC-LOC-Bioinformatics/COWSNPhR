#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path
from vsnp.methods import Methods
from vsnp.vcf import VCF
import pytest
import os
__author__ = 'adamkoziol'

testpath = os.path.abspath(os.path.dirname(__file__))
filepath = os.path.join(testpath, 'files')


def test_invalid_path():
    with pytest.raises(AssertionError):
        assert VCF(path='not_a_real_path',
                   threads=1)


def test_no_threads():
    with pytest.raises(TypeError):
        assert VCF(path=testpath)


def test_valid_path():
    global vcf_object
    vcf_object = VCF(path=testpath,
                     threads=1)
    assert vcf_object


def test_empty_filer():
    fileset = filer(filelist=list())
    assert fileset == set()


def test_empty_filer_dict():
    filedict = filer(filelist=list(),
                     returndict=True)
    assert filedict == dict()


def test_normal_filer_dict():
    filedict = filer(filelist=['03-1057_S10_L001_R1_001.fastq.gz', '03-1057_S10_L001_R2_001.fastq.gz',
                               '13-1941_S4_L001_R1_001.fastq.gz', '13-1941_S4_L001_R2_001.fastq.gz'],
                     returndict=True)
    assert [file_name for file_name in filedict] == ['03-1057', '13-1941']


def test_normal_filer():
    fileset = filer(filelist=['03-1057_S10_L001_R1_001.fastq.gz', '03-1057_S10_L001_R2_001.fastq.gz',
                              '13-1941_S4_L001_R1_001.fastq.gz', '13-1941_S4_L001_R2_001.fastq.gz'])
    assert fileset == {'03-1057', '13-1941'}


def test_missing_file_filer():
    fileset = filer(filelist=['03-1057_S10_L001_R1_001.fastq.gz', '03-1057_S10_L001_R2_001.fastq.gz',
                              '13-1941_S4_L001_R1_001.fastq.gz'])
    assert fileset == {'03-1057', '13-1941'}


def test_non_illumina_filer():
    fileset = filer(filelist=['03-1057_R1.fastq.gz', '03-1057_R2.fastq.gz',
                              '13-1941_R1.fastq.gz', '13-1941_R2.fastq.gz'])
    assert fileset == {'03-1057', '13-1941'}


def test_multiple_differences_filer():
    fileset = filer(filelist=['03-1057_1.fastq', '03-1057_2.fastq',
                              '13-1941_1.fastq.gz', '13-1941_2.fastq'])
    assert fileset == {'03-1057', '13-1941'}


def test_no_directions_filer():
    fileset = filer(filelist=['03-1057.fastq.gz', '13-1941_S4_L001.fastq'])
    assert fileset == {'03-1057', '13-1941'}


def test_vcf_file_list_no_files():
    with pytest.raises(AssertionError):
        Methods.file_list(path=testpath)


def test_vcf_file_list():
    global file_list
    file_list = Methods.file_list(path=filepath)
    assert len(file_list) == 7


def test_strain_dict():
    global strain_folder_dict
    strain_folder_dict = Methods.strain_list(fastq_files=file_list)
    for strain_folder, fastq_files in strain_folder_dict.items():
        if os.path.basename(strain_folder) == '14-2093':
            assert len(fastq_files) == 1
        else:
            assert len(fastq_files) == 2


def test_strain_namer_no_input():
    strain_names = Methods.strain_namer(strain_folders=str())
    assert len(strain_names) == 0


def test_strain_namer_working():
    global strain_name_dict
    strain_name_dict = Methods.strain_namer(strain_folders=strain_folder_dict)
    assert [strain for strain in strain_name_dict] == ['03-1057', '13-1941', '13-1950', '14-2093']


def test_make_path():
    global make_path_folder
    make_path_folder = os.path.join(testpath, 'test_folder')
    make_path(make_path_folder)
    assert os.path.isdir(make_path_folder)


def test_rm_path():
    os.rmdir(make_path_folder)
    assert os.path.isdir(make_path_folder) is False


def test_strain_linker():
    global strain_fastq_dict
    strain_fastq_dict = Methods.file_link(strain_folder_dict=strain_folder_dict,
                                          strain_name_dict=strain_name_dict)
    assert [strain for strain in strain_fastq_dict] == ['03-1057', '13-1941', '13-1950', '14-2093']
    for strain_name, fastq_files in strain_fastq_dict.items():
        if strain_name == '14-2093':
            assert len(fastq_files) == 1
        else:
            assert len(fastq_files) == 2
        for symlink in fastq_files:
            assert os.path.islink(symlink)


def test_mash():
    Methods.call_mash(strain_fastq_dict=strain_fastq_dict)
