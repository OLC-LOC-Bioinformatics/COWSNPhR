#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path
from vsnp.vsnp_tree_methods import VSNPTreeMethods
from vsnp.vsnp_tree_run import VSNPTree
from datetime import datetime
import multiprocessing
from glob import glob
import pytest
import shutil
import os

__author__ = 'adamkoziol'

testpath = os.path.abspath(os.path.dirname(__file__))
filepath = os.path.join(testpath, 'files', 'vcf')
dependencypath = os.path.join(os.path.dirname(testpath), 'dependencies')
report_path = os.path.join(filepath, 'reports')
threads = multiprocessing.cpu_count() - 1
# Define the start time
start_time = datetime.now()


def test_import_vsnp():
    with pytest.raises(SystemExit):
        import vsnp.vSNP


def test_invalid_path():
    with pytest.raises(AssertionError):
        assert VSNPTree(path='not_a_real_path',
                        threads=threads)


def test_no_threads():
    with pytest.raises(TypeError):
        assert VSNPTree(path=filepath)


def test_valid_path():
    global vcf_object
    vcf_object = VSNPTree(path=filepath,
                          threads=threads)
    assert vcf_object


def test_invalid_tilde_path():
    VSNPTree(path='~',
             threads=threads)


def test_empty_filer():
    fileset = filer(filelist=list())
    assert fileset == set()


def test_empty_filer_dict():
    filedict = filer(filelist=list(),
                     returndict=True)
    assert filedict == dict()


def test_normal_filer_dict():
    filedict = filer(filelist=['03-1057.vcf', '13-1941.vcf'],
                     extension='vcf',
                     returndict=True)
    assert [file_name for file_name in filedict] == ['03-1057', '13-1941']


def test_vcf_file_list_no_files():
    with pytest.raises(AssertionError):
        VSNPTreeMethods.file_list(path=testpath)


def test_vcf_file_list():
    global file_list
    file_list = VSNPTreeMethods.file_list(path=filepath)
    assert len(file_list) == 7


def test_strain_dict():
    global strain_folder_dict
    strain_folder_dict = VSNPTreeMethods.strain_list(vcf_files=file_list)
    assert len(strain_folder_dict) == 7


def test_strain_namer_no_input():
    strain_names = VSNPTreeMethods.strain_namer(strain_folders=str())
    assert len(strain_names) == 0


def test_strain_namer_working():
    global strain_name_dict
    strain_name_dict = VSNPTreeMethods.strain_namer(strain_folders=strain_folder_dict)
    assert [strain for strain in sorted(strain_name_dict)] == \
           ['13-1941', '13-1950', '13-1950_legacy', '13-1951',
            'B13-0234', 'B13-0235', 'B13-0238']


def test_make_path():
    global make_path_folder
    make_path_folder = os.path.join(testpath, 'test_folder')
    make_path(make_path_folder)
    assert os.path.isdir(make_path_folder)


def test_rm_path():
    os.rmdir(make_path_folder)
    assert os.path.isdir(make_path_folder) is False


def test_strain_linker():
    global strain_vcf_dict
    strain_vcf_dict = VSNPTreeMethods.file_link(strain_folder_dict=strain_folder_dict,
                                                strain_name_dict=strain_name_dict)
    assert [strain for strain in sorted(strain_vcf_dict)] == \
           ['13-1941', '13-1950', '13-1950_legacy', '13-1951',
            'B13-0234', 'B13-0235', 'B13-0238']


def test_symlink():
    for strain_name, vcf_link in strain_vcf_dict.items():
        assert os.path.islink(vcf_link)


def test_accession_species():
    global accession_species_dict
    accession_species_dict = VSNPTreeMethods.parse_accession_species(ref_species_file=os.path.join(
        dependencypath, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_vcf_load():
    global strain_vcf_object_dict, strain_best_ref_dict
    strain_vcf_object_dict, strain_best_ref_dict = VSNPTreeMethods.load_vcf(strain_vcf_dict=strain_vcf_dict)
    for strain_name, best_ref in strain_best_ref_dict.items():
        assert best_ref in ['NC_002945.4', 'NC_017250.1', 'NC_017251.1']


def test_determine_ref_species():
    global strain_species_dict, strain_best_ref_dict
    strain_species_dict, strain_best_ref_dict = VSNPTreeMethods.determine_ref_species(strain_best_ref_dict=strain_best_ref_dict,
                                                                accession_species_dict=accession_species_dict)
    assert strain_species_dict['13-1941'] == 'af'
    assert strain_species_dict['B13-0234'] == 'suis1'
    assert strain_best_ref_dict['13-1941'] == 'NC_002945v4.fasta'


def test_reference_path():
    global reference_link_path_dict, reference_link_dict
    reference_link_path_dict, \
        reference_link_dict = VSNPTreeMethods.reference_folder(strain_best_ref_dict=strain_best_ref_dict,
                                                               dependency_path=dependencypath)
    assert reference_link_path_dict['13-1950'] == 'mycobacterium/tbc/af2122/script_dependents/NC_002945v4.fasta'


def test_vsn_tree_run():
    vsnp_tree = VSNPTree(path=filepath,
                         threads=threads)
    vsnp_tree.main()


def test_remove_working_dir():
    for strain_name, working_dir in strain_name_dict.items():
        shutil.rmtree(working_dir)
