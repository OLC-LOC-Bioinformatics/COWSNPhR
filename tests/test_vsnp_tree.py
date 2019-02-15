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

test_path = os.path.abspath(os.path.dirname(__file__))
file_path = os.path.join(test_path, 'files', 'vcf')
dependency_path = os.path.join(os.path.dirname(test_path), 'dependencies')
report_path = os.path.join(file_path, 'reports')
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
        assert VSNPTree(path=file_path)


def test_valid_path():
    global vcf_object
    vcf_object = VSNPTree(path=file_path,
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
        VSNPTreeMethods.file_list(file_path=test_path)


def test_vcf_file_list():
    global file_list
    file_list = VSNPTreeMethods.file_list(file_path=file_path)
    assert len(file_list) == 7


def test_strain_dict():
    global strain_vcf_dict
    strain_vcf_dict = VSNPTreeMethods.strain_list(vcf_files=file_list)
    assert [strain for strain in sorted(strain_vcf_dict)] == \
           ['13-1941', '13-1950', '13-1950_legacy', '13-1951',
            'B13-0234', 'B13-0235', 'B13-0238']


def test_accession_species():
    global accession_species_dict
    accession_species_dict = VSNPTreeMethods.parse_accession_species(ref_species_file=os.path.join(
        dependency_path, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_vcf_load():
    global strain_vcf_object_dict, strain_best_ref_dict
    strain_vcf_object_dict, strain_best_ref_dict = VSNPTreeMethods.load_vcf(strain_vcf_dict=strain_vcf_dict)
    for strain_name, best_ref in strain_best_ref_dict.items():
        assert best_ref in ['NC_002945.4', 'NC_017250.1', 'NC_017251.1']


def test_determine_ref_species():
    global strain_species_dict, strain_best_ref_fasta_dict
    strain_species_dict, strain_best_ref_fasta_dict = \
        VSNPTreeMethods.determine_ref_species(strain_best_ref_dict=strain_best_ref_dict,
                                              accession_species_dict=accession_species_dict)
    assert strain_species_dict['13-1941'] == 'af'
    assert strain_species_dict['B13-0234'] == 'suis1'
    assert strain_best_ref_fasta_dict['13-1941'] == 'NC_002945v4.fasta'


def test_reference_path():
    global reference_link_path_dict, reference_link_dict
    reference_link_path_dict, \
        reference_link_dict = VSNPTreeMethods.reference_folder(strain_best_ref_fasta_dict=strain_best_ref_fasta_dict,
                                                               dependency_path=dependency_path)
    assert reference_link_path_dict['13-1950'] == 'mycobacterium/tbc/af2122/script_dependents'


def test_extract_defining_snps():
    global defining_snp_dict
    defining_snp_dict = VSNPTreeMethods.extract_defining_snps(reference_link_path_dict=reference_link_path_dict,
                                                              strain_species_dict=strain_species_dict,
                                                              dependency_path=dependency_path)
    for species, snp_dict in defining_snp_dict.items():
        if species == 'af':
            assert snp_dict['Mbovis-01']['NC_002945.4'] == '2138896'
            assert snp_dict['Group_9_15']['NC_002945.4'] == '2145868!'
        elif species == 'suis1':
            assert snp_dict['Bsuis1-01']['NC_017251.1'] == '213522'


def test_load_snp_positions():
    global ref_snp_positions, strain_snp_positions, strain_snp_sequence
    ref_snp_positions, strain_snp_positions, strain_snp_sequence = \
        VSNPTreeMethods.load_snp_positions(strain_vcf_object_dict=strain_vcf_object_dict,
                                           strain_species_dict=strain_species_dict,
                                           strain_best_ref_dict=strain_best_ref_dict)
    assert ref_snp_positions['NC_002945.4'][4480] == 'T'
    print(ref_snp_positions['NC_017250.1'])
    assert ref_snp_positions['NC_017250.1'][1197913] == 'T'
    assert ref_snp_positions['NC_017251.1'][50509] == 'A'
    assert sorted(strain_snp_positions['13-1941'])[0] == 1057
    assert sorted(strain_snp_positions['B13-0234'])[0] == 8810
    assert strain_snp_sequence['13-1941'][1057] == 'G'
    assert strain_snp_sequence['13-1941'][2071827] == 'R'


def test_determine_groups():
    global strain_groups
    strain_groups = VSNPTreeMethods.determine_groups(strain_snp_positions, defining_snp_dict)
    assert strain_groups['13-1941'] == ['Mbovis-01', 'Mbovis-01A', 'Group_1-5', 'Group_9_15']
    assert strain_groups['B13-0234'] == ['Bsuis1-09', 'Bsuis1-09B']


def test_load_filter_file():
    global filter_dict
    filter_dict = VSNPTreeMethods.load_filter_file(reference_link_path_dict=reference_link_path_dict,
                                                   strain_best_ref_dict=strain_best_ref_dict,
                                                   dependency_path=dependency_path)
    assert filter_dict['NC_002945.4']['Mbovis-All'][0] == 79513
    assert filter_dict['NC_017250.1']['Bsuis1-All'][0] == 282708
    assert len(filter_dict['NC_017251.1']['Bsuis1-All']) == 427


def test_filter_positions():
    global strain_filtered_sequences, group_positions_dict
    strain_filtered_sequences, group_positions_dict = \
        VSNPTreeMethods.filter_positions(strain_snp_positions=strain_snp_positions,
                                                                 strain_groups=strain_groups,
                                                                 strain_best_ref_dict=strain_best_ref_dict,
                                                                 filter_dict=filter_dict,
                                                                 strain_snp_sequence=strain_snp_sequence)
    assert strain_filtered_sequences['13-1941']['Mbovis-All'][1057] == 'G'
    with pytest.raises(KeyError):
        assert strain_filtered_sequences['13-1941']['Bsuis1-All']
    assert strain_filtered_sequences['B13-0234']['Bsuis1-All'][29269] == 'G'
    with pytest.raises(KeyError):
        assert strain_filtered_sequences['B13-0234']['Bsuis1-All'][1057]
    assert 364560 in group_positions_dict['Mbovis-All']
    print(strain_filtered_sequences[''])
    assert False


def test_create_multifasta():
    global group_folders, species_folders,group_fasta_dict
    group_folders, species_folders, group_fasta_dict = VSNPTreeMethods.create_multifasta(strain_filtered_sequences, strain_species_dict, group_positions_dict, file_path)
    assert len(group_folders) == 9
    assert len(species_folders) == 2
    for group, fasta in group_fasta_dict.items():
        assert os.path.getsize(fasta) > 100


def test_load_genbank_file():
    global best_ref_gbk_dict
    best_ref_gbk_dict = VSNPTreeMethods.load_genbank_file(reference_link_path_dict=reference_link_path_dict,
                                                          strain_best_ref_dict=strain_best_ref_dict,
                                                          dependency_path=dependency_path)
    for key, value in best_ref_gbk_dict['NC_002945.4'].items():
        for feature in value.features:
            assert feature.type == 'source'
            assert str(feature.location) == '[0:4349904](+)'
            # Only test the first feature
            break


# def test_species_specific_working_dir():
#     global strain_specific_dir
#     strain_specific_dir = \
#         VSNPTreeMethods.species_specific_working_dir(strain_species_dict=strain_species_dict,
#                                                      reference_link_path_dict=reference_link_path_dict,
#                                                      path=file_path)


def test_vsnp_tree_run():
    vsnp_tree = VSNPTree(path=file_path,
                         threads=threads)
    vsnp_tree.main()


# def test_remove_species_folders():
#     for species_folder in species_folders:
#         shutil.rmtree(species_folder)