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
    assert len(file_list) == 10


def test_strain_dict():
    global strain_vcf_dict
    strain_vcf_dict = VSNPTreeMethods.strain_list(vcf_files=file_list)
    assert [strain for strain in sorted(strain_vcf_dict)] == \
           ['03-1057', '13-1941', '13-1950', '13-3082', '14-2093', 'B13-0234', 'B13-0235', 'B13-0237', 'B13-0238',
            'B13-0239']


def test_accession_species():
    global accession_species_dict
    accession_species_dict = VSNPTreeMethods.parse_accession_species(ref_species_file=os.path.join(
        dependency_path, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_vcf_load():
    global strain_parsed_vcf_dict, strain_best_ref_dict
    strain_parsed_vcf_dict, strain_best_ref_dict = VSNPTreeMethods.load_vcf(strain_vcf_dict=strain_vcf_dict)
    for strain_name, best_ref in strain_best_ref_dict.items():
        assert best_ref in ['NC_002945.4', 'NC_017250.1', 'NC_017251.1']
    assert strain_parsed_vcf_dict['13-1941'][1057]['CHROM'] == 'NC_002945.4'
    assert strain_parsed_vcf_dict['13-1950'][29470]['FILTER'] == 'DELETION'
    assert strain_parsed_vcf_dict['13-1950'][1057]['STATS']['GT'] == '1/1'
    assert strain_parsed_vcf_dict['13-1950'][55517]['FILTER'] == 'INSERTION'
    assert strain_parsed_vcf_dict['13-1950'][714775]['ALT'] == 'AG'
    assert strain_parsed_vcf_dict['13-1950'][714775]['STATS']['VAF'] == '0.75,0'
    with pytest.raises(KeyError):
        assert strain_parsed_vcf_dict['13-1950_legacy'][55517]['FILTER'] == 'INSERTION'
    assert strain_parsed_vcf_dict['B13-0235'][8810]['QUAL'] == '70.1'
    assert strain_parsed_vcf_dict['B13-0235'][78796]['FILTER'] == 'DELETION'


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


def test_consolidate_group_ref_genomes():
    global strain_consolidated_ref_dict
    strain_consolidated_ref_dict = \
        VSNPTreeMethods.consolidate_group_ref_genomes(reference_link_dict=reference_link_dict,
                                                      strain_best_ref_dict=strain_best_ref_dict)
    assert strain_consolidated_ref_dict['13-1941'] == 'NC_002945v4'
    assert strain_consolidated_ref_dict['13-1950'] == 'NC_002945v4'
    assert strain_consolidated_ref_dict['B13-0234'] == 'NC_017251-NC_017250'
    assert strain_consolidated_ref_dict['B13-0238'] == 'NC_017251-NC_017250'


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
            assert snp_dict['Bsuis1-02']['NC_017250.1'] == '1173757'


def test_load_snp_positions():
    global ref_snp_positions, strain_snp_positions
    ref_snp_positions, strain_snp_positions = \
        VSNPTreeMethods.load_snp_positions(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
                                           strain_consolidated_ref_dict=strain_consolidated_ref_dict)
    assert ref_snp_positions['NC_002945v4'][4480] == 'T'
    assert ref_snp_positions['NC_017251-NC_017250'][1197913] == 'T'
    assert ref_snp_positions['NC_017251-NC_017250'][16405] == 'C'
    assert ref_snp_positions['NC_017251-NC_017250'][50509] == 'A'
    assert sorted(strain_snp_positions['13-1941'])[0] == 1057
    assert sorted(strain_snp_positions['B13-0234'])[0] == 8810
    assert sorted(strain_snp_positions['13-1950'])[1] == 4480


def test_determine_groups():
    global strain_groups
    strain_groups = VSNPTreeMethods.determine_groups(strain_snp_positions=strain_snp_positions,
                                                     defining_snp_dict=defining_snp_dict)
    assert strain_groups['13-1941'] == ['All', 'Mbovis-01', 'Mbovis-01A', 'Mbovis-01A3', 'Group_1-5', 'Group_9_15']
    assert strain_groups['13-1950'] == ['All', 'Mbovis-01', 'Mbovis-01A', 'Mbovis-01A3', 'Group_1-5', 'Group_9_15']
    assert strain_groups['B13-0234'] == ['All', 'Bsuis1-09', 'Bsuis1-09B']


def test_determine_group_snp_positions():
    global group_positions_set
    group_positions_set = VSNPTreeMethods.determine_group_snp_positions(strain_snp_positions=strain_snp_positions,
                                                                        strain_groups=strain_groups,
                                                                        strain_species_dict=strain_species_dict)
    assert len(group_positions_set['af']['All']) == 516
    assert len(group_positions_set['af']['Mbovis-01A3']) == 516
    assert len(group_positions_set['suis1']['Bsuis1-09']) == 108
    assert sorted(group_positions_set['suis1']['Bsuis1-09'])[0] == 8810


def test_load_snp_sequence():
    global group_strain_snp_sequence
    group_strain_snp_sequence = \
        VSNPTreeMethods.load_snp_sequence(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
                                          strain_consolidated_ref_dict=strain_consolidated_ref_dict,
                                          group_positions_set=group_positions_set,
                                          strain_groups=strain_groups,
                                          strain_species_dict=strain_species_dict,
                                          ref_snp_positions=ref_snp_positions)
    assert group_strain_snp_sequence['af']['All']['13-1941'][1057] == 'G'
    assert len(group_strain_snp_sequence['af']['All']['13-1941']) == 516
    assert group_strain_snp_sequence['af']['All']['13-1941'][948023] == 'M'
    assert group_strain_snp_sequence['suis1']['All']['B13-0234'][1195206] == 'R'
    assert group_strain_snp_sequence['suis1']['All']['B13-0234'][623503] == '-'


def test_create_multifasta():
    global group_folders, species_folders,group_fasta_dict, fasta_path
    fasta_path = os.path.join(file_path, 'alignments')
    group_folders, species_folders, group_fasta_dict = \
        VSNPTreeMethods.create_multifasta(group_strain_snp_sequence=group_strain_snp_sequence,
                                          fasta_path=fasta_path)
    assert len(group_folders) == 9
    assert len(species_folders) == 2
    for species, group_dict in group_fasta_dict.items():
        for group, fasta in group_dict.items():
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


def test_vsnp_tree_run():
    vsnp_tree = VSNPTree(path=file_path,
                         threads=threads)
    vsnp_tree.main()


def test_remove_species_folders():
    for species_folder in species_folders:
        shutil.rmtree(species_folder)


def test_remove_alignments_folder():
    shutil.rmtree(fasta_path)
