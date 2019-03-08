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
summary_path = os.path.join(file_path, 'summary_tables')
fasta_path = os.path.join(file_path, 'alignments')
report_path = os.path.join(file_path, 'reports')
logfile = os.path.join(file_path, 'log')
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
    global strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict
    strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict = \
        VSNPTreeMethods.load_vcf(strain_vcf_dict=strain_vcf_dict,
                                 threads=threads)
    for strain_name, best_ref in strain_best_ref_dict.items():
        assert best_ref in ['NC_002945.4', 'NC_017250.1', 'NC_017251.1']
    assert strain_best_ref_set_dict['13-1941'] == {'NC_002945.4'}
    assert strain_best_ref_set_dict['B13-0234'] == {'NC_017250.1', 'NC_017251.1'}
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


def test_summarise_vcf_outputs():
    pass_dict, insertion_dict, deletion_dict = \
        VSNPTreeMethods.summarise_vcf_outputs(strain_parsed_vcf_dict=strain_parsed_vcf_dict)
    assert pass_dict['13-1941'] == 493
    assert insertion_dict['B13-0234'] == 148
    assert deletion_dict['13-1950'] == 86911


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
    reference_link_path_dict, reference_link_dict = \
        VSNPTreeMethods.reference_folder(strain_best_ref_fasta_dict=strain_best_ref_fasta_dict,
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
    global consolidated_ref_snp_positions, strain_snp_positions, ref_snp_positions
    consolidated_ref_snp_positions, strain_snp_positions, ref_snp_positions = \
        VSNPTreeMethods.load_snp_positions(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
                                           strain_consolidated_ref_dict=strain_consolidated_ref_dict)
    assert consolidated_ref_snp_positions['NC_002945v4']['NC_002945.4'][4480] == 'T'
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017250.1'][1197913] == 'T'
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017250.1'][16405] == 'C'
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017251.1'][50509] == 'A'
    assert sorted(strain_snp_positions['13-1941']['NC_002945.4'])[0] == 1057
    assert sorted(strain_snp_positions['B13-0234']['NC_017251.1'])[0] == 29269
    assert sorted(strain_snp_positions['B13-0234']['NC_017250.1'])[0] == 8810
    with pytest.raises(KeyError):
        assert sorted(strain_snp_positions['B13-0234']['NC_002945.4'])[0] == 8810
    assert sorted(strain_snp_positions['13-1950']['NC_002945.4'])[1] == 4480
    assert ref_snp_positions['NC_002945.4'][4480] == 'T'
    assert ref_snp_positions['NC_017250.1'][1197913] == 'T'
    with pytest.raises(KeyError):
        assert ref_snp_positions['NC_017251.1'][16405]
    assert ref_snp_positions['NC_017250.1'][16405] == 'C'
    assert ref_snp_positions['NC_017251.1'][50509] == 'A'


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
    assert len(group_positions_set['af']['All']['NC_002945.4']) == 516
    assert len(group_positions_set['af']['Mbovis-01A3']['NC_002945.4']) == 516
    assert len(group_positions_set['suis1']['Bsuis1-09']['NC_017251.1']) == 63
    assert len(group_positions_set['suis1']['Bsuis1-09']['NC_017250.1']) == 45
    assert sorted(group_positions_set['suis1']['Bsuis1-09']['NC_017250.1'])[0] == 8810
    with pytest.raises(AssertionError):
        assert sorted(group_positions_set['suis1']['Bsuis1-09']['NC_017251.1'])[0] == 8810


def test_load_snp_sequence():
    global group_strain_snp_sequence, species_group_best_ref
    group_strain_snp_sequence, species_group_best_ref = \
        VSNPTreeMethods.load_snp_sequence(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
                                          strain_consolidated_ref_dict=strain_consolidated_ref_dict,
                                          group_positions_set=group_positions_set,
                                          strain_groups=strain_groups,
                                          strain_species_dict=strain_species_dict,
                                          consolidated_ref_snp_positions=consolidated_ref_snp_positions)
    assert group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][1057] == 'G'
    assert len(group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4']) == 516
    assert group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][948023] == 'M'
    with pytest.raises(KeyError):
        assert group_strain_snp_sequence['af']['All']['13-1941']['NC_017250.1']
    assert group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017250.1'][1195206] == 'R'
    assert group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017250.1'][623503] == '-'
    with pytest.raises(KeyError):
        assert group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017251.1'][623503]
    assert species_group_best_ref['af']['All'] == 'NC_002945v4'
    assert species_group_best_ref['suis1']['All'] == 'NC_017251-NC_017250'


def test_create_multifasta():
    global group_folders, species_folders, group_fasta_dict
    group_folders, species_folders, group_fasta_dict = \
        VSNPTreeMethods.create_multifasta(group_strain_snp_sequence=group_strain_snp_sequence,
                                          fasta_path=fasta_path)
    assert len(group_folders) == 9
    assert len(species_folders) == 2
    for species, group_dict in group_fasta_dict.items():
        for group, fasta in group_dict.items():
            assert os.path.getsize(fasta) > 100


def test_run_raxml():
    global species_group_trees
    species_group_trees = VSNPTreeMethods.run_raxml(group_fasta_dict=group_fasta_dict,
                                                    strain_consolidated_ref_dict=strain_consolidated_ref_dict,
                                                    strain_groups=strain_groups,
                                                    threads=threads,
                                                    logfile=logfile)
    for species, group_dict in species_group_trees.items():
        for group, options_dict in group_dict.items():
            assert os.path.getsize(options_dict['best_tree']) > 100
            assert os.path.getsize(options_dict['bootstrap_tree']) > 100


def test_parse_tree_order():
    global species_group_order_dict
    species_group_order_dict = VSNPTreeMethods.parse_tree_order(species_group_trees=species_group_trees)
    assert species_group_order_dict['af']['All'] == ['13-1941', '13-1950', '14-2093', '03-1057', '13-3082',
                                                     'NC_002945v4']

    assert species_group_order_dict['suis1']['All'] == ['B13-0235', 'B13-0238', 'B13-0234', 'B13-0237', 'B13-0239',
                                                        'NC_017251-NC_017250']


def test_copy_trees():
    global tree_path
    tree_path = os.path.join(file_path, 'tree_files')
    VSNPTreeMethods.copy_trees(species_group_trees=species_group_trees,
                               tree_path=tree_path)
    assert len(glob(os.path.join(tree_path, '*_All*'))) == 4
    assert len(glob(os.path.join(tree_path, '*Mbovis-01*'))) == 6
    assert len(glob(os.path.join(tree_path, '*suis1*'))) == 6
    assert len(glob(os.path.join(tree_path, '*'))) == 18


def test_load_genbank_file():
    global full_best_ref_gbk_dict
    full_best_ref_gbk_dict = VSNPTreeMethods.load_genbank_file(reference_link_path_dict=reference_link_path_dict,
                                                               strain_best_ref_set_dict=strain_best_ref_set_dict,
                                                               dependency_path=dependency_path)
    assert full_best_ref_gbk_dict['NC_002945.4'][0].type == 'CDS'
    assert full_best_ref_gbk_dict['NC_017251.1'][28315].type == 'tRNA'
    assert full_best_ref_gbk_dict['NC_017250.1'][1104133].qualifiers['gene'] == ['rrf']


def test_annotate_snps():
    global species_group_annotated_snps_dict
    species_group_annotated_snps_dict = \
        VSNPTreeMethods.annotate_snps(group_strain_snp_sequence=group_strain_snp_sequence,
                                      full_best_ref_gbk_dict=full_best_ref_gbk_dict,
                                      strain_best_ref_set_dict=strain_best_ref_set_dict,
                                      ref_snp_positions=ref_snp_positions)
    assert species_group_annotated_snps_dict['af']['Mbovis-01']['NC_002945.4'][1057]['locus'] == 'BQ2027_MB0001'
    assert species_group_annotated_snps_dict['af']['Mbovis-01A3']['NC_002945.4'][4480]['product'] == \
        'hypothetical protein'
    assert species_group_annotated_snps_dict['suis1']['All']['NC_017250.1'][579895]['gene'] == 'gatA'
    assert species_group_annotated_snps_dict['suis1']['All']['NC_017251.1'][633244]['product'] == 'porin'


def test_determine_snp_number():
    global species_group_snp_num_dict
    species_group_snp_num_dict = \
        VSNPTreeMethods.determine_snp_number(group_strain_snp_sequence=group_strain_snp_sequence,
                                             species_group_best_ref=species_group_best_ref)
    assert species_group_snp_num_dict['af']['All']['NC_002945.4'][1057] == 5
    assert species_group_snp_num_dict['suis1']['All']['NC_017250.1'][8810] == 5
    with pytest.raises(KeyError):
        assert species_group_snp_num_dict['suis1']['All']['NC_017251.1'][8810]


def test_rank_snps():
    global species_group_snp_rank, species_group_num_snps
    species_group_snp_rank, species_group_num_snps \
        = VSNPTreeMethods.rank_snps(species_group_snp_num_dict=species_group_snp_num_dict)
    assert species_group_snp_rank['af']['All'][5]['NC_002945.4'][0] == 1057
    assert species_group_snp_rank['suis1']['All'][5]['NC_017250.1'][0] == 8810
    with pytest.raises(AssertionError):
        assert species_group_snp_rank['suis1']['All'][5]['NC_017251.1'][0] == 8810
    assert len(species_group_snp_rank['af']['All'][5]['NC_002945.4']) == 473
    assert len(species_group_snp_rank['suis1']['All'][5]['NC_017251.1']) == 46
    assert len(species_group_snp_rank['suis1']['All'][1]['NC_017251.1']) == 17
    assert len(species_group_snp_rank['suis1']['All'][5]['NC_017250.1']) == 21
    assert species_group_num_snps['af']['All'] == 516
    assert species_group_num_snps['suis1']['All'] == 108


def test_sort_snps():
    global species_group_sorted_snps
    species_group_sorted_snps = \
        VSNPTreeMethods.sort_snps(species_group_order_dict=species_group_order_dict,
                                  species_group_snp_rank=species_group_snp_rank,
                                  species_group_best_ref=species_group_best_ref,
                                  group_strain_snp_sequence=group_strain_snp_sequence)
    assert species_group_sorted_snps['af']['All']['NC_002945.4'][0] == 1057
    assert species_group_sorted_snps['af']['All']['NC_002945.4'][-1] == 2523406
    assert species_group_sorted_snps['suis1']['All']['NC_017250.1'][0] == 8810
    assert species_group_sorted_snps['suis1']['All']['NC_017250.1'][-1] == 1183517
    assert species_group_sorted_snps['suis1']['All']['NC_017251.1'][-1] == 2065515
    assert len(species_group_sorted_snps['af']['All']['NC_002945.4']) == 516
    assert len(species_group_sorted_snps['suis1']['All']['NC_017251.1']) == 63
    assert len(species_group_sorted_snps['suis1']['All']['NC_017250.1']) == 45


def test_create_summary_table():
    VSNPTreeMethods.create_summary_table(species_group_sorted_snps=species_group_sorted_snps,
                                         species_group_order_dict=species_group_order_dict,
                                         species_group_best_ref=species_group_best_ref,
                                         group_strain_snp_sequence=group_strain_snp_sequence,
                                         species_group_annotated_snps_dict=species_group_annotated_snps_dict,
                                         species_group_num_snps=species_group_num_snps,
                                         summary_path=summary_path)
    assert len(glob(os.path.join(summary_path, '*.xlsx'))) == 9


def test_vsnp_tree_run():
    vsnp_tree = VSNPTree(path=file_path,
                         threads=threads)
    vsnp_tree.main()


def test_remove_species_folders():
    for species_folder in species_folders:
        shutil.rmtree(species_folder)


def test_remove_alignments_folder():
    shutil.rmtree(fasta_path)


def test_remove_tree_folder():
    shutil.rmtree(tree_path)


def test_remove_report_folder():
    shutil.rmtree(summary_path)


def test_remove_logs():
    logs = glob(os.path.join(file_path, '*.txt'))
    for log in logs:
        os.remove(log)
