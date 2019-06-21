#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import filer, make_path
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


def test_invalid_path():
    with pytest.raises(AssertionError):
        assert VSNPTree(path='not_a_real_path',
                        threads=threads,
                        debug=False,
                        variant_caller='deepvariant',
                        filter_positions=False)


def test_no_threads():
    with pytest.raises(TypeError):
        assert VSNPTree(path=file_path,
                        debug=False,
                        variant_caller='deepvariant',
                        filter_positions=False)


def test_valid_path():
    global vcf_object
    vcf_object = VSNPTree(path=file_path,
                          threads=threads,
                          debug=False,
                          variant_caller='deepvariant',
                          filter_positions=False)
    assert vcf_object


def test_invalid_tilde_path():
    VSNPTree(path='~',
             threads=threads,
             debug=False,
             variant_caller='deepvariant',
             filter_positions=False)


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


def test_gvcf_load():
    global gvcf_parsed_dict, gvcf_best_ref_dict, gvcf_best_ref_set_dict
    gvcf_vcf_dict = dict()
    for strain_name, vcf_file in strain_vcf_dict.items():
        if strain_name != '13-1950':
            gvcf_vcf_dict[strain_name] = vcf_file
    gvcf_parsed_dict, gvcf_best_ref_dict, gvcf_best_ref_set_dict = \
        VSNPTreeMethods.load_gvcf(strain_vcf_dict=gvcf_vcf_dict,
                                  threads=threads,
                                  qual_cutoff=20)
    for strain_name, best_ref in gvcf_best_ref_dict.items():
        assert best_ref in ['NC_002945.4', 'NC_017250.1', 'NC_017251.1']
    assert gvcf_best_ref_set_dict['13-1941'] == {'NC_002945.4'}
    assert gvcf_best_ref_set_dict['B13-0234'] == {'NC_017250.1', 'NC_017251.1'}
    assert gvcf_parsed_dict['13-1941']['NC_002945.4'][1057]['CHROM'] == 'NC_002945.4'
    with pytest.raises(KeyError):
        assert gvcf_parsed_dict['13-1950_legacy'][55517]['FILTER'] == 'INSERTION'
    assert gvcf_parsed_dict['B13-0235']['NC_017250.1'][8810]['QUAL'] == '70.1'
    assert gvcf_parsed_dict['B13-0235']['NC_017251.1'][78796]['FILTER'] == 'DELETION'


def test_vcf_load():
    global vcf_parsed_dict, vcf_best_ref_dict, vcf_best_ref_set_dict
    vcf_vcf_dict = dict()
    for strain_name, vcf_file in strain_vcf_dict.items():
        if strain_name == '13-1950':
            vcf_vcf_dict[strain_name] = vcf_file
    vcf_parsed_dict, vcf_best_ref_dict, vcf_best_ref_set_dict = \
        VSNPTreeMethods.load_vcf(strain_vcf_dict=vcf_vcf_dict)
    for strain_name, best_ref in vcf_best_ref_dict.items():
        assert best_ref in ['NC_002945.4']
    assert vcf_parsed_dict['13-1950']['NC_002945.4'][29470]['FILTER'] == 'DELETION'
    assert vcf_best_ref_set_dict['13-1950'] == {'NC_002945.4'}
    for strain_name, best_ref in vcf_best_ref_dict.items():
        assert best_ref in ['NC_002945.4']
    assert vcf_parsed_dict['13-1950']['NC_002945.4'][55517]['FILTER'] == 'INSERTION'
    assert vcf_parsed_dict['13-1950']['NC_002945.4'][714775]['ALT'] == 'A'
    # assert vcf_parsed_dict['13-1950']['NC_002945.4'][714775]['STATS']['VAF'] == '0.75,0'


def test_strain_parsed_vcf_dict_consolidate():
    global strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict
    strain_parsed_vcf_dict = dict()
    strain_best_ref_dict = dict()
    strain_best_ref_set_dict = dict()
    for strain_name, vcf_dict in gvcf_parsed_dict.items():
        best_ref = gvcf_best_ref_dict[strain_name]
        best_ref_set = gvcf_best_ref_set_dict[strain_name]
        strain_parsed_vcf_dict[strain_name] = vcf_dict
        strain_best_ref_dict[strain_name] = best_ref
        strain_best_ref_set_dict[strain_name] = best_ref_set
    for strain_name, vcf_dict in vcf_parsed_dict.items():
        best_ref = vcf_best_ref_dict[strain_name]
        best_ref_set = vcf_best_ref_set_dict[strain_name]
        strain_parsed_vcf_dict[strain_name] = vcf_dict
        strain_best_ref_dict[strain_name] = best_ref
        strain_best_ref_set_dict[strain_name] = best_ref_set
    assert strain_best_ref_set_dict['13-1941'] == {'NC_002945.4'}
    assert strain_best_ref_set_dict['B13-0234'] == {'NC_017250.1', 'NC_017251.1'}
    assert strain_parsed_vcf_dict['13-1941']['NC_002945.4'][1057]['CHROM'] == 'NC_002945.4'
    assert strain_parsed_vcf_dict['13-1950']['NC_002945.4'][55517]['FILTER'] == 'INSERTION'
    assert strain_parsed_vcf_dict['13-1950']['NC_002945.4'][714775]['ALT'] == 'A'


def test_summarise_vcf_outputs():
    pass_dict, insertion_dict, deletion_dict = \
        VSNPTreeMethods.summarise_gvcf_outputs(strain_parsed_vcf_dict=strain_parsed_vcf_dict)
    assert pass_dict['13-1941'] == 520
    assert insertion_dict['B13-0234'] == 148
    assert deletion_dict['13-1950'] == 44802


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
        VSNPTreeMethods.load_gvcf_snp_positions(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
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
    assert len(group_positions_set['af']['All']['NC_002945.4']) == 552
    assert len(group_positions_set['af']['Mbovis-01A3']['NC_002945.4']) == 552
    assert len(group_positions_set['suis1']['Bsuis1-09']['NC_017251.1']) == 70
    assert len(group_positions_set['suis1']['Bsuis1-09']['NC_017250.1']) == 48
    assert sorted(group_positions_set['suis1']['Bsuis1-09']['NC_017250.1'])[0] == 8810
    with pytest.raises(AssertionError):
        assert sorted(group_positions_set['suis1']['Bsuis1-09']['NC_017251.1'])[0] == 8810


def test_filter_snps():
    global filtered_group_positions
    filtered_group_positions = VSNPTreeMethods.filter_snps(group_positions_set)
    assert len(filtered_group_positions['af']['All']['NC_002945.4']) == 516
    assert len(filtered_group_positions['af']['Mbovis-01A3']['NC_002945.4']) == 516
    assert len(filtered_group_positions['suis1']['Bsuis1-09']['NC_017251.1']) == 70
    assert len(filtered_group_positions['suis1']['Bsuis1-09']['NC_017250.1']) == 48
    assert sorted(filtered_group_positions['suis1']['Bsuis1-09']['NC_017250.1'])[0] == 8810
    with pytest.raises(AssertionError):
        assert sorted(filtered_group_positions['suis1']['Bsuis1-09']['NC_017251.1'])[0] == 8810


def test_load_snp_sequence():
    global group_strain_snp_sequence, species_group_best_ref
    group_strain_snp_sequence, species_group_best_ref = \
        VSNPTreeMethods.load_snp_sequence(strain_parsed_vcf_dict=strain_parsed_vcf_dict,
                                          strain_consolidated_ref_dict=strain_consolidated_ref_dict,
                                          group_positions_set=filtered_group_positions,
                                          strain_groups=strain_groups,
                                          strain_species_dict=strain_species_dict,
                                          consolidated_ref_snp_positions=consolidated_ref_snp_positions)
    assert group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][1057] == 'G'
    assert group_strain_snp_sequence['af']['All']['13-1950']['NC_002945.4'][1057] == 'G'
    assert group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][714775] == 'G'
    assert group_strain_snp_sequence['af']['All']['13-1950']['NC_002945.4'][714775] == 'R'
    assert len(group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4']) == 516
    assert group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][1467394] == 'Y'
    with pytest.raises(KeyError):
        assert group_strain_snp_sequence['af']['All']['13-1941']['NC_017250.1']
    assert group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017250.1'][623503] == '-'
    with pytest.raises(KeyError):
        assert group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017251.1'][623503]
    assert species_group_best_ref['af']['All'] == 'NC_002945v4'
    assert species_group_best_ref['suis1']['All'] == 'NC_017251-NC_017250'


def test_remove_identical_calls():
    global non_identical_group_strain_snp_sequence, non_identical_group_positions
    non_identical_group_strain_snp_sequence, non_identical_group_positions = \
        VSNPTreeMethods.remove_identical_calls(group_strain_snp_sequence=group_strain_snp_sequence,
                                               consolidated_ref_snp_positions=consolidated_ref_snp_positions)
    assert non_identical_group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][103811] == 'G'
    assert non_identical_group_strain_snp_sequence['af']['All']['14-2093']['NC_002945.4'][103811] == 'C'
    assert len(non_identical_group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4']) == 66
    assert non_identical_group_strain_snp_sequence['af']['All']['13-1941']['NC_002945.4'][1467394] == 'Y'
    with pytest.raises(KeyError):
        assert non_identical_group_strain_snp_sequence['af']['All']['13-1941']['NC_017250.1']
    assert non_identical_group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017250.1'][623503] == '-'
    with pytest.raises(KeyError):
        assert non_identical_group_strain_snp_sequence['suis1']['All']['B13-0234']['NC_017251.1'][623503]
    assert len(non_identical_group_positions['af']['All']['NC_002945.4']) == 66
    assert len(non_identical_group_positions['af']['Mbovis-01A3']['NC_002945.4']) == 66
    assert len(non_identical_group_positions['suis1']['Bsuis1-09']['NC_017251.1']) == 24
    assert len(non_identical_group_positions['suis1']['Bsuis1-09']['NC_017250.1']) == 27
    assert sorted(non_identical_group_positions['suis1']['Bsuis1-09']['NC_017250.1'])[0] == 16405
    with pytest.raises(AssertionError):
        assert sorted(non_identical_group_positions['suis1']['Bsuis1-09']['NC_017251.1'])[0] == 16405


def test_create_multifasta():
    global group_folders, species_folders, group_fasta_dict
    group_folders, species_folders, group_fasta_dict = \
        VSNPTreeMethods.create_multifasta(group_strain_snp_sequence=non_identical_group_strain_snp_sequence,
                                          group_positions_set=non_identical_group_positions,
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
    assert species_group_order_dict['af']['All'] == ['13-1950', '14-2093', '13-1941', '03-1057', '13-3082',
                                                     'NC_002945v4']

    assert species_group_order_dict['suis1']['All'] == ['B13-0237', 'B13-0238', 'B13-0239', 'B13-0235', 'B13-0234',
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
        VSNPTreeMethods.annotate_snps(group_strain_snp_sequence=non_identical_group_strain_snp_sequence,
                                      full_best_ref_gbk_dict=full_best_ref_gbk_dict,
                                      strain_best_ref_set_dict=strain_best_ref_set_dict,
                                      ref_snp_positions=ref_snp_positions)
    assert species_group_annotated_snps_dict['af']['Mbovis-01']['NC_002945.4'][103811]['locus'] == 'BQ2027_MB0097C'
    assert species_group_annotated_snps_dict['suis1']['All']['NC_017250.1'][388552]['gene'] == 'dppD'
    assert species_group_annotated_snps_dict['suis1']['All']['NC_017251.1'][69252]['product'] == \
        'nucleoside/nucleotide kinase family protein'


def test_determine_snp_number():
    global species_group_snp_num_dict
    species_group_snp_num_dict = \
        VSNPTreeMethods.determine_snp_number(group_strain_snp_sequence=non_identical_group_strain_snp_sequence,
                                             species_group_best_ref=species_group_best_ref)
    assert species_group_snp_num_dict['af']['All']['NC_002945.4'][103811] == 4
    assert species_group_snp_num_dict['suis1']['All']['NC_017250.1'][16405] == 1
    with pytest.raises(KeyError):
        assert species_group_snp_num_dict['suis1']['All']['NC_017251.1'][8810]


def test_rank_snps():
    global species_group_snp_rank, species_group_num_snps
    species_group_snp_rank, species_group_num_snps \
        = VSNPTreeMethods.rank_snps(species_group_snp_num_dict=species_group_snp_num_dict)
    assert species_group_snp_rank['af']['All'][5]['NC_002945.4'][0] == 362818
    assert species_group_snp_rank['suis1']['All'][4]['NC_017250.1'][0] == 1195206
    with pytest.raises(KeyError):
        assert species_group_snp_rank['suis1']['All'][4]['NC_017251.1'][0]
    assert len(species_group_snp_rank['af']['All'][5]['NC_002945.4']) == 9
    assert len(species_group_snp_rank['suis1']['All'][1]['NC_017251.1']) == 24
    assert len(species_group_snp_rank['suis1']['All'][1]['NC_017251.1']) == 24
    assert len(species_group_snp_rank['suis1']['All'][4]['NC_017250.1']) == 1
    assert species_group_num_snps['af']['All'] == 66
    assert species_group_num_snps['suis1']['All'] == 51


def test_sort_snps():
    global species_group_sorted_snps
    species_group_sorted_snps = \
        VSNPTreeMethods.sort_snps(species_group_order_dict=species_group_order_dict,
                                  species_group_snp_rank=species_group_snp_rank,
                                  species_group_best_ref=species_group_best_ref,
                                  group_strain_snp_sequence=non_identical_group_strain_snp_sequence)
    assert species_group_sorted_snps['af']['All'][4]['NC_002945.4'][0] == 103811
    assert species_group_sorted_snps['af']['All'][1]['NC_002945.4'][-1] == 4279531
    assert species_group_sorted_snps['suis1']['All'][4]['NC_017250.1'][0] == 1195206
    assert species_group_sorted_snps['suis1']['All'][1]['NC_017250.1'][-1] == 1183517
    assert species_group_sorted_snps['suis1']['All'][1]['NC_017251.1'][-1] == 2065515
    assert len(species_group_sorted_snps['af']['All'][5]['NC_002945.4']) == 9
    assert len(species_group_sorted_snps['af']['All'][1]['NC_002945.4']) == 18
    assert len(species_group_sorted_snps['suis1']['All'][1]['NC_017251.1']) == 24
    assert len(species_group_sorted_snps['suis1']['All'][4]['NC_017250.1']) == 1


def test_create_summary_table():
    VSNPTreeMethods.create_summary_table(species_group_sorted_snps=species_group_sorted_snps,
                                         species_group_order_dict=species_group_order_dict,
                                         species_group_best_ref=species_group_best_ref,
                                         group_strain_snp_sequence=non_identical_group_strain_snp_sequence,
                                         species_group_annotated_snps_dict=species_group_annotated_snps_dict,
                                         species_group_num_snps=species_group_num_snps,
                                         summary_path=summary_path)
    assert len(glob(os.path.join(summary_path, '*.xlsx'))) == 9


def test_folder_prep():
    global deep_variant_path
    # Set the name, and create folders to hold VCF files for the test run of the pipeline
    deep_variant_path = os.path.join(file_path, 'deepvariant')
    make_path(deep_variant_path)
    # Copy the test files to the appropriate folders
    shutil.copyfile(src=os.path.join(file_path, 'B13-0234.gvcf.gz'),
                    dst=os.path.join(deep_variant_path, 'B13-0234.gvcf.gz'))
    shutil.copyfile(src=os.path.join(file_path, 'B13-0235.gvcf.gz'),
                    dst=os.path.join(deep_variant_path, 'B13-0235.gvcf.gz'))
    shutil.copyfile(src=os.path.join(file_path, 'B13-0237.gvcf.gz'),
                    dst=os.path.join(deep_variant_path, 'B13-0237.gvcf.gz'))
    shutil.copyfile(src=os.path.join(file_path, 'B13-0238.gvcf.gz'),
                    dst=os.path.join(deep_variant_path, 'B13-0238.gvcf.gz'))
    assert os.path.isfile(os.path.join(deep_variant_path, 'B13-0234.gvcf.gz'))


def test_vsnp_tree_run_deepvariant():
    vsnp_tree = VSNPTree(path=deep_variant_path,
                         threads=threads,
                         debug=False,
                         variant_caller='deepvariant',
                         filter_positions=False)
    vsnp_tree.main()
    assert os.path.isfile(os.path.join(deep_variant_path, 'summary_tables', 'suis1_All_sorted_table.xlsx'))


def test_remove_species_folders():
    for species_folder in species_folders:
        shutil.rmtree(species_folder)


def test_remove_alignments_folder():
    shutil.rmtree(fasta_path)


def test_remove_tree_folder():
    shutil.rmtree(tree_path)


def test_remove_report_folder():
    shutil.rmtree(summary_path)


def test_remove_run_folder():
    shutil.rmtree(deep_variant_path)


def test_remove_logs():
    logs = glob(os.path.join(file_path, '*.txt'))
    for log in logs:
        os.remove(log)
