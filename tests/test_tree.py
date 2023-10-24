#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import filer, make_path
from cowsnphr_src.tree_methods import TreeMethods
from cowsnphr_src.cowsnphr import COWSNPhR
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
tree_path = os.path.join(file_path, 'tree_files')
deep_variant_path = os.path.join(file_path, 'deepvariant')
matrix_path = os.path.join(file_path, 'snp_matrix')
logfile = os.path.join(file_path, 'log')
threads = multiprocessing.cpu_count() - 1
# Define the start time
start_time = datetime.now()
iupac = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
    '-': ['-']
}
# Initialise variables
vcf_object = object
file_list = list()
strain_vcf_dict = dict()
accession_species_dict = dict()
strain_parsed_vcf_dict = dict()
strain_best_ref_dict = dict()
strain_best_ref_set_dict = dict()
strain_species_dict = dict()
strain_best_ref_fasta_dict = dict()
reference_link_path_dict = dict()
reference_link_dict = dict()
reference_strain_dict = dict()
strain_consolidated_ref_dict = dict()
defining_snp_dict = dict()
consolidated_ref_snp_positions = dict()
strain_snp_positions = dict()
ref_snp_positions = dict()
group_positions_set = dict()
strain_groups = dict()
filtered_group_positions = dict()
coords_dict = dict()
mask_pos_dict = dict()
supplied_mask_pos_dict = dict()
filter_reasons = dict()
group_strain_snp_sequence = dict()
species_group_best_ref = dict()
group_folders = list()
species_folders = list()
group_fasta_dict = dict()
ident_group_positions = dict()
full_best_ref_gbk_dict = dict()
species_group_trees = dict()
species_group_order_dict = dict()
species_group_annotated_snps_dict = dict()
species_group_snp_num_dict = dict()
species_group_snp_rank = dict()
species_group_num_snps = dict()
species_group_sorted_snps = dict()
translated_snp_residue_dict = dict()
ref_translated_snp_residue_dict = dict()


def test_invalid_path():
    with pytest.raises(AssertionError):
        COWSNPhR(
            seq_path='not a real path',
            ref_path=dependency_path,
            threads=1,
            working_path=file_path,
            mask_file=None,
            gpu=None,
            platform='WGS',
            container_platform='docker',
            debug=False
        )


def test_no_threads():
    with pytest.raises(TypeError):
        assert COWSNPhR(
            seq_path='not a real path',
            ref_path=dependency_path,
            working_path=file_path,
            mask_file=None,
            platform='WGS',
            container_platform='docker',
            gpu=None,
            debug=False
        )


def test_valid_path():
    global vcf_object
    vcf_object = COWSNPhR(
        seq_path=file_path,
        ref_path=dependency_path,
        threads=1,
        working_path=file_path,
        mask_file=None,
        gpu=None,
        platform='WGS',
        container_platform='docker',
        debug=False
    )
    assert vcf_object


def test_invalid_tilde_path():
    COWSNPhR(
        seq_path='~',
        ref_path=dependency_path,
        threads=1,
        working_path=file_path,
        mask_file=None,
        gpu=None,
        platform='WGS',
        container_platform='docker',
        debug=False
    )


def test_empty_filer():
    fileset = filer(filelist=list())
    assert fileset == set()


def test_empty_filer_dict():
    filedict = filer(
        filelist=list(),
        returndict=True
    )
    assert filedict == dict()


def test_normal_filer_dict():
    filedict = filer(
        filelist=['03-1057.vcf', '13-1941.vcf'],
        extension='vcf',
        returndict=True
    )
    assert [file_name for file_name in filedict] == ['03-1057', '13-1941']


def test_vcf_file_list_no_files():
    with pytest.raises(AssertionError):
        TreeMethods.file_list(file_path=test_path)


def test_vcf_file_list():
    global file_list
    file_list = TreeMethods.file_list(file_path=file_path)
    assert len(file_list) == 5


def test_strain_dict():
    global strain_vcf_dict
    strain_vcf_dict = TreeMethods.strain_list(vcf_files=file_list)
    assert [strain for strain in sorted(strain_vcf_dict)] == \
           ['B13-0234', 'B13-0235', 'B13-0237', 'B13-0238', 'B13-0239']


def test_accession_species():
    global accession_species_dict
    accession_species_dict = TreeMethods.parse_accession_species(ref_species_file=os.path.join(
        dependency_path, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_gvcf_load():
    global strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict
    gvcf_vcf_dict = dict()
    for strain_name, vcf_file in strain_vcf_dict.items():
        if strain_name != '13-1950':
            gvcf_vcf_dict[strain_name] = vcf_file
    strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict = \
        TreeMethods.load_vcf(strain_vcf_dict=gvcf_vcf_dict)
    for strain_name, best_ref in strain_best_ref_dict.items():
        assert best_ref in ['NC_017250.1', 'NC_017251.1']
    assert strain_best_ref_set_dict['B13-0234'] == {'NC_017250.1', 'NC_017251.1'}
    with pytest.raises(KeyError):
        assert strain_parsed_vcf_dict['13-1941']
    assert strain_parsed_vcf_dict['B13-0235']['NC_017250.1'][8810]['QUAL'] == '70.1'


def test_summarise_vcf_outputs():
    pass_dict, insertion_dict, deletion_dict = \
        TreeMethods.summarise_gvcf_outputs(strain_parsed_vcf_dict=strain_parsed_vcf_dict)
    assert insertion_dict['B13-0234'] == 55
    assert pass_dict['B13-0235'] == 373
    assert deletion_dict['B13-0237'] == 34364
    with pytest.raises(KeyError):
        assert deletion_dict['13-1950']


def test_determine_ref_species():
    global strain_species_dict, strain_best_ref_fasta_dict
    strain_species_dict, strain_best_ref_fasta_dict = \
        TreeMethods.determine_ref_species(
            strain_best_ref_dict=strain_best_ref_dict,
            accession_species_dict=accession_species_dict
        )
    assert strain_species_dict['B13-0234'] == 'suis1'
    assert strain_best_ref_fasta_dict['B13-0234'] == 'NC_017251-NC_017250.fasta'


def test_reference_path():
    global reference_link_path_dict, reference_link_dict, reference_strain_dict
    reference_link_path_dict, reference_link_dict, reference_strain_dict = \
        TreeMethods.reference_folder(
            strain_best_ref_fasta_dict=strain_best_ref_fasta_dict,
            dependency_path=dependency_path
        )
    assert 'brucella/suis1/script_dependents' in reference_link_path_dict['B13-0234']
    assert 'brucella/suis1/script_dependents/NC_017251-NC_017250.fasta' in reference_strain_dict['B13-0234']
    assert 'brucella/suis1/script_dependents/NC_017251-NC_017250.fasta' in \
           reference_link_dict['NC_017251-NC_017250.fasta']
    with pytest.raises(KeyError):
        assert reference_link_path_dict['13-1950']


def test_consolidate_group_ref_genomes():
    global strain_consolidated_ref_dict
    strain_consolidated_ref_dict = \
        TreeMethods.consolidate_group_ref_genomes(
            reference_link_dict=reference_link_dict,
            strain_best_ref_dict=strain_best_ref_dict
        )
    assert strain_consolidated_ref_dict['B13-0234'] == 'NC_017251-NC_017250'
    assert strain_consolidated_ref_dict['B13-0238'] == 'NC_017251-NC_017250'
    with pytest.raises(KeyError):
        assert strain_consolidated_ref_dict['13-1950']


def test_extract_defining_snps():
    global defining_snp_dict
    defining_snp_dict = TreeMethods.extract_defining_snps(
        reference_strain_dict=reference_strain_dict,
        strain_species_dict=strain_species_dict
    )
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
        TreeMethods.load_gvcf_snp_positions(
            strain_parsed_vcf_dict=strain_parsed_vcf_dict,
            strain_consolidated_ref_dict=strain_consolidated_ref_dict
        )
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017250.1'][1197913] == 'T'
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017250.1'][16405] == 'C'
    assert consolidated_ref_snp_positions['NC_017251-NC_017250']['NC_017251.1'][50509] == 'A'
    assert sorted(strain_snp_positions['B13-0234']['NC_017251.1'])[0] == 1068
    assert sorted(strain_snp_positions['B13-0234']['NC_017250.1'])[0] == 8810
    with pytest.raises(KeyError):
        assert strain_snp_positions['B13-0234']['NC_002945.4']
        assert strain_snp_positions['13-1950']
    assert ref_snp_positions['NC_017250.1'][1197913] == 'T'
    with pytest.raises(KeyError):
        assert ref_snp_positions['NC_017251.1'][16405]
    assert ref_snp_positions['NC_017250.1'][16405] == 'C'
    assert ref_snp_positions['NC_017251.1'][50509] == 'A'


def test_determine_groups():
    global group_positions_set, strain_groups, strain_species_dict
    group_positions_set, strain_groups, strain_species_dict = \
        TreeMethods.group_strains(strain_snp_positions=strain_snp_positions)
    assert strain_groups['B13-0234'] == ['group']
    assert 694277 in group_positions_set['species']['group']['NC_017250.1']
    assert strain_species_dict['B13-0237'] == 'species'
    with pytest.raises(KeyError):
        assert strain_groups['13-1950']


def test_density_filter_snps():
    global filtered_group_positions
    filtered_group_positions = TreeMethods.density_filter_snps(group_positions_set)
    assert 2816 in filtered_group_positions['species']['group']['NC_017250.1']
    with pytest.raises(AssertionError):
        assert 364560 in filtered_group_positions['species']['group']['NC_017250.1']


def test_mask_ref_genome():
    global coords_dict
    coords_dict = TreeMethods.mask_ref_genome(reference_strain_dict, logfile)
    assert 'NC_017251-NC_017250_mummer_maskfile.coords' in [os.path.basename(value) for value in coords_dict.values()]


def test_filter_coordinates():
    global mask_pos_dict
    mask_pos_dict = TreeMethods.determine_coordinates(strain_groups, coords_dict)
    assert 294908 in mask_pos_dict['B13-0234']['group']['NC_017250.1']
    with pytest.raises(KeyError):
        assert 43351 in mask_pos_dict['13-1941']['Mbovis-01']['NC_017250.1']


def test_load_supplied_mask():
    global supplied_mask_pos_dict
    supplied_mask_pos_dict = \
        TreeMethods.load_supplied_mask(
            strain_groups=strain_groups,
            maskfile=os.path.join(dependency_path, 'brucella/suis1/script_dependents/maskfile.bed')
        )
    print('\n')
    print(supplied_mask_pos_dict)
    assert 43351 in supplied_mask_pos_dict['species']['group']['NC_017250.1']
    assert 43355 not in supplied_mask_pos_dict['group']['group']['NC_017250.1']


def test_filter_masked_snps():
    global filtered_group_positions, filter_reasons
    filtered_masked_group_positions, filter_reasons = \
        TreeMethods.filter_masked_snp_positions(
            group_positions_set=group_positions_set,
            filtered_group_positions=filtered_group_positions,
            mask_pos_dict=mask_pos_dict,
            supplied_mask_pos_dict=supplied_mask_pos_dict
        )
    assert filter_reasons['species']['group']['NC_017250.1'][5306] == ['density']
    assert len(filtered_group_positions['species']['group']['NC_017251.1']) == 609
    assert len(filtered_group_positions['species']['group']['NC_017250.1']) == 366
    assert sorted(filtered_group_positions['species']['group']['NC_017250.1'])[0] == 2816
    with pytest.raises(KeyError):
        assert sorted(filtered_group_positions['suis1']['Bsuis1-09']['NC_017251.1'])[0] == 8810


def test_load_snp_sequence():
    global group_strain_snp_sequence, species_group_best_ref
    group_strain_snp_sequence, species_group_best_ref = \
        TreeMethods.load_snp_sequence(
            strain_parsed_vcf_dict=strain_parsed_vcf_dict,
            strain_consolidated_ref_dict=strain_consolidated_ref_dict,
            group_positions_set=filtered_group_positions,
            strain_groups=strain_groups,
            strain_species_dict=strain_species_dict,
            consolidated_ref_snp_positions=consolidated_ref_snp_positions,
            iupac=iupac
        )
    assert group_strain_snp_sequence['species']['group']['B13-0234']['NC_017250.1'][2816] == 'T'
    assert len(group_strain_snp_sequence['species']['group']['B13-0234']['NC_017250.1']) == 256
    with pytest.raises(KeyError):
        assert group_strain_snp_sequence['af']['All']['13-1941']['NC_017250.1']
    assert species_group_best_ref['species']['group'] == 'NC_017251-NC_017250'


def test_remove_identical_calls():
    global ident_group_positions
    ident_group_positions = \
        TreeMethods.find_identical_calls(group_strain_snp_sequence=group_strain_snp_sequence)
    assert 12689 in ident_group_positions['species']['group']['NC_017250.1']


def test_create_multifasta():
    global group_folders, species_folders, group_fasta_dict
    base_dict = dict()
    for key, value in reference_link_dict.items():
        base_dict[os.path.splitext(key)[0]] = value
    group_folders, species_folders, group_fasta_dict = \
        TreeMethods.create_multifasta(
            group_strain_snp_sequence=group_strain_snp_sequence,
            fasta_path=fasta_path,
            group_positions_set=group_positions_set,
            strain_parsed_vcf_dict=strain_parsed_vcf_dict,
            species_group_best_ref=species_group_best_ref,
            reference_strain_dict=base_dict,
            ident_group_positions=ident_group_positions
        )
    assert len(group_folders) == 1
    assert len(species_folders) == 1
    for species, group_dict in group_fasta_dict.items():
        for group, fasta in group_dict.items():
            assert os.path.getsize(fasta) > 100


def test_run_fasttree():
    global species_group_trees
    species_group_trees = TreeMethods.run_fasttree(
        group_fasta_dict=group_fasta_dict,
        strain_consolidated_ref_dict=strain_consolidated_ref_dict,
        strain_groups=strain_groups,
        logfile=logfile
    )
    for species, group_dict in species_group_trees.items():
        for group, options_dict in group_dict.items():
            assert os.path.getsize(options_dict['best_tree']) > 100


def test_parse_tree_order():
    global species_group_order_dict
    species_group_order_dict = TreeMethods.parse_tree_order(species_group_trees=species_group_trees)
    assert species_group_order_dict['species']['group'] == \
           ['B13-0234', 'B13-0238', 'NC_017251-NC_017250', 'B13-0239', 'B13-0235', 'B13-0237']


def test_copy_trees():
    TreeMethods.copy_trees(
        species_group_trees=species_group_trees,
        tree_path=tree_path
    )
    assert len(glob(os.path.join(tree_path, '*_All*'))) == 0
    assert len(glob(os.path.join(tree_path, '*best_tree*'))) == 1
    assert len(glob(os.path.join(tree_path, '*'))) == 1


def test_prokka():
    TreeMethods.prokka(reference_strain_dict=reference_strain_dict,
                       logfile=logfile)
    for strain_name, ref_fasta in reference_strain_dict.items():
        gbk_file = ref_fasta.replace('.fasta', '.gbk')
        assert os.path.isfile(gbk_file)


def test_load_genbank_file():
    global full_best_ref_gbk_dict
    full_best_ref_gbk_dict = TreeMethods.load_genbank_file_single(reference_strain_dict=reference_strain_dict)
    assert full_best_ref_gbk_dict['NC_017251.1'][28315].type == 'tRNA'
    assert full_best_ref_gbk_dict['NC_017250.1'][185334].qualifiers['gene'] == ['ileS']


def test_annotate_snps():
    global species_group_annotated_snps_dict
    species_group_annotated_snps_dict = \
        TreeMethods.annotate_snps(
            group_strain_snp_sequence=group_strain_snp_sequence,
            full_best_ref_gbk_dict=full_best_ref_gbk_dict,
            strain_best_ref_set_dict=strain_best_ref_set_dict,
            ref_snp_positions=ref_snp_positions
        )
    assert species_group_annotated_snps_dict['species']['group']['NC_017251.1'][1631902]['locus'] == 'IBCBLKLE_01571'
    assert species_group_annotated_snps_dict['species']['group']['NC_017251.1'][1451059]['gene'] == 'yfcG_2'
    assert species_group_annotated_snps_dict['species']['group']['NC_017251.1'][187595]['product'] == \
           'Heme chaperone HemW'


def test_determine_snp_number():
    global species_group_snp_num_dict
    species_group_snp_num_dict = \
        TreeMethods.determine_snp_number(
            group_strain_snp_sequence=group_strain_snp_sequence,
            species_group_best_ref=species_group_best_ref
        )
    assert species_group_snp_num_dict['species']['group']['NC_017251.1'][1263210] == 1
    assert species_group_snp_num_dict['species']['group']['NC_017251.1'][29269] == 5
    assert species_group_snp_num_dict['species']['group']['NC_017251.1'][1066226] == 4
    with pytest.raises(KeyError):
        assert species_group_snp_num_dict['species']['group']['NC_017251.1'][8810]


def test_determine_aa_sequence():
    global translated_snp_residue_dict, ref_translated_snp_residue_dict
    # Add the best reference genome to the reference strain dictionary
    best_ref = species_group_best_ref['species']['group']
    reference_strain_dict[best_ref] = reference_strain_dict['B13-0234']
    translated_snp_residue_dict, ref_translated_snp_residue_dict = \
        TreeMethods.determine_aa_sequence(
            group_strain_snp_sequence=group_strain_snp_sequence,
            species_group_best_ref=species_group_best_ref,
            strain_parsed_vcf_dict=strain_parsed_vcf_dict,
            species_group_annotated_snps_dict=species_group_annotated_snps_dict,
            reference_strain_dict=reference_strain_dict,
            species_group_snp_num_dict=species_group_snp_num_dict,
            iupac=iupac
        )
    assert translated_snp_residue_dict['species']['group']['B13-0234']['NC_017251.1'][1263210]['snp_nt_seq_raw'] == 'T'
    assert translated_snp_residue_dict['species']['group']['B13-0234']['NC_017251.1'][1263210]['snp_nt_seq_alt'] == 'A'
    assert translated_snp_residue_dict['species']['group']['B13-0234']['NC_017251.1'][1263210]['ref_aa_seq_cds'] == 'L'

    assert ref_translated_snp_residue_dict['species']['group']['NC_017251.1'][2055257]['ref_nt_seq_raw'] == 'A'
    assert ref_translated_snp_residue_dict['species']['group']['NC_017251.1'][2055257]['ref_nt_seq_cds'] == 'G'
    assert ref_translated_snp_residue_dict['species']['group']['NC_017251.1'][2055257]['cds_strand'] == -1


def test_create_snp_matrix():
    TreeMethods.create_snp_matrix(
        species_group_best_ref=species_group_best_ref,
        group_strain_snp_sequence=group_strain_snp_sequence,
        matrix_path=matrix_path
    )
    assert os.path.isfile(os.path.join(matrix_path, 'snv_matrix.tsv'))


def test_rank_snps():
    global species_group_snp_rank, species_group_num_snps
    species_group_snp_rank, species_group_num_snps \
        = TreeMethods.rank_snps(species_group_snp_num_dict=species_group_snp_num_dict)
    assert species_group_snp_rank['species']['group'][5]['NC_017250.1'][0] == 8810
    assert species_group_snp_rank['species']['group'][4]['NC_017250.1'][0] == 432146
    with pytest.raises(KeyError):
        assert species_group_snp_rank['suis1']['All'][4]['NC_017251.1'][0]
    assert len(species_group_snp_rank['species']['group'][5]['NC_017250.1']) == 12
    assert len(species_group_snp_rank['species']['group'][1]['NC_017251.1']) == 21
    assert species_group_num_snps['species']['group'] == 102


def test_sort_snps():
    global species_group_sorted_snps
    species_group_sorted_snps = \
        TreeMethods.sort_snps(
            species_group_order_dict=species_group_order_dict,
            species_group_snp_rank=species_group_snp_rank,
            species_group_best_ref=species_group_best_ref,
            group_strain_snp_sequence=group_strain_snp_sequence
        )
    assert species_group_sorted_snps['species']['group'][4]['NC_017250.1'][0] == 432146
    assert species_group_sorted_snps['species']['group'][1]['NC_017250.1'][-1] == 1183517
    assert species_group_sorted_snps['species']['group'][1]['NC_017251.1'][-1] == 2065515
    assert len(species_group_sorted_snps['species']['group'][1]['NC_017251.1']) == 21
    assert len(species_group_sorted_snps['species']['group'][4]['NC_017250.1']) == 3


def test_create_nt_summary_table():
    global summary_path
    make_path(summary_path)
    TreeMethods.create_summary_table(
        species_group_sorted_snps=species_group_sorted_snps,
        species_group_order_dict=species_group_order_dict,
        species_group_best_ref=species_group_best_ref,
        group_strain_snp_sequence=group_strain_snp_sequence,
        species_group_annotated_snps_dict=species_group_annotated_snps_dict,
        translated_snp_residue_dict=translated_snp_residue_dict,
        ref_translated_snp_residue_dict=ref_translated_snp_residue_dict,
        species_group_num_snps=species_group_num_snps,
        summary_path=summary_path,
        molecule='nt'
    )
    assert os.path.isfile(os.path.join(summary_path, 'nt_snv_sorted_table.xlsx')) == 1


def test_create_aa_summary_table():
    TreeMethods.create_summary_table(
        species_group_sorted_snps=species_group_sorted_snps,
        species_group_order_dict=species_group_order_dict,
        species_group_best_ref=species_group_best_ref,
        group_strain_snp_sequence=group_strain_snp_sequence,
        species_group_annotated_snps_dict=species_group_annotated_snps_dict,
        translated_snp_residue_dict=translated_snp_residue_dict,
        ref_translated_snp_residue_dict=ref_translated_snp_residue_dict,
        species_group_num_snps=species_group_num_snps,
        summary_path=summary_path,
        molecule='aa'
    )
    assert os.path.isfile(os.path.join(summary_path, 'aa_snv_sorted_table.xlsx')) == 1


def test_folder_prep():
    global deep_variant_path
    # Set the name, and create folders to hold VCF files for the test run of the pipeline
    make_path(deep_variant_path)
    # Copy the test files to the appropriate folders
    shutil.copyfile(
        src=os.path.join(file_path, 'B13-0234.gvcf.gz'),
        dst=os.path.join(deep_variant_path, 'B13-0234.gvcf.gz')
    )
    shutil.copyfile(
        src=os.path.join(file_path, 'B13-0235.gvcf.gz'),
        dst=os.path.join(deep_variant_path, 'B13-0235.gvcf.gz')
    )
    shutil.copyfile(
        src=os.path.join(file_path, 'B13-0237.gvcf.gz'),
        dst=os.path.join(deep_variant_path, 'B13-0237.gvcf.gz')
    )
    shutil.copyfile(
        src=os.path.join(file_path, 'B13-0238.gvcf.gz'),
        dst=os.path.join(deep_variant_path, 'B13-0238.gvcf.gz')
    )
    assert os.path.isfile(os.path.join(deep_variant_path, 'B13-0234.gvcf.gz'))


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


def test_remove_matrix_folder():
    shutil.rmtree(matrix_path)


def test_remove_logs():
    logs = glob(os.path.join(file_path, '*.txt'))
    for log in logs:
        os.remove(log)
