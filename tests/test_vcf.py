#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import filer, make_path
from cowsnphr_src.vcf_methods import VCFMethods
from datetime import datetime
from pathlib import Path
import multiprocessing
from glob import glob
import subprocess
import pytest
import shutil
import os

__author__ = 'adamkoziol'

test_path = os.path.abspath(os.path.dirname(__file__))
file_path = os.path.join(test_path, 'files', 'fastq')
dependency_path = os.path.join(os.path.dirname(test_path), 'dependencies')
summary_path = os.path.join(file_path, 'summary_tables')
vcf_path = os.path.join(file_path, 'vcf_files')
logfile = os.path.join(file_path, 'log')
threads = multiprocessing.cpu_count() - 1
home = str(Path.home())
# Define the start time
start_time = datetime.now()
deepvariant_version = '1.5.0'
file_list = list()
strain_folder_dict = dict()
strain_name_dict = dict()
make_path_folder = str()
strain_fastq_dict = dict()
strain_qhist_dict = dict()
strain_lhist_dict = dict()
strain_average_quality_dict = dict()
strain_qual_over_thirty_dict = dict()
strain_avg_read_lengths = dict()
strain_fastq_size_dict = dict()
fastq_sketch_dict = dict()
accession_species_dict = dict()
strain_best_ref_dict = dict()
strain_ref_matches_dict = dict()
strain_species_dict = dict()
mash_dist_dict = dict()
reference_link_path_dict = dict()
reference_link_dict = dict()
bowtie2_mapper_index_dict = dict()
bowtie2_reference_abs_path_dict = dict()
bowtie2_reference_dep_path_dict = dict()
bowtie2_reference_link_path_dict = dict()
strain_mapper_index_dict = dict()
strain_reference_abs_path_dict = dict()
strain_reference_dep_path_dict = dict()
bowtie2_sorted_bam_dict = dict()
strain_sorted_bam_dict = dict()
strain_unmapped_reads_dict = dict()
strain_skesa_output_fasta_dict = dict()
quast_report_dict = dict()
strain_examples_dict = dict()
strain_variant_path_dict = dict()
strain_call_variants_dict = dict()
strain_gvcf_tfrecords_dict = dict()
strain_vcf_dict = dict()
strain_ref_regions_dict = dict()
gvcf_num_high_quality_snps_dict = dict()
strain_num_high_quality_snps_dict = dict()


def test_deepvariant_version():
    global deepvariant_version
    cmd = 'docker run --rm google/deepvariant:{dvv} ' \
          '/opt/deepvariant/bin/call_variants --help'.format(dvv=deepvariant_version)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    assert out.decode().startswith('Code')


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
        filelist=[
            '03-1057_S10_L001_R1_001.fastq.gz',
            '03-1057_S10_L001_R2_001.fastq.gz',
            '13-1941_S4_L001_R1_001.fastq.gz',
            '13-1941_S4_L001_R2_001.fastq.gz'
        ],
        returndict=True
    )
    assert [file_name for file_name in filedict] == ['03-1057', '13-1941']


def test_normal_filer():
    fileset = filer(
        filelist=[
            '03-1057_S10_L001_R1_001.fastq.gz',
            '03-1057_S10_L001_R2_001.fastq.gz',
            '13-1941_S4_L001_R1_001.fastq.gz',
            '13-1941_S4_L001_R2_001.fastq.gz'
        ]
    )
    assert fileset == {'03-1057', '13-1941'}


def test_missing_file_filer():
    fileset = filer(
        filelist=[
            '03-1057_S10_L001_R1_001.fastq.gz',
            '03-1057_S10_L001_R2_001.fastq.gz',
            '13-1941_S4_L001_R1_001.fastq.gz'
        ]
    )
    assert fileset == {'03-1057', '13-1941'}


def test_non_illumina_filer():
    fileset = filer(
        filelist=[
            '03-1057_R1.fastq.gz',
            '03-1057_R2.fastq.gz',
            '13-1941_R1.fastq.gz',
            '13-1941_R2.fastq.gz'
        ]
    )
    assert fileset == {'03-1057', '13-1941'}


def test_multiple_differences_filer():
    fileset = filer(
        filelist=[
            '03-1057_1.fastq',
            '03-1057_2.fastq',
            '13-1941_1.fastq.gz',
            '13-1941_2.fastq'
        ]
    )
    assert fileset == {'03-1057', '13-1941'}


def test_no_directions_filer():
    fileset = filer(
        filelist=[
            '03-1057.fastq.gz',
            '13-1941_S4_L001.fastq'
        ]
    )
    assert fileset == {'03-1057', '13-1941'}


def test_vcf_file_list_no_files():
    with pytest.raises(AssertionError):
        VCFMethods.file_list(path=test_path)


def test_vcf_file_list():
    global file_list
    file_list = VCFMethods.file_list(path=file_path)
    assert len(file_list) == 7


def test_strain_dict():
    global strain_folder_dict
    strain_folder_dict = VCFMethods.strain_list(fastq_files=file_list)
    for strain_folder, fastq_files in strain_folder_dict.items():
        if '13-1941' in strain_folder:
            assert len(fastq_files) == 1
        else:
            assert len(fastq_files) == 2


def test_strain_namer_no_input():
    strain_names = VCFMethods.strain_namer(strain_folders=str())
    assert len(strain_names) == 0


def test_strain_namer_working():
    global strain_name_dict
    strain_name_dict = VCFMethods.strain_namer(strain_folders=strain_folder_dict)
    assert [strain for strain in strain_name_dict] == ['13-1941', '13-1950', 'B13-0234', 'NC_002695']


def test_make_path():
    global make_path_folder
    make_path_folder = os.path.join(test_path, 'test_folder')
    make_path(make_path_folder)
    assert os.path.isdir(make_path_folder)


def test_rm_path():
    os.rmdir(make_path_folder)
    assert os.path.isdir(make_path_folder) is False


def test_strain_linker():
    global strain_fastq_dict
    strain_fastq_dict = VCFMethods.file_link(
        strain_folder_dict=strain_folder_dict,
        strain_name_dict=strain_name_dict
    )
    assert [strain for strain in strain_fastq_dict] == ['13-1941', '13-1950', 'B13-0234', 'NC_002695']
    for strain_name, fastq_files in strain_fastq_dict.items():
        if strain_name == '13-1941':
            assert len(fastq_files) == 1
        else:
            assert len(fastq_files) == 2
        for symlink in fastq_files:
            assert os.path.islink(symlink)


def test_reformat_quality():
    global strain_qhist_dict, strain_lhist_dict
    strain_qhist_dict, strain_lhist_dict = VCFMethods.run_reformat_reads(
        strain_fastq_dict=strain_fastq_dict,
        strain_name_dict=strain_name_dict,
        logfile=logfile
    )
    for strain_name, qhist_paths in strain_qhist_dict.items():
        for strain_qhist_file in qhist_paths:
            assert os.path.basename(strain_qhist_file).startswith(strain_name)
            assert strain_qhist_file.endswith('_qchist.csv')


def test_parse_reformat_quality():
    global strain_average_quality_dict, strain_qual_over_thirty_dict
    strain_average_quality_dict, strain_qual_over_thirty_dict = VCFMethods. \
        parse_quality_histogram(strain_qhist_dict=strain_qhist_dict)
    assert strain_average_quality_dict['13-1950'] == [33.82559209616877, 28.64100810052621]
    assert strain_qual_over_thirty_dict['13-1950'] == [84.5724480421427, 61.085547466494404]


def test_parse_reformat_length():
    global strain_avg_read_lengths
    strain_avg_read_lengths = VCFMethods.parse_length_histograms(strain_lhist_dict=strain_lhist_dict)
    assert strain_avg_read_lengths['13-1950'] == 230.9919625


def test_file_size():
    global strain_fastq_size_dict
    strain_fastq_size_dict = VCFMethods.find_fastq_size(strain_fastq_dict)
    assert strain_fastq_size_dict['13-1950'] == [32.82019233703613, 37.25274848937988]


def test_mash_sketch():
    global fastq_sketch_dict
    fastq_sketch_dict = VCFMethods.call_mash_sketch(
        strain_fastq_dict=strain_fastq_dict,
        strain_name_dict=strain_name_dict,
        logfile=logfile
    )
    for strain, sketch_file in fastq_sketch_dict.items():
        assert os.path.isfile(sketch_file)


def test_mash_dist():
    global mash_dist_dict
    mash_dist_dict = VCFMethods.call_mash_dist(
        strain_fastq_dict=strain_fastq_dict,
        strain_name_dict=strain_name_dict,
        fastq_sketch_dict=fastq_sketch_dict,
        ref_sketch_file=os.path.join(
            dependency_path, 'mash', 'reference.msh'),
        logfile=logfile
    )
    for strain, tab_output in mash_dist_dict.items():
        assert os.path.isfile(tab_output)


def test_mash_accession_species():
    global accession_species_dict
    accession_species_dict = VCFMethods.parse_mash_accession_species(mash_species_file=os.path.join(
        dependency_path, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_mash_best_ref():
    global strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict
    strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict = \
        VCFMethods.mash_best_ref(
            mash_dist_dict=mash_dist_dict,
            accession_species_dict=accession_species_dict,
            min_matches=500
        )
    assert strain_best_ref_dict['13-1950'] == 'NC_002945v4.fasta'
    assert strain_ref_matches_dict['13-1950'] == 916
    assert strain_species_dict['13-1950'] == 'af'


def test_reference_file_paths():
    global reference_link_path_dict, reference_link_dict
    reference_link_path_dict, reference_link_dict = VCFMethods.reference_folder(
        strain_best_ref_dict=strain_best_ref_dict,
        dependency_path=dependency_path
    )
    assert reference_link_path_dict['13-1950'] == 'mycobacterium/tbc/af2122/script_dependents/NC_002945v4.fasta'


def test_bowtie2_build():
    global bowtie2_mapper_index_dict, bowtie2_reference_abs_path_dict, bowtie2_reference_dep_path_dict, \
        bowtie2_reference_link_path_dict
    bowtie2_reference_link_path_dict = dict()
    for strain_name, ref_link in reference_link_path_dict.items():
        bowtie2_reference_link_path_dict[strain_name] = ref_link
    bowtie2_mapper_index_dict, bowtie2_reference_abs_path_dict, bowtie2_reference_dep_path_dict = \
        VCFMethods.index_ref_genome(
            reference_link_path_dict=bowtie2_reference_link_path_dict,
            dependency_path=dependency_path,
            logfile=logfile,
            reference_mapper='bowtie2'
        )
    assert os.path.isfile(os.path.join(
        dependency_path, 'brucella', 'suis1', 'script_dependents', 'NC_017251-NC_017250.1.bt2')
    )
    assert os.path.split(bowtie2_mapper_index_dict['B13-0234'])[-1] == 'NC_017251-NC_017250'


def test_faidx_reference():
    VCFMethods.faidx_ref_genome(
        reference_link_path_dict=bowtie2_reference_link_path_dict,
        dependency_path=dependency_path,
        logfile=logfile
    )
    for strain_name, ref_link in reference_link_path_dict.items():
        # Set the absolute path, and strip off the file extension for use in the build call
        ref_fasta = os.path.abspath(os.path.join(dependency_path, ref_link))
        assert os.path.isfile(ref_fasta + '.fai')


def test_consolidate_index_dicts():
    global strain_mapper_index_dict, strain_reference_abs_path_dict, strain_reference_dep_path_dict
    for strain_name, index_file in bowtie2_mapper_index_dict.items():
        strain_mapper_index_dict[strain_name] = index_file
        strain_reference_abs_path_dict[strain_name] = bowtie2_reference_abs_path_dict[strain_name]
        strain_reference_dep_path_dict[strain_name] = bowtie2_reference_dep_path_dict[strain_name]
    assert os.path.split(strain_mapper_index_dict['13-1950'])[-1] == 'NC_002945v4'
    assert os.path.split(strain_mapper_index_dict['B13-0234'])[-1] == 'NC_017251-NC_017250'


def test_bowtie2_map():
    global bowtie2_sorted_bam_dict
    # Create a dictionary of the subset of strains to process with bowtie2
    bowtie2_strain_fastq_dict = dict()
    for strain_name, fastq_files in strain_fastq_dict.items():
        if strain_name in ['B13-0234', 'NC_002695', '13-1950', '13-1941']:
            bowtie2_strain_fastq_dict[strain_name] = fastq_files
    bowtie2_sorted_bam_dict = VCFMethods.map_ref_genome(
        strain_fastq_dict=bowtie2_strain_fastq_dict,
        strain_name_dict=strain_name_dict,
        strain_mapper_index_dict=strain_mapper_index_dict,
        threads=threads,
        logfile=logfile,
        reference_mapper='bowtie2'
    )
    for strain_name, sorted_bam in bowtie2_sorted_bam_dict.items():
        assert os.path.isfile(sorted_bam)


def test_merge_bowtie_bam_dict():
    global strain_sorted_bam_dict
    strain_sorted_bam_dict = dict()
    for strain_name, sorted_bam in bowtie2_sorted_bam_dict.items():
        strain_sorted_bam_dict[strain_name] = sorted_bam
    assert strain_sorted_bam_dict['13-1950']
    assert strain_sorted_bam_dict['B13-0234']


def test_unmapped_reads_extract():
    global strain_unmapped_reads_dict
    strain_unmapped_reads_dict = VCFMethods.extract_unmapped_reads(
        strain_sorted_bam_dict=strain_sorted_bam_dict,
        strain_name_dict=strain_name_dict,
        threads=threads,
        logfile=logfile
    )
    for strain_name, unmapped_reads_fastq in strain_unmapped_reads_dict.items():
        assert os.path.getsize(unmapped_reads_fastq) > 0


def test_skesa_assembled_unmapped():
    global strain_skesa_output_fasta_dict
    strain_skesa_output_fasta_dict = VCFMethods.assemble_unmapped_reads(
        strain_unmapped_reads_dict=strain_unmapped_reads_dict,
        strain_name_dict=strain_name_dict,
        threads=threads,
        logfile=logfile
    )
    assert os.path.getsize(strain_skesa_output_fasta_dict['13-1941']) > 100


def test_quast():
    global quast_report_dict
    quast_report_dict = VCFMethods.quast(
        strain_skesa_output_fasta_dict=strain_skesa_output_fasta_dict,
        strain_unmapped_reads_dict=strain_unmapped_reads_dict,
        strain_sorted_bam_dict=strain_sorted_bam_dict,
        threads=threads,
        logfile=logfile
    )

    assembly_file = strain_skesa_output_fasta_dict['13-1941']
    output_dir = os.path.dirname(assembly_file)
    quast_report = os.path.join(output_dir, 'report.tsv')
    assert os.path.isfile(quast_report)


def test_parse_quast_report():
    VCFMethods.parse_quast_report(
        quast_report_dict=quast_report_dict,
        summary_path=summary_path
    )
    assert os.path.isfile(os.path.join(summary_path, 'assembly_report.tsv'))
    with open(os.path.join(summary_path, 'assembly_report.tsv'), 'r') as quast:
        header = quast.readline().split('\t')
        data = quast.readline().split('\t')
    assert header[0] == 'Assembly'
    assert data[1] == '37'


def test_samtools_index():
    VCFMethods.samtools_index(
        strain_sorted_bam_dict=strain_sorted_bam_dict,
        strain_name_dict=strain_name_dict,
        threads=threads,
        logfile=logfile
    )
    for strain_name, sorted_bam in strain_sorted_bam_dict.items():
        assert os.path.isfile(sorted_bam + '.bai')


# def test_deepvariant_make_examples():
#     global strain_examples_dict, strain_variant_path_dict, strain_gvcf_tfrecords_dict, vcf_path
#     reduced_strain_sorted_bam_dict = dict()
#     reduced_strain_sorted_bam_dict['13-1941'] = strain_sorted_bam_dict['13-1941']
#     strain_examples_dict, strain_variant_path_dict, strain_gvcf_tfrecords_dict = \
#         VCFMethods.deepvariant_make_examples(
#             strain_sorted_bam_dict=reduced_strain_sorted_bam_dict,
#             strain_name_dict=strain_name_dict,
#             strain_reference_abs_path_dict=strain_reference_abs_path_dict,
#             vcf_path=vcf_path,
#             home=home,
#             threads=threads,
#             deepvariant_version=deepvariant_version,
#             logfile=logfile
#         )
#     assert len(strain_examples_dict['13-1941']) == threads
#     for strain_name, gvcf_tfrecord in strain_gvcf_tfrecords_dict.items():
#         gvcf_tfrecord = gvcf_tfrecord.split('@')[0]
#         if strain_name != 'NC_002695':
#             assert len(glob('{gvcf_tfrecord}*.gz'.format(gvcf_tfrecord=gvcf_tfrecord))) == threads
#         else:
#             assert len(glob('{gvcf_tfrecord}*.gz'.format(gvcf_tfrecord=gvcf_tfrecord))) == 0
#
#
# def test_deepvariant_call_variants():
#     global strain_call_variants_dict
#     strain_call_variants_dict = \
#         VCFMethods.deepvariant_call_variants(
#             strain_variant_path_dict=strain_variant_path_dict,
#             strain_name_dict=strain_name_dict,
#             home=home,
#             vcf_path=vcf_path,
#             threads=threads,
#             logfile=logfile,
#             variant_caller='deepvariant',
#             deepvariant_version=deepvariant_version
#         )
#     assert os.path.getsize(strain_call_variants_dict['13-1941']) > 100
#     with pytest.raises(KeyError):
#         assert strain_call_variants_dict['NC_002695']
#
#
# def test_deepvariant_postprocess_variants():
#     global strain_vcf_dict
#     strain_vcf_dict = \
#         VCFMethods.deepvariant_postprocess_variants_multiprocessing(
#             strain_call_variants_dict=strain_call_variants_dict,
#             strain_variant_path_dict=strain_variant_path_dict,
#             strain_name_dict=strain_name_dict,
#             strain_reference_abs_path_dict=strain_reference_abs_path_dict,
#             strain_gvcf_tfrecords_dict=strain_gvcf_tfrecords_dict,
#             vcf_path=vcf_path,
#             home=home,
#             logfile=logfile,
#             deepvariant_version=deepvariant_version,
#             threads=threads
#         )
#     assert os.path.getsize(strain_vcf_dict['13-1941']) > 100
#     with pytest.raises(KeyError):
#         assert strain_vcf_dict['NC_002695']


def test_combined_deepvariant_call():
    global strain_vcf_dict
    strain_vcf_dict = VCFMethods.deepvariant_run_container(
        strain_sorted_bam_dict=strain_sorted_bam_dict,
        strain_reference_abs_path_dict=strain_reference_abs_path_dict,
        strain_name_dict=strain_name_dict,
        gpu=False,
        container_platform='docker',
        home=home,
        working_path=None,
        version=deepvariant_version,
        platform='WGS',
        threads=threads,
        logfile=logfile
    )
    assert os.path.getsize(strain_vcf_dict['13-1941']) > 100
    with pytest.raises(KeyError):
        assert strain_vcf_dict['NC_002695']


def test_copy_test_vcf_files():
    """
    Copy VCF files from test folder to supplement the lone deepvariant-created VCF file. Populate the strain_vcf_dict
    dictionary with these VCF files
    """
    # Set the absolute path of the test folder containing the VCF files
    vcf_test_path = os.path.join(test_path, 'files', 'vcf')
    # Create a list of all the VCF files
    vcf_files = glob(os.path.join(vcf_test_path, '*.gvcf.gz'))
    for strain_name, strain_folder in strain_name_dict.items():
        if strain_name not in strain_vcf_dict:
            # Set the name of the output .vcf file
            vcf_base_name = '{sn}.gvcf.gz'.format(sn=strain_name)
            out_vcf = os.path.join(strain_folder, vcf_base_name)
            # Don't try to copy the file if the original exists
            for vcf_file in vcf_files:
                if os.path.basename(vcf_file) == vcf_base_name:
                    shutil.copyfile(vcf_file, out_vcf)
            # Update the dictionary
            strain_vcf_dict[strain_name] = out_vcf
            if strain_name not in ['NC_002695', '13-1950']:
                assert os.path.isfile(out_vcf)


def test_parse_gvcf():
    global gvcf_num_high_quality_snps_dict
    gvcf_strain_vcf_dict = dict()
    for strain_name, vcf_file in strain_vcf_dict.items():
        if strain_name != '13-1950':
            gvcf_strain_vcf_dict[strain_name] = vcf_file
    gvcf_num_high_quality_snps_dict = VCFMethods.parse_gvcf(strain_vcf_dict=gvcf_strain_vcf_dict)
    assert gvcf_num_high_quality_snps_dict['13-1941'] == 398
    assert gvcf_num_high_quality_snps_dict['B13-0234'] == 63
    assert gvcf_num_high_quality_snps_dict['NC_002695'] == 0


def test_copy_vcf_files():
    VCFMethods.copy_vcf_files(strain_vcf_dict=strain_vcf_dict,
                              vcf_path=vcf_path)
    assert os.path.isdir(vcf_path)
    assert len(glob(os.path.join(vcf_path, '*.gvcf.gz'))) == 3


def test_consolidate_high_quality_snp_dicts():
    global strain_num_high_quality_snps_dict
    strain_num_high_quality_snps_dict = dict()
    for strain_name, num_high_quality_snps in gvcf_num_high_quality_snps_dict.items():
        strain_num_high_quality_snps_dict[strain_name] = num_high_quality_snps
    assert strain_num_high_quality_snps_dict['13-1941'] == 398
    assert strain_num_high_quality_snps_dict['B13-0234'] == 63
    assert strain_num_high_quality_snps_dict['NC_002695'] == 0
    with pytest.raises(KeyError):
        assert strain_num_high_quality_snps_dict['13-1950']


def test_remove_bt2_indexes():
    for strain_name, ref_link in reference_link_path_dict.items():
        # Set the absolute path, and strip off the file extension for use in the build call
        ref_abs_path = os.path.dirname(os.path.abspath(os.path.join(dependency_path, ref_link)))
        bt2_files = glob(os.path.join(ref_abs_path, '*.bt2'))
        for bt2_index in bt2_files:
            os.remove(bt2_index)
        bt2_files = glob(os.path.join(ref_abs_path, '*.bt2'))
        assert not bt2_files


def test_remove_fai_files():
    for strain_name, ref_link in reference_link_path_dict.items():
        # Set the absolute path, and strip off the file extension for use in the build call
        ref_abs_path = os.path.dirname(os.path.abspath(os.path.join(dependency_path, ref_link)))
        fai_files = glob(os.path.join(ref_abs_path, '*.fai'))
        for fai in fai_files:
            os.remove(fai)
        fai_files = glob(os.path.join(ref_abs_path, '*.fai'))
        assert not fai_files


def test_remove_logs():
    logs = glob(os.path.join(file_path, '*.txt'))
    for log in logs:
        os.remove(log)


def test_remove_vcf_path():
    shutil.rmtree(vcf_path)


def test_remove_report_folder():
    shutil.rmtree(summary_path)


def test_remove_working_dir():
    for strain_name, working_dir in strain_name_dict.items():
        shutil.rmtree(working_dir)
