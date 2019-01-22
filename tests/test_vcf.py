#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path
from vsnp.methods import Methods
from vsnp.vcf import VCF
import pytest
import os
__author__ = 'adamkoziol'

testpath = os.path.abspath(os.path.dirname(__file__))
filepath = os.path.join(testpath, 'files')
dependencypath = os.path.join(os.path.dirname(testpath), 'dependencies')

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
    assert len(file_list) == 6


def test_strain_dict():
    global strain_folder_dict
    strain_folder_dict = Methods.strain_list(fastq_files=file_list)
    for strain_folder, fastq_files in strain_folder_dict.items():
        if os.path.basename(strain_folder) in ['14-2093', '13-1941']:
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
        if strain_name in ['14-2093', '13-1941']:
            assert len(fastq_files) == 1
        else:
            assert len(fastq_files) == 2
        for symlink in fastq_files:
            assert os.path.islink(symlink)


def test_reformat_quality():
    global logfile, strain_qhist_dict, strain_lhist_dict
    logfile = os.path.join(filepath, 'log')
    strain_qhist_dict, strain_lhist_dict = Methods.run_reformat_reads(strain_fastq_dict=strain_fastq_dict,
                                                                      strain_name_dict=strain_name_dict,
                                                                      logfile=logfile)
    for strain_name, qhist_paths in strain_qhist_dict.items():
        for strain_qhist_file in qhist_paths:
            assert os.path.basename(strain_qhist_file).startswith(strain_name)
            assert strain_qhist_file.endswith('_qchist.csv')


def test_parse_reformat_quality():
    global strain_average_quality_dict, strain_qual_over_thirty_dict
    strain_average_quality_dict, strain_qual_over_thirty_dict = Methods.\
        parse_quality_histogram(strain_qhist_dict=strain_qhist_dict)
    assert strain_average_quality_dict['03-1057'] == [34.742922778562175, 30.805370650837197]


def test_parse_reformat_length():
    global strain_avg_read_lengths
    strain_avg_read_lengths = Methods.parse_length_histograms(strain_lhist_dict=strain_lhist_dict)
    assert strain_avg_read_lengths['03-1057'] == 225.23133


def test_file_size():
    global strain_fastq_size_dict
    strain_fastq_size_dict = Methods.find_fastq_size(strain_fastq_dict)
    assert strain_fastq_size_dict['13-1950'] == [32.82019233703613, 37.25274848937988]


def test_mash_sketch():
    global fastq_sketch_dict
    fastq_sketch_dict = Methods.call_mash_sketch(strain_fastq_dict=strain_fastq_dict,
                                                 strain_name_dict=strain_name_dict,
                                                 logfile=logfile)
    for strain, sketch_file in fastq_sketch_dict.items():
        assert os.path.isfile(sketch_file)


def test_mash_dist():
    global mash_dist_dict
    mash_dist_dict = Methods.call_mash_dist(strain_fastq_dict=strain_fastq_dict,
                                            strain_name_dict=strain_name_dict,
                                            fastq_sketch_dict=fastq_sketch_dict,
                                            ref_sketch_file=os.path.join(
                                                  dependencypath, 'mash', 'vsnp_reference.msh'),
                                            logfile=logfile)
    for strain, tab_output in mash_dist_dict.items():
        assert os.path.isfile(tab_output)


def test_mash_accession_species():
    global accession_species_dict
    accession_species_dict = Methods.parse_mash_accession_species(mash_species_file=os.path.join(
                                                  dependencypath, 'mash', 'species_accessions.csv'))
    assert accession_species_dict['NC_002945v4.fasta'] == 'af'


def test_mash_best_ref():
    global strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict
    strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict = \
        Methods.mash_best_ref(mash_dist_dict=mash_dist_dict,
                              accession_species_dict=accession_species_dict)
    assert strain_best_ref_dict['03-1057'] == 'NC_002945v4.fasta'
    assert strain_ref_matches_dict['03-1057'] == 968
    assert strain_species_dict['03-1057'] == 'af'


def test_reference_file_paths():
    global reference_link_path_dict
    reference_link_path_dict = Methods.reference_folder(strain_best_ref_dict=strain_best_ref_dict,
                                                        dependency_path=dependencypath)
    assert reference_link_path_dict['03-1057'] == 'mycobacterium/tbc/af2122/script_dependents/NC_002945v4.fasta'


def test_bowtie2_build():
    global strain_bowtie2_index_dict
    strain_bowtie2_index_dict = Methods.bowtie2_build(reference_link_path_dict=reference_link_path_dict,
                                                      dependency_path=dependencypath,
                                                      logfile=logfile)
    assert os.path.isfile(os.path.join(dependencypath, 'mycobacterium', 'tbc', 'af2122', 'script_dependents',
                                       'NC_002945v4.1.bt2'))
    assert os.path.split(strain_bowtie2_index_dict['03-1057'])[-1] == 'NC_002945v4'


def test_bowtie2_map():
    global strain_sorted_bam_dict
    strain_sorted_bam_dict = Methods.bowtie2_map(strain_fastq_dict=strain_fastq_dict,
                                                 strain_name_dict=strain_name_dict,
                                                 strain_bowtie2_index_dict=strain_bowtie2_index_dict,
                                                 threads=4,
                                                 logfile=logfile)
    for strain_name, sorted_bam in strain_sorted_bam_dict.items():
        assert os.path.isfile(sorted_bam)


def test_unmapped_reads_extract():
    global strain_unmapped_reads_dict
    strain_unmapped_reads_dict = Methods.extract_unmapped_reads(strain_sorted_bam_dict=strain_sorted_bam_dict,
                                                                strain_name_dict=strain_name_dict,
                                                                threads=4,
                                                                logfile=logfile)
    for strain_name, unmapped_reads_fastq in strain_unmapped_reads_dict.items():
        assert os.path.getsize(unmapped_reads_fastq) > 0


def test_skesa_assembled_unmapped():
    global strain_skesa_output_fasta_dict
    strain_skesa_output_fasta_dict = Methods.assemble_unmapped_reads(
        strain_unmapped_reads_dict=strain_unmapped_reads_dict,
        strain_name_dict=strain_name_dict,
        threads=4,
        logfile=logfile)
    assert os.path.getsize(strain_skesa_output_fasta_dict['03-1057']) == 0
    assert os.path.getsize(strain_skesa_output_fasta_dict['13-1941']) == 45462


def test_number_unmapped_contigs():
    global strain_unmapped_contigs_dict
    strain_unmapped_contigs_dict = Methods.assembly_stats(strain_skesa_output_fasta_dict=strain_skesa_output_fasta_dict)
    assert strain_unmapped_contigs_dict['13-1950'] == 0
    assert strain_unmapped_contigs_dict['13-1941'] == 37


def test_samtools_index():
    Methods.samtools_index(strain_sorted_bam_dict=strain_sorted_bam_dict,
                           strain_name_dict=strain_name_dict,
                           threads=4,
                           logfile=logfile)
    for strain_name, sorted_bam in strain_sorted_bam_dict.items():
        assert os.path.isfile(sorted_bam + '.bai')


def test_qualimap():
    global strain_qualimap_report_dict
    strain_qualimap_report_dict = Methods.run_qualimap(strain_sorted_bam_dict=strain_sorted_bam_dict,
                                                       strain_name_dict=strain_name_dict,
                                                       logfile=logfile)
    for strain_name, qualimap_report in strain_qualimap_report_dict.items():
        assert os.path.isfile(qualimap_report)


def test_qualimap_parse():
    global strain_qualimap_outputs_dict
    strain_qualimap_outputs_dict = Methods.parse_qualimap(strain_qualimap_report_dict=strain_qualimap_report_dict)
    assert strain_qualimap_outputs_dict['13-1950']['MappedReads'] == '374103(93.53%)'

