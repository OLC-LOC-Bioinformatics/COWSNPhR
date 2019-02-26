#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import SetupLogging
from vsnp.vsnp_vcf_methods import VCFMethods
from datetime import datetime
from pathlib import Path
import logging
import os

__author__ = 'adamkoziol'


class VCF(object):
    def main(self):
        """
        Run all the VCF-specific methods
        """
        self.fastq_manipulation()
        self.best_reference_calculation()
        self.reference_mapping()
        self.stat_calculation()
        self.snp_calling()
        self.typing()
        self.report()

    def fastq_manipulation(self):
        """
        Determine the number of strains to process. Create strain-specific working directories with relative symlinks
        to FASTQ files
        """
        logging.info('Locating FASTQ files, creating strain-specific working directories and symlinks to files')
        fastq_files = VCFMethods.file_list(path=self.path)
        logging.debug('FASTQ files: \n{fastq_files}'.format(fastq_files='\n'.join(fastq_files)))
        strain_folder_dict = VCFMethods.strain_list(fastq_files=fastq_files)
        self.strain_name_dict = VCFMethods.strain_namer(strain_folders=strain_folder_dict)
        logging.debug('Strain names: \n{strain_names}'.format(strain_names='\n'.join(sorted(self.strain_name_dict))))
        self.strain_fastq_dict = VCFMethods.file_link(strain_folder_dict=strain_folder_dict,
                                                      strain_name_dict=self.strain_name_dict)
        logging.debug(
            'Strain-specific symlinked FASTQ files: \n{symlinks}'.format(
                symlinks='\n'.join(['{strain_name}: {fastq_files}'.format(strain_name=sn, fastq_files=ff)
                                    for sn, ff in self.strain_fastq_dict.items()])))

    def best_reference_calculation(self):
        """
        Determine the best reference using MASH
        """
        logging.info('Running MASH analyses')
        fastq_sketch_dict = VCFMethods.call_mash_sketch(strain_fastq_dict=self.strain_fastq_dict,
                                                        strain_name_dict=self.strain_name_dict,
                                                        logfile=self.logfile)
        logging.info('Parsing MASH outputs to determine closest reference genomes')
        mash_dist_dict = VCFMethods.call_mash_dist(strain_fastq_dict=self.strain_fastq_dict,
                                                   strain_name_dict=self.strain_name_dict,
                                                   fastq_sketch_dict=fastq_sketch_dict,
                                                   ref_sketch_file=os.path.join(
                                                       self.dependency_path, 'mash', 'vsnp_reference.msh'),
                                                   logfile=self.logfile)
        logging.info('Loading reference genome: species dictionary')

        accession_species_dict = VCFMethods.parse_mash_accession_species(mash_species_file=os.path.join(
            self.dependency_path, 'mash', 'species_accessions.csv'))

        logging.info('Determining closest reference genome and extracting corresponding species from MASH outputs')
        self.strain_best_ref_dict, self.strain_ref_matches_dict, self.strain_species_dict = \
            VCFMethods.mash_best_ref(mash_dist_dict=mash_dist_dict,
                                     accession_species_dict=accession_species_dict)

    def reference_mapping(self):
        """
        Perform reference mapping with bowtie2, and attempt to assemble any unmapped reads into contigs with SKESA
        """
        logging.info('Extracting paths to reference genomes')
        reference_link_path_dict, reference_link_dict \
            = VCFMethods.reference_folder(strain_best_ref_dict=self.strain_best_ref_dict,
                                          dependency_path=self.dependency_path)
        logging.info('Running bowtie2 build')
        strain_bowtie2_index_dict, self.strain_reference_abs_path_dict, self.strain_reference_dep_path_dict = \
            VCFMethods.bowtie2_build(reference_link_path_dict=reference_link_path_dict,
                                     dependency_path=self.dependency_path,
                                     logfile=self.logfile)
        logging.info('Running bowtie2 reference mapping')
        self.strain_sorted_bam_dict = VCFMethods.bowtie2_map(
            strain_fastq_dict=self.strain_fastq_dict,
            strain_name_dict=self.strain_name_dict,
            strain_bowtie2_index_dict=strain_bowtie2_index_dict,
            threads=self.threads,
            logfile=self.logfile)
        logging.info('Extracting unmapped reads')
        strain_unmapped_reads_dict = VCFMethods.extract_unmapped_reads(
            strain_sorted_bam_dict=self.strain_sorted_bam_dict,
            strain_name_dict=self.strain_name_dict,
            threads=self.threads,
            logfile=self.logfile)
        logging.info('Attempting to assemble unmapped reads with SKESA')
        self.strain_skesa_output_fasta_dict = VCFMethods.assemble_unmapped_reads(
            strain_unmapped_reads_dict=strain_unmapped_reads_dict,
            strain_name_dict=self.strain_name_dict,
            threads=self.threads,
            logfile=self.logfile)

    def stat_calculation(self):
        """
        Calculate raw stats on FASTQ file size, quality and length distributions of FASTQ reads, and qualimap-generated
        statistics on reference mapping
        """
        logging.info('Calculating quality and length distributions of FASTQ reads')
        self.strain_qhist_dict, \
            self.strain_lhist_dict = VCFMethods.run_reformat_reads(strain_fastq_dict=self.strain_fastq_dict,
                                                                   strain_name_dict=self.strain_name_dict,
                                                                   logfile=self.logfile)
        self.strain_average_quality_dict, self.strain_qual_over_thirty_dict = \
            VCFMethods.parse_quality_histogram(strain_qhist_dict=self.strain_qhist_dict)
        self.strain_avg_read_lengths = VCFMethods.parse_length_histograms(strain_lhist_dict=self.strain_lhist_dict)
        logging.info('Calculating size of FASTQ files')
        self.strain_fastq_size_dict = VCFMethods.find_fastq_size(self.strain_fastq_dict)
        logging.info('Counting of contigs in assemblies of unmapped reads')
        self.strain_unmapped_contigs_dict = VCFMethods.assembly_stats(
            strain_skesa_output_fasta_dict=self.strain_skesa_output_fasta_dict)
        VCFMethods.samtools_index(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                  strain_name_dict=self.strain_name_dict,
                                  threads=self.threads,
                                  logfile=self.logfile)
        logging.info('Running qualimap analyses on sorted BAM files')
        strain_qualimap_report_dict = VCFMethods.run_qualimap(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                                              strain_name_dict=self.strain_name_dict,
                                                              logfile=self.logfile)
        self.strain_qualimap_outputs_dict = VCFMethods.parse_qualimap(
            strain_qualimap_report_dict=strain_qualimap_report_dict)

    def snp_calling(self):
        """
        Prep files for SNP calling. Use deepvariant to call SNPs. Parse the outputs from deepvariant
        """
        logging.info('Preparing files for SNP calling with deepvariant make_examples')
        strain_examples_dict, strain_variant_path_dict, strain_gvcf_tfrecords_dict = \
            VCFMethods.deepvariant_make_examples(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                                 strain_name_dict=self.strain_name_dict,
                                                 strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
                                                 vcf_path=os.path.join(self.path, 'vcf_files'),
                                                 home=self.home,
                                                 threads=self.threads)
        logging.info('Calling variants with deepvariant call_variants')
        strain_call_variants_dict = \
            VCFMethods.deepvariant_call_variants_multiprocessing(strain_variant_path_dict=strain_variant_path_dict,
                                                                 strain_name_dict=self.strain_name_dict,
                                                                 dependency_path=self.dependency_path,
                                                                 vcf_path=os.path.join(self.path, 'vcf_files'),
                                                                 home=self.home,
                                                                 threads=self.threads,
                                                                 logfile=self.logfile)
        logging.info('Creating VCF files with deepvariant postprocess_variants')
        strain_vcf_dict = \
            VCFMethods.deepvariant_postprocess_variants(
                strain_call_variants_dict=strain_call_variants_dict,
                strain_variant_path_dict=strain_variant_path_dict,
                strain_name_dict=self.strain_name_dict,
                strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
                strain_gvcf_tfrecords_dict=strain_gvcf_tfrecords_dict,
                vcf_path=os.path.join(self.path, 'vcf_files'),
                home=self.home,
                logfile=self.logfile)
        # logging.info('Quality filtering VCF files with vcftools')
        # strain_vcf_dict = VCFMethods.filter_vcf(strain_unfiltered_vcf_dict=strain_unfiltered_vcf_dict,
        #                                         strain_name_dict=self.strain_name_dict,
        #                                         logfile=self.logfile)
        logging.info('Parsing VCF outputs')
        self.strain_num_high_quality_snps_dict = VCFMethods.parse_variants(strain_vcf_dict=strain_vcf_dict)
        VCFMethods.copy_vcf_files(strain_vcf_dict=strain_vcf_dict,
                                  vcf_path=os.path.join(self.path, 'vcf_files'))

    def typing(self):
        """
        Perform typing analyses including spoligotyping, and the subsequence extraction of binary, octal, hexadecimal,
        and sb codes. Also determine MLST profiles of samples
        """
        logging.info('Searching for matches to spoligo sequences')
        strain_spoligo_stats_dict = VCFMethods.bait_spoligo(strain_fastq_dict=self.strain_fastq_dict,
                                                            strain_name_dict=self.strain_name_dict,
                                                            spoligo_file=os.path.join(self.dependency_path,
                                                                                      'mycobacterium',
                                                                                      'spacers.fasta'),
                                                            threads=self.threads,
                                                            logfile=self.logfile,
                                                            kmer=25)
        logging.info('Calculating binary, octal, and hexadecimal codes')
        self.strain_binary_code_dict, \
            self.strain_octal_code_dict, \
            self.strain_hexadecimal_code_dict = \
            VCFMethods.parse_spoligo(strain_spoligo_stats_dict=strain_spoligo_stats_dict)
        logging.info('Extracting sb codes')
        self.strain_sbcode_dict = VCFMethods.extract_sbcode(
            strain_reference_dep_path_dict=self.strain_reference_dep_path_dict,
            strain_octal_code_dict=self.strain_octal_code_dict)
        logging.info('Performing MLST analyses')
        VCFMethods.brucella_mlst(seqpath=self.path,
                                 mlst_db_path=os.path.join(self.dependency_path, 'brucella', 'MLST'),
                                 logfile=self.logfile)
        logging.info('Parsing MLST outputs')
        self.strain_mlst_dict = VCFMethods.parse_mlst_report(strain_name_dict=self.strain_name_dict,
                                                             mlst_report=os.path.join(self.path, 'reports', 'mlst.csv'))

    def report(self):
        """
        Create the .xlsx report consistent with the legacy vSNP format
        """
        logging.info('Creating report')
        VCFMethods.create_vcf_report(
            start_time=self.start_time,
            strain_species_dict=self.strain_species_dict,
            strain_best_ref_dict=self.strain_best_ref_dict,
            strain_fastq_size_dict=self.strain_fastq_size_dict,
            strain_average_quality_dict=self.strain_average_quality_dict,
            strain_qual_over_thirty_dict=self.strain_qual_over_thirty_dict,
            strain_qualimap_outputs_dict=self.strain_qualimap_outputs_dict,
            strain_avg_read_lengths=self.strain_avg_read_lengths,
            strain_unmapped_contigs_dict=self.strain_unmapped_contigs_dict,
            strain_num_high_quality_snps_dict=self.strain_num_high_quality_snps_dict,
            strain_mlst_dict=self.strain_mlst_dict,
            strain_octal_code_dict=self.strain_octal_code_dict,
            strain_sbcode_dict=self.strain_sbcode_dict,
            strain_hexadecimal_code_dict=self.strain_hexadecimal_code_dict,
            strain_binary_code_dict=self.strain_binary_code_dict,
            report_path=self.report_path)

    def __init__(self, path, threads, debug=False):
        """
        :param path: type STR: Path of folder containing FASTQ files
        :param threads: type INT: Number of threads to use in the analyses
        :param debug: type BOOL: Boolean of whether debug level logs are printed to terminal
        """
        SetupLogging(debug=debug)
        # Determine the path in which the sequence files are located. Allow for ~ expansion
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        # Ensure that the path exists
        assert os.path.isdir(self.path), 'Invalid path specified: {path}'.format(path=self.path)
        logging.debug('Supplied sequence path: \n{path}'.format(path=self.path))
        # Initialise class variables
        self.threads = threads
        self.report_path = os.path.join(self.path, 'reports')
        # Extract the path of the folder containing this script
        self.script_path = os.path.abspath(os.path.dirname(__file__))
        # Use the script path to set the absolute path of the dependencies folder
        self.dependency_path = os.path.join(os.path.dirname(self.script_path), 'dependencies')
        assert os.path.isdir(self.dependency_path), 'Something went wrong with the install. Cannot locate the ' \
                                                   'dependencies folder in: {sp}'.format(sp=self.script_path)
        self.logfile = os.path.join(self.path, 'log')
        self.start_time = datetime.now()
        self.home = str(Path.home())
        self.strain_name_dict = dict()
        self.strain_fastq_dict = dict()
        self.strain_best_ref_dict = dict()
        self.strain_ref_matches_dict = dict()
        self.strain_species_dict = dict()
        self.strain_sorted_bam_dict = dict()
        self.strain_reference_abs_path_dict = dict()
        self.strain_reference_dep_path_dict = dict()
        self.strain_skesa_output_fasta_dict = dict()
        self.strain_qhist_dict = dict()
        self.strain_lhist_dict = dict()
        self.strain_average_quality_dict = dict()
        self.strain_qual_over_thirty_dict = dict()
        self.strain_avg_read_lengths = dict()
        self.strain_fastq_size_dict = dict()
        self.strain_unmapped_contigs_dict = dict()
        self.strain_qualimap_outputs_dict = dict()
        self.strain_num_high_quality_snps_dict = dict()
        self.strain_binary_code_dict = dict()
        self.strain_octal_code_dict = dict()
        self.strain_hexadecimal_code_dict = dict()
        self.strain_sbcode_dict = dict()
        self.strain_mlst_dict = dict()
