#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from vsnp.vcf_methods import Methods
from datetime import datetime
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
        fastq_files = Methods.file_list(path=self.path)
        logging.debug('FASTQ files: \n{fastq_files}'.format(fastq_files='\n'.join(fastq_files)))
        strain_folder_dict = Methods.strain_list(fastq_files=fastq_files)
        self.strain_name_dict = Methods.strain_namer(strain_folders=strain_folder_dict)
        logging.debug('Strain names: \n{strain_names}'.format(strain_names='\n'.join(sorted(self.strain_name_dict))))
        self.strain_fastq_dict = Methods.file_link(strain_folder_dict=strain_folder_dict,
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
        fastq_sketch_dict = Methods.call_mash_sketch(strain_fastq_dict=self.strain_fastq_dict,
                                                     strain_name_dict=self.strain_name_dict,
                                                     logfile=self.logfile)
        logging.info('Parsing MASH outputs to determine closest reference genomes')
        mash_dist_dict = Methods.call_mash_dist(strain_fastq_dict=self.strain_fastq_dict,
                                                strain_name_dict=self.strain_name_dict,
                                                fastq_sketch_dict=fastq_sketch_dict,
                                                ref_sketch_file=os.path.join(
                                                    self.dependencypath, 'mash', 'vsnp_reference.msh'),
                                                logfile=self.logfile)
        logging.info('Loading reference genome: species dictionary')

        accession_species_dict = Methods.parse_mash_accession_species(mash_species_file=os.path.join(
            self.dependencypath, 'mash', 'species_accessions.csv'))

        logging.info('Determining closest reference genome and extracting corresponding species from MASH outputs')
        self.strain_best_ref_dict, self.strain_ref_matches_dict, self.strain_species_dict = \
            Methods.mash_best_ref(mash_dist_dict=mash_dist_dict,
                                  accession_species_dict=accession_species_dict)

    def reference_mapping(self):
        """
        Perform reference mapping with bowtie2, and attempt to assemble any unmapped reads into contigs with SKESA
        """
        logging.info('Extracting paths to reference genomes')
        reference_link_path_dict, reference_link_dict \
            = Methods.reference_folder(strain_best_ref_dict=self.strain_best_ref_dict,
                                       dependency_path=self.dependencypath)
        logging.info('Running bowtie2 build')
        strain_bowtie2_index_dict, self.strain_reference_abs_path_dict, self.strain_reference_dep_path_dict = \
            Methods.bowtie2_build(reference_link_path_dict=reference_link_path_dict,
                                  dependency_path=self.dependencypath,
                                  logfile=self.logfile)
        logging.info('Running bowtie2 reference mapping')
        self.strain_sorted_bam_dict = Methods.bowtie2_map(
            strain_fastq_dict=self.strain_fastq_dict,
            strain_name_dict=self.strain_name_dict,
            strain_bowtie2_index_dict=strain_bowtie2_index_dict,
            threads=self.threads,
            logfile=self.logfile)
        logging.info('Extracting unmapped reads')
        strain_unmapped_reads_dict = Methods.extract_unmapped_reads(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                                                    strain_name_dict=self.strain_name_dict,
                                                                    threads=self.threads,
                                                                    logfile=self.logfile)
        logging.info('Attempting to assemble unmapped reads with SKESA')
        self.strain_skesa_output_fasta_dict = Methods.assemble_unmapped_reads(
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
            self.strain_lhist_dict = Methods.run_reformat_reads(strain_fastq_dict=self.strain_fastq_dict,
                                                                strain_name_dict=self.strain_name_dict,
                                                                logfile=self.logfile)
        self.strain_average_quality_dict, self.strain_qual_over_thirty_dict = \
            Methods.parse_quality_histogram(strain_qhist_dict=self.strain_qhist_dict)
        self.strain_avg_read_lengths = Methods.parse_length_histograms(strain_lhist_dict=self.strain_lhist_dict)
        logging.info('Calculating size of FASTQ files')
        self.strain_fastq_size_dict = Methods.find_fastq_size(self.strain_fastq_dict)
        logging.info('Counting of contigs in assemblies of unmapped reads')
        self.strain_unmapped_contigs_dict = Methods.assembly_stats(
            strain_skesa_output_fasta_dict=self.strain_skesa_output_fasta_dict)
        Methods.samtools_index(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                               strain_name_dict=self.strain_name_dict,
                               threads=self.threads,
                               logfile=self.logfile)
        logging.info('Running qualimap analyses on sorted BAM files')
        strain_qualimap_report_dict = Methods.run_qualimap(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                                           strain_name_dict=self.strain_name_dict,
                                                           logfile=self.logfile)
        self.strain_qualimap_outputs_dict = Methods.parse_qualimap(
            strain_qualimap_report_dict=strain_qualimap_report_dict)

    def snp_calling(self):
        """
        Prep files for SNP calling. Use FreeBayes to call SNPs. Parse the outputs from FreeBayes
        """
        logging.info('Creating regions file for reference genomes')
        strain_ref_regions_dict = Methods.reference_regions(
            strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
            logfile=self.logfile)
        logging.info('Running FreeBayes on samples')
        strain_vcf_dict = Methods.freebayes(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                                            strain_name_dict=self.strain_name_dict,
                                            strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
                                            strain_ref_regions_dict=strain_ref_regions_dict,
                                            threads=self.threads,
                                            logfile=self.logfile)
        logging.info('Parsing VCF outputs')
        self.strain_num_high_quality_snps_dict, strain_filtered_vcf_dict = \
            Methods.parse_vcf(strain_vcf_dict=strain_vcf_dict)
        Methods.copy_vcf_files(strain_filtered_vcf_dict=strain_filtered_vcf_dict,
                               vcf_path=os.path.join(self.path, 'vcf_files'))

    def typing(self):
        """
        Perform typing analyses including spoligotyping, and the subsequence extraction of binary, octal, hexadecimal,
        and sb codes. Also determine MLST profiles of samples
        """
        logging.info('Searching for matches to spoligo sequences')
        strain_spoligo_stats_dict = Methods.bait_spoligo(strain_fastq_dict=self.strain_fastq_dict,
                                                         strain_name_dict=self.strain_name_dict,
                                                         spoligo_file=os.path.join(self.dependencypath,
                                                                                   'mycobacterium',
                                                                                   'spacers.fasta'),
                                                         threads=self.threads,
                                                         logfile=self.logfile,
                                                         kmer=25)
        logging.info('Calculating binary, octal, and hexadecimal codes')
        self.strain_binary_code_dict, \
            self.strain_octal_code_dict, \
            self.strain_hexadecimal_code_dict = \
            Methods.parse_spoligo(strain_spoligo_stats_dict=strain_spoligo_stats_dict)
        logging.info('Extracting sb codes')
        self.strain_sbcode_dict = Methods.extract_sbcode(
            strain_reference_dep_path_dict=self.strain_reference_dep_path_dict,
            strain_octal_code_dict=self.strain_octal_code_dict)
        logging.info('Performing MLST analyses')
        Methods.brucella_mlst(seqpath=self.path,
                              mlst_db_path=os.path.join(self.dependencypath, 'brucella', 'MLST'),
                              logfile=self.logfile)
        logging.info('Parsing MLST outputs')
        self.strain_mlst_dict = Methods.parse_mlst_report(strain_name_dict=self.strain_name_dict,
                                                          mlst_report=os.path.join(self.path, 'reports', 'mlst.csv'))

    def report(self):
        """
        Create the .xlsx report consistent with the legacy vSNP format
        """
        logging.info('Creating report')
        Methods.create_vcf_report(
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
        # make_path(self.report_path)
        # assert os.path.isdir(self.report_path), 'Could not create report path as requested: {rp}' \
        #     .format(rp=self.report_path)
        # Extract the path of the folder containing this script
        self.scriptpath = os.path.abspath(os.path.dirname(__file__))
        # Use the script path to set the absolute path of the dependencies folder
        self.dependencypath = os.path.join(os.path.dirname(self.scriptpath), 'dependencies')
        assert os.path.isdir(self.dependencypath), 'Something went wrong with the install. Cannot locate the ' \
                                                   'dependencies folder in: {sp}'.format(sp=self.scriptpath)
        self.logfile = os.path.join(self.path, 'log')
        self.start_time = datetime.now()
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
