#!/usr/bin/env python3

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import multiprocessing
from pathlib import Path
import os

# Local imports
from cowsnphr_src.version import __version__
from olctools.accessoryFunctions.accessoryFunctions import SetupLogging
from cowsnphr_src.vcf_methods import VCFMethods
from cowsnphr_src.tree_methods import TreeMethods


__author__ = 'adamkoziol'


class COWSNPhR(object):

    def main(self):
        self.fastq_manipulation()
        self.reference_mapping()
        self.snp_calling()
        self.load_snps()
        self.phylogenetic_trees()
        self.annotate_snps()
        self.order_snps()
        self.create_report()

    def fastq_manipulation(self):
        """
        Determine the number of strains to process. Create strain-specific working directories with relative symlinks
        to FASTQ files
        """
        logging.info('Locating FASTQ files, creating strain-specific working directories and symlinks to files')
        fastq_files = VCFMethods.file_list(path=self.seq_path)
        logging.info('FASTQ files: \n{fastq_files}'.format(fastq_files='\n'.join(fastq_files)))
        strain_folder_dict = VCFMethods.strain_list(fastq_files=fastq_files)
        self.strain_name_dict = VCFMethods.strain_namer(strain_folders=strain_folder_dict)
        if self.debug:
            logging.info('Strain names: \n{strain_names}'.format(strain_names='\n'.join(sorted(self.strain_name_dict))))
        self.strain_fastq_dict = VCFMethods.file_link(
            strain_folder_dict=strain_folder_dict,
            strain_name_dict=self.strain_name_dict
        )
        if self.debug:
            logging.info(
                'Strain-specific symlinked FASTQ files: \n{symlinks}'.format(
                    symlinks='\n'.join(['{strain_name}: {fastq_files}'.format(strain_name=sn, fastq_files=ff)
                                        for sn, ff in self.strain_fastq_dict.items()])))

    def reference_mapping(self):
        """
        Perform reference mapping with bowtie2
        """
        logging.info('Extracting paths to reference genomes')
        self.ref_file()
        logging.info('Running bowtie2 build')
        strain_bowtie2_index_dict, self.strain_reference_abs_path_dict, self.strain_reference_dep_path_dict = \
            VCFMethods.index_ref_genome(
                reference_link_path_dict=self.reference_strain_dict,
                dependency_path=self.ref_path,
                logfile=self.logfile,
                reference_mapper='bowtie2'
            )
        logging.info('Creating .fai index file of {ref}'.format(ref=self.ref_strain))
        VCFMethods.faidx_ref_genome(
            reference_link_path_dict=self.reference_strain_dict,
            dependency_path=self.ref_path,
            logfile=self.logfile
        )
        logging.info('Running bowtie2 reference mapping')
        self.strain_sorted_bam_dict = VCFMethods.map_ref_genome(
            strain_fastq_dict=self.strain_fastq_dict,
            strain_name_dict=self.strain_name_dict,
            strain_mapper_index_dict=strain_bowtie2_index_dict,
            threads=self.threads,
            logfile=self.logfile,
            reference_mapper='bowtie2'
        )
        logging.debug('Sorted BAM files: \n{files}'.format(
            files='\n'.join(['{strain_name}: {bam_file}'.format(strain_name=sn, bam_file=bf)
                             for sn, bf in self.strain_sorted_bam_dict.items()])))
        logging.info('Indexing sorted BAM files')
        VCFMethods.samtools_index(
            strain_sorted_bam_dict=self.strain_sorted_bam_dict,
            strain_name_dict=self.strain_name_dict,
            threads=self.threads,
            logfile=self.logfile
        )
        logging.info('Extracting unmapped reads')
        strain_unmapped_reads_dict = VCFMethods.extract_unmapped_reads(
            strain_sorted_bam_dict=self.strain_sorted_bam_dict,
            strain_name_dict=self.strain_name_dict,
            threads=self.threads,
            logfile=self.logfile
        )
        logging.info('Attempting to assemble unmapped reads with SKESA')
        strain_skesa_output_fasta_dict = VCFMethods.assemble_unmapped_reads(
            strain_unmapped_reads_dict=strain_unmapped_reads_dict,
            strain_name_dict=self.strain_name_dict,
            threads=self.threads,
            logfile=self.logfile
        )
        logging.debug('SKESA assemblies: \n{files}'.format(
            files='\n'.join(['{strain_name}: {assembly}'.format(strain_name=sn, assembly=af)
                             for sn, af in strain_skesa_output_fasta_dict.items()])))
        logging.info('Running Quast on SKESA assemblies')
        quast_report_dict = VCFMethods \
            .quast(
                strain_skesa_output_fasta_dict=strain_skesa_output_fasta_dict,
                strain_unmapped_reads_dict=strain_unmapped_reads_dict,
                strain_sorted_bam_dict=self.strain_sorted_bam_dict,
                threads=self.threads,
                logfile=self.logfile
            )
        VCFMethods.parse_quast_report(
            quast_report_dict=quast_report_dict,
            summary_path=self.summary_path
        )

    def ref_file(self):
        """
        Populate dictionaries with appropriate strain name: reference file details
        reference_strain_dict: Dictionary of strain name: absolute path to symlinked reference genome
        """
        self.ref_fasta = glob(os.path.join(self.ref_path, '*.fasta'))[0]
        self.ref_strain = os.path.basename(os.path.splitext(self.ref_fasta)[0])
        for strain_name in self.strain_name_dict:
            self.reference_strain_dict[strain_name] = self.ref_fasta
            self.strain_consolidated_ref_dict[strain_name] = self.ref_strain
            self.reference_strain_dict[self.ref_strain] = self.ref_fasta

    def snp_calling(self):
        """
        Prep files for SNP calling. Use deepvariant to call SNPs. Parse the outputs from deepvariant
        """
        logging.info('Preparing files for SNP calling with deepvariant make_examples')
        self.strain_vcf_dict = VCFMethods.deepvariant_run_container(
            strain_sorted_bam_dict=self.strain_sorted_bam_dict,
            strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
            strain_name_dict=self.strain_name_dict,
            gpu=self.gpu,
            container_platform=self.container_platform,
            home=self.home,
            working_path=self.working_path,
            version=self.deepvariant_version,
            platform=self.platform,
            threads=self.threads,
            logfile=self.logfile
        )
        # strain_examples_dict, strain_variant_path_dict, strain_gvcf_tfrecords_dict = \
        #     VCFMethods.deepvariant_make_examples(strain_sorted_bam_dict=self.strain_sorted_bam_dict,
        #                                          strain_name_dict=self.strain_name_dict,
        #                                          strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
        #                                          vcf_path=os.path.join(self.seq_path, 'vcf_files'),
        #                                          home=self.home,
        #                                          logfile=self.logfile,
        #                                          threads=self.threads,
        #                                          working_path=self.working_path,
        #                                          deepvariant_version=self.deepvariant_version)
        # logging.info('Calling variants with deepvariant call_variants')
        # strain_call_variants_dict = \
        #     VCFMethods.deepvariant_call_variants(strain_variant_path_dict=strain_variant_path_dict,
        #                                          strain_name_dict=self.strain_name_dict,
        #                                          vcf_path=os.path.join(self.seq_path, 'vcf_files'),
        #                                          home=self.home,
        #                                          threads=self.threads,
        #                                          logfile=self.logfile,
        #                                          working_path=self.working_path,
        #                                          deepvariant_version=self.deepvariant_version,
        #                                          variant_caller='deepvariant')
        # logging.info('Creating VCF files with deepvariant postprocess_variants')
        # self.strain_vcf_dict = \
        #     VCFMethods.deepvariant_postprocess_variants_multiprocessing(
        #         strain_call_variants_dict=strain_call_variants_dict,
        #         strain_variant_path_dict=strain_variant_path_dict,
        #         strain_name_dict=self.strain_name_dict,
        #         strain_reference_abs_path_dict=self.strain_reference_abs_path_dict,
        #         strain_gvcf_tfrecords_dict=strain_gvcf_tfrecords_dict,
        #         vcf_path=os.path.join(self.seq_path, 'vcf_files'),
        #         home=self.home,
        #         logfile=self.logfile,
        #         deepvariant_version=self.deepvariant_version,
        #         threads=self.threads)
        # logging.info('Processing samples with deepvariant')
        # logging.info('Copying gVCF files to common folder')
        VCFMethods.copy_vcf_files(
            strain_vcf_dict=self.strain_vcf_dict,
            vcf_path=os.path.join(self.seq_path, 'vcf_files')
        )

    def load_snps(self):
        logging.info('Parsing gVCF files')
        self.strain_parsed_vcf_dict, self.strain_best_ref_dict, self.strain_best_ref_set_dict = \
            TreeMethods.load_vcf(strain_vcf_dict=self.strain_vcf_dict)
        if self.debug:
            logging.info('Parsed gVCF summaries:')
            pass_dict, insertion_dict, deletion_dict = \
                TreeMethods.summarise_gvcf_outputs(strain_parsed_vcf_dict=self.strain_parsed_vcf_dict)

            logging.info('SNP bases: \n{results}'.format(
                results='\n'.join(['{strain_name}: {pass_filter}'.format(strain_name=sn, pass_filter=ps)
                                   for sn, ps in pass_dict.items()])))
            logging.info('Inserted bases: \n{results}'.format(
                results='\n'.join(['{strain_name}: {insertion_calls}'.format(strain_name=sn, insertion_calls=ic)
                                   for sn, ic in insertion_dict.items()])))
            logging.info('Deleted bases: \n{results}'.format(
                results='\n'.join(['{strain_name}: {deletion_calls}'.format(strain_name=sn, deletion_calls=dc)
                                   for sn, dc in deletion_dict.items()])))
        logging.info('Loading SNP positions')
        consolidated_ref_snp_positions, strain_snp_positions, self.ref_snp_positions = \
            TreeMethods.load_gvcf_snp_positions(
                strain_parsed_vcf_dict=self.strain_parsed_vcf_dict,
                strain_consolidated_ref_dict=self.strain_consolidated_ref_dict
            )
        group_positions_set, self.strain_groups, self.strain_species_dict = \
            TreeMethods.group_strains(strain_snp_positions=strain_snp_positions)
        if self.debug:
            logging.info('Number of SNPs per contig:')
            for species_code, group_dict in group_positions_set.items():
                for group, ref_dict in group_dict.items():
                    for ref_chrom, pos_set in ref_dict.items():
                        if pos_set:
                            print(ref_chrom, len(pos_set))
        logging.info("Performing SNP density filtering")
        filtered_group_positions = TreeMethods.density_filter_snps(
            group_positions_set=group_positions_set,
            threshold=0
        )
        if self.debug:
            logging.info('Number of SNPs per contig following density filtering:')
            for species_code, group_dict in filtered_group_positions.items():
                for group, ref_dict in group_dict.items():
                    for ref_chrom, pos_set in ref_dict.items():
                        if pos_set:
                            print(ref_chrom, len(pos_set), sorted(list(pos_set)))
        logging.info('Masking low complexity and repeat regions in reference genomes')
        coords_dict = TreeMethods.mask_ref_genome(
            reference_strain_dict=self.reference_strain_dict,
            logfile=self.logfile
        )
        logging.info('Extracting coordinates to mask')
        mask_pos_dict = TreeMethods.determine_coordinates(
            strain_groups=self.strain_groups,
            coords_dict=coords_dict
        )
        if self.mask_file:
            logging.info(f'Loading masked regions from supplied mask file {self.mask_file}')
            supplied_mask_pos_dict = TreeMethods.load_supplied_mask(
                strain_groups=self.strain_groups,
                maskfile=self.mask_file
            )
        else:
            supplied_mask_pos_dict = dict()
        logging.info('Filtering SNPs in masked regions')
        filtered_masked_group_positions, filter_reasons = TreeMethods \
            .filter_masked_snp_positions(
                group_positions_set=group_positions_set,
                filtered_group_positions=filtered_group_positions,
                mask_pos_dict=mask_pos_dict,
                supplied_mask_pos_dict=supplied_mask_pos_dict
            )
        if self.debug:
            logging.info('Number of SNPs per contig following masking:')
            for species_code, group_dict in filtered_masked_group_positions.items():
                for group, ref_dict in group_dict.items():
                    for ref_chrom, pos_set in ref_dict.items():
                        if pos_set:
                            print(ref_chrom, len(pos_set), sorted(list(pos_set)))
        logging.info('Loading SNP sequences')
        self.group_strain_snp_sequence, self.species_group_best_ref = \
            TreeMethods.load_snp_sequence(
                strain_parsed_vcf_dict=self.strain_parsed_vcf_dict,
                strain_consolidated_ref_dict=self.strain_consolidated_ref_dict,
                group_positions_set=filtered_masked_group_positions,
                strain_groups=self.strain_groups,
                strain_species_dict=self.strain_species_dict,
                consolidated_ref_snp_positions=consolidated_ref_snp_positions,
                iupac=self.iupac
            )
        logging.info('Removing identical SNP positions from group')
        ident_group_positions = \
            TreeMethods.find_identical_calls(group_strain_snp_sequence=self.group_strain_snp_sequence)
        logging.info('Creating multi-FASTA files of core SNPs')
        group_folders, species_folders, self.group_fasta_dict = \
            TreeMethods.create_multifasta(
                group_strain_snp_sequence=self.group_strain_snp_sequence,
                fasta_path=self.fasta_path,
                group_positions_set=filtered_masked_group_positions,
                strain_parsed_vcf_dict=self.strain_parsed_vcf_dict,
                species_group_best_ref=self.species_group_best_ref,
                reference_strain_dict=self.reference_strain_dict,
                ident_group_positions=ident_group_positions,
                nested=False
            )
        if self.debug:
            logging.info('Multi-FASTA alignment output:')
            for species_code, group_dict in self.group_fasta_dict.items():
                for group, fasta_file in group_dict.items():
                    print(fasta_file)
        logging.info('Summarising SNPs')
        self.reference_strain_dict[self.ref_strain] = self.ref_fasta
        TreeMethods.snp_summary(
            group_strain_snp_sequence=self.group_strain_snp_sequence,
            species_group_best_ref=self.species_group_best_ref,
            reference_strain_dict=self.reference_strain_dict,
            group_positions_set=group_positions_set,
            filter_reasons=filter_reasons,
            strain_parsed_vcf_dict=self.strain_parsed_vcf_dict,
            filtered_group_positions=filtered_group_positions,
            mask_pos_dict=mask_pos_dict,
            supplied_mask_pos_dict=supplied_mask_pos_dict,
            ident_group_positions=ident_group_positions,
            summary_path=self.summary_path
        )

    def phylogenetic_trees(self):
        """
        Create, parse, and copy phylogenetic trees
        """
        logging.info('Creating phylogenetic trees with FastTree')
        species_group_trees = TreeMethods \
            .run_fasttree(
                group_fasta_dict=self.group_fasta_dict,
                strain_consolidated_ref_dict=self.strain_consolidated_ref_dict,
                strain_groups=self.strain_groups,
                logfile=self.logfile
            )
        logging.info('Parsing strain order from phylogenetic trees')
        self.species_group_order_dict = TreeMethods.parse_tree_order(species_group_trees=species_group_trees)
        logging.info('Copying phylogenetic trees to {tree_path}'.format(tree_path=self.tree_path))
        TreeMethods.copy_trees(
            species_group_trees=species_group_trees,
            tree_path=self.tree_path
        )

    def annotate_snps(self):
        """
        Load GenBank files, and annotate SNPs
        """
        logging.info('Creating GenBank file for {ref} as required'.format(ref=self.ref_strain))
        TreeMethods.prokka(reference_strain_dict=self.reference_strain_dict,
                           logfile=self.logfile)
        logging.info('Loading GenBank files for closest reference genomes')
        self.full_best_ref_gbk_dict = TreeMethods \
            .load_genbank_file_single(reference_strain_dict=self.reference_strain_dict)
        logging.info('Annotating SNPs')
        self.species_group_annotated_snps_dict = \
            TreeMethods.annotate_snps(
                group_strain_snp_sequence=self.group_strain_snp_sequence,
                full_best_ref_gbk_dict=self.full_best_ref_gbk_dict,
                strain_best_ref_set_dict=self.strain_best_ref_set_dict,
                ref_snp_positions=self.ref_snp_positions
            )

    def order_snps(self):
        """
        Order the SNPs based on prevalence and phylogeny
        """
        logging.info('Counting prevalence of SNPs')
        species_group_snp_num_dict = \
            TreeMethods.determine_snp_number(
                group_strain_snp_sequence=self.group_strain_snp_sequence,
                species_group_best_ref=self.species_group_best_ref)
        if self.debug:
            logging.info('SNP prevalence')
            for ref_chrom, pos_dict in species_group_snp_num_dict['species']['group'].items():
                if pos_dict:
                    print(ref_chrom, pos_dict)
        logging.info('Determining amino acid sequence at SNP locations')
        self.translated_snp_residue_dict, self.ref_translated_snp_residue_dict = \
            TreeMethods.determine_aa_sequence(
                group_strain_snp_sequence=self.group_strain_snp_sequence,
                species_group_best_ref=self.species_group_best_ref,
                strain_parsed_vcf_dict=self.strain_parsed_vcf_dict,
                species_group_annotated_snps_dict=self.species_group_annotated_snps_dict,
                reference_strain_dict=self.reference_strain_dict,
                species_group_snp_num_dict=species_group_snp_num_dict,
                iupac=self.iupac
            )
        logging.info('Creating SNP matrix')
        TreeMethods.create_snp_matrix(
            species_group_best_ref=self.species_group_best_ref,
            group_strain_snp_sequence=self.group_strain_snp_sequence,
            matrix_path=self.matrix_path
        )
        logging.info('Ranking SNPs based on prevalence')
        species_group_snp_rank, self.species_group_num_snps = \
            TreeMethods.rank_snps(species_group_snp_num_dict=species_group_snp_num_dict)
        if self.debug:
            logging.info('Ranked SNPs')
            for num_snps, ref_dict in sorted(species_group_snp_rank['species']['group'].items(), reverse=True):
                for ref_chrom, pos_dict in ref_dict.items():
                    print(num_snps, ref_chrom, pos_dict)
        logging.info('Sorting SNPs based on order of strains in phylogenetic trees')
        self.species_group_sorted_snps = \
            TreeMethods.sort_snps(
                species_group_order_dict=self.species_group_order_dict,
                species_group_snp_rank=species_group_snp_rank,
                species_group_best_ref=self.species_group_best_ref,
                group_strain_snp_sequence=self.group_strain_snp_sequence
            )
        if self.debug:
            logging.info('Sorted SNPs')
            for num_snps, ref_dict in self.species_group_sorted_snps['species']['group'].items():
                for ref_chrom, pos_dict in ref_dict.items():
                    print(num_snps, ref_chrom, pos_dict)

    def create_report(self):
        """
        Create the summary report of the analyses
        """
        logging.info('Creating summary tables')
        TreeMethods.create_summary_table(
            species_group_sorted_snps=self.species_group_sorted_snps,
            species_group_order_dict=self.species_group_order_dict,
            species_group_best_ref=self.species_group_best_ref,
            group_strain_snp_sequence=self.group_strain_snp_sequence,
            species_group_annotated_snps_dict=self.species_group_annotated_snps_dict,
            translated_snp_residue_dict=self.translated_snp_residue_dict,
            ref_translated_snp_residue_dict=self.ref_translated_snp_residue_dict,
            species_group_num_snps=self.species_group_num_snps,
            summary_path=self.summary_path,
            molecule='nt'
        )
        # Amino acid summary table
        TreeMethods.create_summary_table(
            species_group_sorted_snps=self.species_group_sorted_snps,
            species_group_order_dict=self.species_group_order_dict,
            species_group_best_ref=self.species_group_best_ref,
            group_strain_snp_sequence=self.group_strain_snp_sequence,
            species_group_annotated_snps_dict=self.species_group_annotated_snps_dict,
            translated_snp_residue_dict=self.translated_snp_residue_dict,
            ref_translated_snp_residue_dict=self.ref_translated_snp_residue_dict,
            species_group_num_snps=self.species_group_num_snps,
            summary_path=self.summary_path,
            molecule='aa'
        )

    def __init__(self, seq_path, ref_path, threads, working_path, mask_file, gpu, platform, container_platform, debug):
        # Determine the path in which the sequence files are located. Allow for ~ expansion
        if seq_path.startswith('~'):
            self.seq_path = os.path.abspath(os.path.expanduser(os.path.join(seq_path)))
        else:
            self.seq_path = os.path.abspath(os.path.join(seq_path))
        self.debug = debug
        SetupLogging(self.debug)
        logging.info('Welcome to {version}'.format(version=__version__))
        # Ensure that the path exists
        assert os.path.isdir(self.seq_path), 'Invalid path specified: {path}'.format(path=self.seq_path)
        logging.info('Supplied sequence path: \n{path}'.format(path=self.seq_path))
        # Initialise class variables
        self.threads = threads
        self.report_path = os.path.join(self.seq_path, 'reports')
        if ref_path.startswith('~'):
            self.ref_path = os.path.abspath(os.path.expanduser(os.path.join(ref_path)))
        else:
            self.ref_path = os.path.abspath(os.path.join(ref_path))
        # Ensure that the path exists
        assert os.path.isdir(self.ref_path), 'Invalid path specified: {path}'.format(path=self.ref_path)
        logging.info('Supplied reference path: \n{path}'.format(path=self.ref_path))
        # Determine if an additional volume
        if working_path:
            self.working_path = os.path.abspath(os.path.join(working_path))
        else:
            self.working_path = str()
        if mask_file:
            if mask_file.startswith('~'):
                self.mask_file = os.path.abspath(os.path.expanduser(os.path.join(mask_file)))
                assert os.path.isfile(self.mask_file), f'Cannot locate supplied mask file {mask_file}'
            elif mask_file.startswith('/'):
                self.mask_file = os.path.abspath(os.path.join(mask_file))
                assert os.path.isfile(self.mask_file), f'Cannot locate supplied mask file {mask_file}'
            elif '/' in mask_file:
                self.mask_file = os.path.join(os.path.dirname(self.ref_path), mask_file)
                assert os.path.isfile(self.mask_file), 'Cannot locate supplied mask file {mask_file} {path}'\
                    .format(
                        mask_file=mask_file,
                        path=os.path.join(os.path.dirname(self.ref_path), mask_file)
                )
            else:
                self.mask_file = os.path.join(self.ref_path, mask_file)
                assert os.path.isfile(self.mask_file), 'Cannot locate supplied mask file {mask_file} {ref_path}'.format(
                    mask_file=mask_file,
                    ref_path=self.ref_path
                )
        else:
            self.mask_file = str()
        self.platform = platform
        self.container_platform = container_platform
        self.home = str(Path.home())
        self.dependency_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'dependencies')
        self.fasta_path = os.path.join(self.seq_path, 'alignments')
        self.tree_path = os.path.join(self.seq_path, 'tree_files')
        self.summary_path = os.path.join(self.seq_path, 'summary_tables')
        self.matrix_path = os.path.join(self.seq_path, 'snv_matrix')
        self.logfile = os.path.join(self.seq_path, 'log')
        # Dictionary of degenerate IUPAC codes
        self.iupac = {
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
        self.ref_fasta = str()
        self.ref_strain = str()
        self.strain_name_dict = dict()
        self.strain_fastq_dict = dict()
        self.strain_consolidated_ref_dict = dict()
        self.strain_reference_abs_path_dict = dict()
        self.strain_reference_dep_path_dict = dict()
        self.strain_sorted_bam_dict = dict()
        self.reference_strain_dict = dict()
        self.strain_vcf_dict = dict()
        self.strain_parsed_vcf_dict = dict()
        self.strain_best_ref_dict = dict()
        self.strain_best_ref_set_dict = dict()
        self.ref_snp_positions = dict()
        self.group_strain_snp_sequence = dict()
        self.species_group_best_ref = dict()
        self.group_fasta_dict = dict()
        self.strain_groups = dict()
        self.strain_species_dict = dict()
        self.species_group_order_dict = dict()
        self.species_group_annotated_snps_dict = dict()
        self.translated_snp_residue_dict = dict()
        self.ref_translated_snp_residue_dict = dict()
        self.full_best_ref_gbk_dict = dict()
        self.species_group_num_snps = dict()
        self.species_group_sorted_snps = dict()
        self.gpu = gpu
        if self.gpu:
            self.deepvariant_version = '1.5.0-gpu'
        else:
            self.deepvariant_version = '1.5.0'


def main():
    # Parser for arguments
    parser = ArgumentParser(description='Finds SNPs between provided sequences and reference genome')
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=__version__)
    parser.add_argument(
        '-s', '--sequence_path',
        metavar='sequence_path',
        required=True,
        help='Path to folder containing sequencing reads')
    parser.add_argument(
        '-r', '--reference_path',
        metavar='reference_path',
        required=True,
        help='Provide the location of the folder containing the reference FASTA and GenBank files')
    parser.add_argument(
        '-t', '--threads',
        metavar='threads',
        type=int,
        default=multiprocessing.cpu_count() - 1,
        help='Number of threads. Default is the number of cores in the system - 1')
    parser.add_argument(
        '-d', '--debug',
        action='store_true',
        help='Enable debugging-level messages to be printed to the terminal')
    parser.add_argument(
        '-w', '--working_path',
        metavar='working_path',
        default=str(),
        help='If you are running these analyses anywhere other than your $HOME directory, you will '
             'need to provide the path to the drive e.g. /mnt/nas. This is necessary for the docker '
             'calls to deepvariant. An additional volume will be mounted in the docker container: '
             'e.g. -v /mnt/nas:/mnt/nas')
    parser.add_argument(
        '-m', '--mask_file',
        metavar='mask_file',
        type=str,
        default=str(),
        help='Supply a BED-formatted file with regions to mask. Generally, the format is: \n'
             'chrom\tchromStart\tchromEnd\n'
             'where chrom is the name of the reference chromosome; chromStart is the Start position '
             'of the feature in standard chromosomal coordinates (i.e. first base is 0); chromEnd is '
             'the End position of the feature in standard chromosomal coordinates')
    parser.add_argument(
        '-g', '--gpu',
        action='store_true',
        help='Enable this flag if your workstation has a GPU compatible with deepvariant. '
             'The program will use the deepvariant-gpu Docker image instead of the regular deepvariant '
             'image. Note that since I do not have a setup with a GPU, this is COMPLETELY UNTESTED!')
    parser.add_argument(
        '-platform', '--platform',
        metavar='platform',
        choices=['WGS', 'WES', 'PACBIO', 'ONT_R104', 'HYBRID_PACBIO_ILLUMINA'],
        default='WGS',
        help='Select the sequencing platform used to create the FASTQ reads. Default is WGS (Illumina)'
    )
    parser.add_argument(
        '-container_platform', '--container_platform',
        metavar='container_platform',
        choices=['docker', 'singularity'],
        default='docker',
        help='Select the container platform to use for running DeepVariant. Choices are "docker" and "singularity". '
             'Default is "docker"'
    )
    args = parser.parse_args()
    cowsnphr = COWSNPhR(
        seq_path=args.sequence_path,
        ref_path=args.reference_path,
        threads=args.threads,
        working_path=args.working_path,
        mask_file=args.mask_file,
        gpu=args.gpu,
        debug=args.debug,
        platform=args.platform,
        container_platform=args.container_platform
    )
    cowsnphr.main()
    logging.info('Analyses complete!')


if __name__ == '__main__':
    main()
