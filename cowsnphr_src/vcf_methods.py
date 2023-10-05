#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, \
    write_to_logfile
import multiprocessing
from glob import glob
import shutil
import gzip
import os
import re

__author__ = 'adamkoziol'


class VCFMethods(object):
    @staticmethod
    def file_list(path):
        """
        Create a list of all FASTQ files present in the supplied path. Accepts .fastq.gz, .fastq, .fq, and .fq.gz
        extensions only
        :param path: type STR: absolute path of folder containing FASTQ files
        :return fastq_files: sorted list of all FASTQ files present in the supplied path
        """
        # Use glob to find the acceptable extensions of FASTQ files in the supplied path
        fastq_files = glob(os.path.join(path, '*.fastq'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fastq.gz'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fq'))
        fastq_files = fastq_files + glob(os.path.join(path, '*.fq.gz'))
        # Sort the list of fastq files
        fastq_files = sorted(fastq_files)
        # Ensure that there are actually files present in the path
        assert fastq_files, 'Cannot find FASTQ files in the supplied path: {path}'.format(path=path)
        return fastq_files

    @staticmethod
    def strain_list(fastq_files):
        """
        Use the filer method to parse a list of absolute paths of FASTQ files to yield the name of the strain.
        e.g. /path/to/03-1057_S10_L001_R1_001.fastq.gz will return /path/to/03-1057:
        /path/to/03-1057_S10_L001_R1_001.fastq.gz
        :param fastq_files: type LIST: List of absolute of paired and/or unpaired FASTQ files
        :return strain_dict: dictionary of the absolute paths to the base strain name: FASTQ files
        """
        # As filer returns a set of the names, transform this set into a sorted list
        strain_dict = filer(
            filelist=fastq_files,
            returndict=True
        )
        return strain_dict

    @staticmethod
    def strain_namer(strain_folders):
        """
        Extract the base strain name from a list of absolute paths with the strain name (this list is usually created
        using the strain_list method above
        e.g. /path/to/03-1057 will yield 03-1057
        :param strain_folders: type iterable: List or dictionary of absolute paths and base strain names
        :return: strain_names: list of strain base names
        """
        # Initialise a dictionary to store the strain names
        strain_name_dict = dict()
        for strain_folder in strain_folders:
            # Extract the base name from the absolute path plus base name
            strain_name = os.path.basename(strain_folder)
            strain_name_dict[strain_name] = strain_folder
        return strain_name_dict

    @staticmethod
    def file_link(
            strain_folder_dict,
            strain_name_dict):
        """
        Create folders for each strain. Create relative symlinks to the original FASTQ files from within the folder
        :param strain_folder_dict: type DICT: Dictionary of strain folder path: FASTQ files
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :return: strain_fastq_dict: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        """
        #
        strain_fastq_dict = dict()
        for strain_name, strain_folder in strain_name_dict.items():
            # Create the strain folder path if required
            make_path(strain_folder)
            # Use the strain_folder value from the strain_name_dict as the key to extract the list of FASTQ files
            # associated with each strain
            for fastq_file in strain_folder_dict[strain_folder]:
                # Create relative symlinks between the original FASTQ files and the strain folder
                symlink_path = relative_symlink(
                    src_file=fastq_file,
                    output_dir=strain_folder,
                    export_output=True
                )
                # Add the absolute path of the symlink to the dictionary
                try:
                    strain_fastq_dict[strain_name].append(symlink_path)
                except KeyError:
                    strain_fastq_dict[strain_name] = [symlink_path]
        return strain_fastq_dict

    @staticmethod
    def run_reformat_reads(
            strain_fastq_dict,
            strain_name_dict,
            logfile):
        """
        Run reformat.sh from the BBMAP suite of tools. This will create histograms of the number of reads with a
        specific phred quality score, as well the number of reads of a specific length
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_qhist_dict: Dictionary of strain name: list of absolute paths to read set-specific
        quality histograms
        :return: strain_lhist_dict: Dictionary of strain name: list of absolute path to read set-specific
        length histograms
        """
        # Initialise dictionaries to store the absolute paths of the quality count, and the length distribution
        # histograms for each set of reads
        strain_qhist_dict = dict()
        strain_lhist_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            # Initialise the strain-specific lists in the dictionaries
            strain_qhist_dict[strain_name] = list()
            strain_lhist_dict[strain_name] = list()
            # Extract the absolute path of the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Initialise a counter to properly name the histogram files
            count = 1
            for fastq_file in fastq_files:
                # Set the absolute path of the quality and length histogram output file. Include the
                # strain name, as well as the file number count
                qual_histo = os.path.join(strain_folder, '{sn}_R{count}_qchist.csv'.format(sn=strain_name,
                                                                                           count=count))
                length_histo = os.path.join(strain_folder, '{sn}_R{count}_lhist.csv'.format(sn=strain_name,
                                                                                            count=count))
                # Create the system call to reformat.sh. qchist: count of bases with each quality value, lhist:
                # read length histogram
                reformat_cmd = 'reformat.sh in={fastq} qchist={qchist} lhist={lhist}' \
                    .format(fastq=fastq_file,
                            qchist=qual_histo,
                            lhist=length_histo)
                # Only run the analyses if the output file doesn't already exist
                if not os.path.isfile(length_histo):
                    out, err = run_subprocess(command=reformat_cmd)
                    # Write the stdout, and stderr to the main logfile, as well as to the strain-specific logs
                    write_to_logfile(
                        out=out,
                        err=err,
                        logfile=logfile,
                        samplelog=os.path.join(strain_folder, 'log.out'),
                        sampleerr=os.path.join(strain_folder, 'log.err')
                    )
                # Increase the file numbering count
                count += 1
                # Append the absolute path to the list of paths
                strain_qhist_dict[strain_name].append(qual_histo)
                strain_lhist_dict[strain_name].append(length_histo)
        return strain_qhist_dict, strain_lhist_dict

    @staticmethod
    def parse_quality_histogram(strain_qhist_dict):
        """
        Parse the quality histograms created by reformat.sh to calculate the average read quality as well as the
        percentage of reads with a Q-score greater than 30
        :param strain_qhist_dict: type: DICT: Dictionary of strain name: list of absolute paths to read set-specific
        quality histograms
        :return: strain_average_quality_dict: Dictionary of strain name: list of read set-specific average quality
        scores
        :return: strain_qual_over_thirty_dict: Dictionary of strain name: list of read set-specific percentage of
        reads with a Phred quality-score greater or equal to 30
        """
        # Initialise dictionaries to store the average read quality, and the percentage of reads with q-scores greater
        # than 30
        strain_average_quality_dict = dict()
        strain_qual_over_thirty_dict = dict()
        for strain_name, qual_histos in strain_qhist_dict.items():
            # Initialise the strain-specific list of outputs
            strain_average_quality_dict[strain_name] = list()
            strain_qual_over_thirty_dict[strain_name] = list()
            for qual_histo in sorted(qual_histos):
                # Read in the quality histogram
                with open(qual_histo, 'r') as histo:
                    # Skip the header line
                    next(histo)
                    # Initialise counts to store necessary integers
                    total_count = 0
                    total_read_quality = 0
                    qual_over_thirty_count = 0
                    for line in histo:
                        # Split each line on tabs into the quality score, the number of reads with that particular
                        # quality score, and the fraction of the total number reads these reads represent
                        quality, count, fraction = line.rstrip().split('\t')
                        # Manipulate the variables to integers
                        quality = int(quality)
                        count = int(count)
                        # Add read count * quality score to the cumulative read * quality score
                        total_read_quality += count * quality
                        # Add the current count to the total read count
                        total_count += count
                        # Determine if the read quality is >= 30. Add only those reads to the cumulative count
                        if quality >= 30:
                            qual_over_thirty_count += count
                    # Calculate the average quality: total quality count / total number of reads
                    average_qual = total_read_quality / total_count
                    # Calculate the % of reads with Q >= 30: number of reads with Q >= 30 / total number of reads
                    perc_reads_over_thirty = qual_over_thirty_count / total_count * 100
                    # Add the calculated values to the appropriate dictionaries
                    strain_average_quality_dict[strain_name].append(average_qual)
                    strain_qual_over_thirty_dict[strain_name].append(perc_reads_over_thirty)
        return strain_average_quality_dict, strain_qual_over_thirty_dict

    @staticmethod
    def parse_length_histograms(strain_lhist_dict):
        """
        Parse the length histogram created by reformat.sh to calculate the strain-specific average read length
        :param strain_lhist_dict: type DICT: Dictionary of strain name: list of absolute path to read set-specific
        length histograms
        :return: strain_avg_read_lengths: Dictionary of strain name: float of calculated strain-specific average
        read length
        """
        # Initialise a dictionary to store the strain-specific average read length
        strain_avg_read_lengths = dict()
        for strain_name, length_histos in strain_lhist_dict.items():
            # Initialise integers to store the total number of reads, as well as the total read count * length int
            total_count = 0
            total_count_length = 0
            # The average read quality is calculated on a per-sample, rather than a per-read set basis. So, the
            # variables are initialised outside the histo loop
            for length_histo in length_histos:
                with open(length_histo, 'r') as histo:
                    # Skip the header
                    next(histo)
                    for line in histo:
                        # Extract the read length and the number of reads of that length from the current line
                        length, count = line.rstrip().split('\t')
                        # Set the length and count to be integers
                        length = int(length)
                        count = int(count)
                        # Increment the total count by the current count
                        total_count += count
                        # Increment the total length * count variable by the current length * count
                        total_count_length += length * count
            # The average read length is calculated by dividing the total number of bases (length * count) by the
            # total number of reads
            avg_read_length = total_count_length / total_count
            # Populate the dictionary with the calculate average read length
            strain_avg_read_lengths[strain_name] = avg_read_length
        return strain_avg_read_lengths

    @staticmethod
    def find_fastq_size(strain_fastq_dict):
        """
        Use os.path.getsize to extract the size of the FASTQ files. Convert the value in bytes to megabytes
        :param strain_fastq_dict: type DICT: Dictionary of strain name: strain-specific FASTQ files
        :return: strain_fastq_size_dict: Dictionary of strain name: list of sizes of FASTQ files in megabytes
        """
        # Initialise a dictionary to store the strain-specific FASTQ file sizes
        strain_fastq_size_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            # Get the strain-specific list ready to be populated
            strain_fastq_size_dict[strain_name] = list()
            for fastq_file in fastq_files:
                # Use os.path.getsize to get the filesize in bytes. Convert this to megabytes by dividing this number
                # 1024*1024.0 (.0 is included, so that the divisor will be a float)
                file_size = os.path.getsize(fastq_file) / (1024 * 1024.0)
                strain_fastq_size_dict[strain_name].append(file_size)
        return strain_fastq_size_dict

    @staticmethod
    def call_mash_sketch(
            strain_fastq_dict,
            strain_name_dict,
            logfile):
        """
        Run MASH sketch on the provided FASTQ files
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :param logfile: type STR: Absolute path to logfile basename
        :return: fastq_sketch_dict: Dictionary of strain name: absolute path to MASH sketch file
        """
        # Initialise a dictionary to store the absolute path of the sketch file
        fastq_sketch_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            # Extract the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute path, and create the the mash folder
            mash_folder = os.path.join(strain_folder, 'mash')
            make_path(mash_folder)
            # Set the absolute paths of the sketch file with and without the .msh extension (used for calling MASH)
            fastq_sketch_no_ext = os.path.join(mash_folder, '{sn}_sketch'.format(sn=strain_name))
            fastq_sketch = fastq_sketch_no_ext + '.msh'
            # Create the system call - cat together the FASTQ files, and pipe them into MASH
            # -p requests the number of desired threads, -m sets the minimum copies of each k-mer required to pass
            # noise filter for reads to 2 (ignores single copy kmers), - indicates that MASH should use stdin
            # -o is the absolute path to the sketch output file
            mash_sketch_command = 'cat {fastq} | mash sketch -m 2 - -o {output_file}' \
                .format(fastq=' '.join(fastq_files),
                        output_file=fastq_sketch_no_ext)
            # Only make the system call if the output sketch file doesn't already exist
            if not os.path.isfile(fastq_sketch):
                out, err = run_subprocess(command=mash_sketch_command)
                # Write the stdout, and stderr to the main logfile, as well as to the strain-specific logs
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err')
                )
            # Populate the dictionary with the absolute path of the sketch file
            fastq_sketch_dict[strain_name] = fastq_sketch
        return fastq_sketch_dict

    @staticmethod
    def call_mash_dist(
            strain_fastq_dict,
            strain_name_dict,
            fastq_sketch_dict,
            ref_sketch_file,
            logfile):
        """
        Run a MASH dist of a pre-sketched set of FASTQ reads against the custom MASH sketch file of the reference
        genomes
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of absolute path(s) of FASTQ file(s)
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :param fastq_sketch_dict: type DICT: Dictionary of strain name: absolute path to MASH sketch file
        :param ref_sketch_file: type STR: Absolute path to the custom sketch file of reference sequences
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_mash_outputs: Dictionary of strain name: absolute path of MASH dist output table
        """
        # Initialise the dictionary to store the absolute path of the MASH dist output file
        strain_mash_outputs = dict()
        for strain_name in strain_fastq_dict:
            # Extract the absolute path of the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Extract the absolute path of the strain-specific sketch file
            fastq_sketch_file = fastq_sketch_dict[strain_name]
            # Set the absolute path of the MASH dist output table
            out_tab = os.path.join(strain_folder, 'mash', '{sn}_mash.tab'.format(sn=strain_name))
            # Create the system call: -p is the number of threads requested, -m sets the minimum copies of each k-mer
            # required to pass noise filter for reads to 2 (ignores single copy kmers). MASH outputs are piped to
            # the sort function, which sorts the data as follows: g: general numeric sort, K: keydef, 5: second column
            # (num matching hashes), r: reversed
            mash_dist_command = 'mash dist -m 2 {ref_sketch_file} {fastq_sketch} > {out}' \
                .format(ref_sketch_file=ref_sketch_file,
                        fastq_sketch=fastq_sketch_file,
                        out=out_tab)
            if not os.path.isfile(out_tab):
                out, err = run_subprocess(command=mash_dist_command)
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err')
                )
            strain_mash_outputs[strain_name] = out_tab
        return strain_mash_outputs

    @staticmethod
    def parse_mash_accession_species(mash_species_file):
        """
        Parse the reference genus accession: species code .csv file included in the mash dependencies path
        :param mash_species_file: type STR: Absolute path to file containing reference file name: species code
        :return: accession_species_dict: Dictionary of reference accession: species code
        """
        # Initialise a dictionary to store the species code
        accession_species_dict = dict()
        with open(mash_species_file, 'r') as species_file:
            for line in species_file:
                # Extract the accession and the species code pair from the line
                accession, species = line.rstrip().split(',')
                # Populate the dictionary with the accession: species pair
                accession_species_dict[accession] = species
        return accession_species_dict

    @staticmethod
    def mash_best_ref(
            mash_dist_dict,
            accession_species_dict,
            min_matches):
        """
        Parse the MASH dist output table to determine the closest reference sequence, as well as the total
        number of matching hashes the strain and that reference genome share
        :param mash_dist_dict: type DICT: Dictionary of strain name: absolute path of MASH dist output table
        :param accession_species_dict: type DICT: Dictionary of reference accession: species code
        :param min_matches: type INT: Minimum number of matching hashes required for a match to pass
        :return: strain_best_ref_dict: Dictionary of strain name: closest MASH-calculated reference genome
        :return: strain_ref_matches_dict: Dictionary of strain name: number of matching hashes between query and
        closest reference genome
        :return: strain_species_dict: Dictionary of strain name: species code
        """
        # Initialise dictionaries to store the strain-specific closest reference genome, number of matching hashes
        # between read sets and the reference genome, as well as the species code
        strain_best_ref_dict = dict()
        strain_ref_matches_dict = dict()
        strain_species_dict = dict()
        for strain_name, mash_dist_table in mash_dist_dict.items():
            # Create a variable to determine the best reference match
            best_matching_hashes = 0
            with open(mash_dist_table, 'r') as mash_dist:
                # Extract all the data included on each line of the table outputs
                for line in mash_dist:
                    # Split the line on tabs
                    best_ref, query_id, mash_distance, p_value, matching_hashes = line.rstrip().split('\t')
                    # Split the total of matching hashes from the total number of hashes
                    matching_hashes = int(matching_hashes.split('/')[0])
                    # Populate the dictionaries appropriately
                    if matching_hashes >= min_matches and matching_hashes > best_matching_hashes:
                        strain_best_ref_dict[strain_name] = best_ref
                        strain_ref_matches_dict[strain_name] = matching_hashes
                        strain_species_dict[strain_name] = accession_species_dict[best_ref]
                        # Update the best number of matching hashes with the current value
                        best_matching_hashes = matching_hashes
        return strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict

    @staticmethod
    def reference_folder(
            strain_best_ref_dict,
            dependency_path):
        """
        Create a dictionary of base strain name to the folder containing all the closest reference genome dependency
        files
        :param dependency_path: type STR: Absolute path to dependency folder
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: closest reference genome
        :return: reference_link_path_dict: Dictionary of strain name: relative path to symlinked reference genome
        :return: reference_link_dict: Dictionary of reference file.fasta: relative path to reference genome
        """
        # Initialise dictionaries
        reference_link_dict = dict()
        reference_link_path_dict = dict()
        # Read in the .csv file with reference file name: relative symbolic link information
        with open(os.path.join(dependency_path, 'reference_links.csv'), 'r') as reference_paths:
            for line in reference_paths:
                # Extract the link information
                reference, linked_file = line.rstrip().split(',')
                reference_link_dict[reference] = linked_file
        # Use the strain-specific best reference genome name to extract the relative symlink information
        for strain_name, best_ref in strain_best_ref_dict.items():
            reference_link_path_dict[strain_name] = reference_link_dict[best_ref]
        return reference_link_path_dict, reference_link_dict

    @staticmethod
    def index_ref_genome(
            reference_link_path_dict,
            dependency_path,
            logfile,
            reference_mapper):
        """
        Use bowtie2-build (or bwa index) to index the reference genomes
        :param reference_link_path_dict: type DICT: Dictionary of base strain name: reference folder path
        :param dependency_path: type STR: Absolute path to dependency folder
        :param logfile: type STR: Absolute path to logfile basename
        :param reference_mapper: type STR: Name of the reference mapping software to use. Choices are bwa and bowtie2
        :return: strain_mapper_index_dict: Dictionary of strain name: Absolute path to reference genome index
        :return: strain_reference_abs_path_dict: Dictionary of strain name: absolute path to reference file
        :return: strain_reference_dep_path_dict: Dictionary of strain name: absolute path to reference dependency folder
        """
        # Initialise a dictionary to store the absolute path to the bowtie2 index and reference genome
        strain_mapper_index_dict = dict()
        strain_reference_abs_path_dict = dict()
        strain_reference_dep_path_dict = dict()
        for strain_name, ref_link in reference_link_path_dict.items():
            # Set the absolute path, and strip off the file extension for use in the build call
            ref_abs_path = os.path.abspath(os.path.join(dependency_path, ref_link))
            base_name = os.path.splitext(ref_abs_path)[0]
            if reference_mapper == 'bowtie2':
                build_cmd = 'bowtie2-build {ref_file} {base_name}'.format(ref_file=ref_abs_path,
                                                                          base_name=base_name)
                index_file = base_name + '.1.bt2'
                strain_mapper_index_dict[strain_name] = base_name
            else:
                build_cmd = 'bwa index {ref_file}'.format(ref_file=ref_abs_path)
                index_file = ref_abs_path + '.bwt'
                strain_mapper_index_dict[strain_name] = ref_abs_path
            # Only run the system call if the index files haven't already been created
            if not os.path.isfile(index_file):
                out, err = run_subprocess(build_cmd)
                # Write the stdout and stderr to the log files
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile)
            # Populate the dictionaries
            strain_reference_abs_path_dict[strain_name] = ref_abs_path
            strain_reference_dep_path_dict[strain_name] = os.path.dirname(ref_abs_path)
        return strain_mapper_index_dict, strain_reference_abs_path_dict, strain_reference_dep_path_dict

    @staticmethod
    def faidx_ref_genome(
            reference_link_path_dict,
            dependency_path,
            logfile):
        """
        Run samtools faidx on the reference file
        :param reference_link_path_dict: type DICT: Dictionary of base strain name: reference folder path
        :param dependency_path: type STR: Absolute path to dependency folder
        :param logfile: type STR: Absolute path to logfile basename
        """
        for strain_name, ref_link in reference_link_path_dict.items():
            # Set the absolute path, and strip off the file extension for use in the build call
            ref_fasta = os.path.abspath(os.path.join(dependency_path, ref_link))
            faidx_command = 'samtools faidx {ref}'.format(ref=ref_fasta)
            if not os.path.isfile(ref_fasta + '.fai'):
                out, err = run_subprocess(command=faidx_command)
                write_to_logfile(
                    out=f'{faidx_command}\n{out}',
                    err=err,
                    logfile=logfile)

    @staticmethod
    def map_ref_genome(
            strain_fastq_dict,
            strain_name_dict,
            strain_mapper_index_dict,
            threads,
            logfile,
            reference_mapper):
        """
        Create a sorted BAM file by mapping the strain-specific FASTQ reads against the closest reference genome with
        bowtie2, converting the SAM outputs from bowtie2 to BAM format with samtools view, and sorting the BAM file
        with samtools sort. The individual commands are piped together to prevent the creation of unnecessary
        intermediate files
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of FASTQ files
        :param strain_name_dict: type DICT: Dictionary of strain name: strain-specific working folder
        :param strain_mapper_index_dict: type DICT: Dictionary of strain name: Absolute path to reference strain index
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to logfile basename
        :param reference_mapper: type STR: Name of the reference mapping software to use. Choices are bwa and bowtie2
        :return: strain_sorted_bam_dict: Dictionary of strain name: absolute path to sorted BAM files
        """
        # Initialise a dictionary to store the absolute path of the sorted BAM files
        strain_sorted_bam_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            try:
                # Extract the required variables from the appropriate dictionaries
                strain_folder = strain_name_dict[strain_name]
                reference_index = strain_mapper_index_dict[strain_name]
                # Set the absolute path of the sorted BAM file
                sorted_bam = os.path.join(strain_folder, '{sn}_sorted.bam'.format(sn=strain_name))
                if reference_mapper == 'bowtie2':
                    # Compound mapping command: bowtie2 (with read groups enabled: --rg-id  and --rg flags)|
                    map_cmd = 'bowtie2 --rg-id {sn} --rg SM:{sn} --rg PL:ILLUMINA --rg PI:250 -x {ref_index} ' \
                              '-U {fastq} -p {threads}'.format(sn=strain_name,
                                                               ref_index=reference_index,
                                                               fastq=','.join(fastq_files),
                                                               threads=threads)
                else:
                    # bwa mem mapping. Set the read group header to include the sample name in the ID and SM fields
                    map_cmd = 'bwa mem -M -R \"@RG\\tID:{sn}\\tSM:{sn}\\tPL:ILLUMINA\\tPI:250\" -t {threads} ' \
                              '{abs_ref_link} {fastq}' \
                        .format(sn=strain_name,
                                fastq=' '.join(fastq_files),
                                threads=threads,
                                abs_ref_link=reference_index)
                # Add the SAM-BAM conversion, duplicate read removal, and sorting commands to the mapping command
                # samtools view (-h: include headers, -b: out BAM, -T: target file)
                # samtools rmdup to remove duplicate reads
                # samtools sort
                map_cmd += ' | samtools view -@ {threads} -h -bT {abs_ref_link} -' \
                           ' | samtools rmdup - -S -' \
                           ' | samtools sort - -@ {threads} -o {sorted_bam}' \
                    .format(threads=threads,
                            abs_ref_link=reference_index,
                            sorted_bam=sorted_bam)
                # Only run the system call if the sorted BAM file doesn't already exist
                if not os.path.isfile(sorted_bam):
                    out, err = run_subprocess(map_cmd)
                    # Write STDOUT and STDERR to the logfile
                    write_to_logfile(
                        out=out,
                        err=err,
                        logfile=logfile,
                        samplelog=os.path.join(strain_folder, 'log.out'),
                        sampleerr=os.path.join(strain_folder, 'log.err'))
                # Populate the dictionary with the absolute path to the sorted BAM file
                strain_sorted_bam_dict[strain_name] = sorted_bam
            except KeyError:
                pass
        return strain_sorted_bam_dict

    @staticmethod
    def extract_unmapped_reads(
            strain_sorted_bam_dict,
            strain_name_dict,
            threads,
            logfile):
        """
        Use samtools bam2fq to extract all unmapped reads from the sorted BAM file into a single FASTQ file
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path to sorted BAM file
        :param strain_name_dict: type DICT: Dictionary of strain name: strain-specific working directory
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to the logfile basename
        :return: strain_unmapped_reads_dict: Dictionary of strain name: absolute path to unmapped reads FASTQ file
        """
        # Initialise a dictionary to store the absolute path of the unmapped reads file
        strain_unmapped_reads_dict = dict()
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            # Extract the absolute path of the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute path of the unmapped reads FASTQ file
            unmapped_reads = os.path.join(strain_folder, '{sn}_unmapped.fastq.gz'.format(sn=strain_name))
            # Create the system call to samtools bam2fq. Use -f4 to specify unmapped reads. Pipe output to gzip
            unmapped_cmd = 'samtools bam2fq -@ {threads} -f4 {sorted_bam} | gzip > {unmapped_reads}' \
                .format(threads=threads,
                        sorted_bam=sorted_bam,
                        unmapped_reads=unmapped_reads)
            # Run the system call if the reads file does not exist
            if not os.path.isfile(unmapped_reads):
                out, err = run_subprocess(unmapped_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err'))
            strain_unmapped_reads_dict[strain_name] = unmapped_reads
        return strain_unmapped_reads_dict

    @staticmethod
    def assemble_unmapped_reads(
            strain_unmapped_reads_dict,
            strain_name_dict,
            threads,
            logfile):
        """
        Run SKESA to attempt to assemble any unmapped reads
        :param strain_unmapped_reads_dict: type DICT: Dictionary of strain name: absolute path to unmapped reads
        FASTQ file
        :param strain_name_dict: type DICT: Dictionary of strain name: strain-specific working directory
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to the logfile basename
        :return: strain_skesa_output_fasta_dict: Dictionary of strain name: absolute path to SKESA assembly
        """
        # Initialise a dictionary to store the absolute path of the assembly
        strain_skesa_output_fasta_dict = dict()
        for strain_name, unmapped_reads in strain_unmapped_reads_dict.items():
            # Extract the strain-specific working directory
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute path, and create the SKESA output directory
            skesa_output_dir = os.path.join(strain_folder, 'skesa')
            make_path(skesa_output_dir)
            # Set the absolute path of the contigs file
            skesa_assembly_file = os.path.join(skesa_output_dir, '{sn}_unmapped.fasta'.format(sn=strain_name))
            # Create the SKESA system call. Use a minimum contig size of 1000
            skesa_cmd = 'skesa --fastq {fastqfiles} --cores {threads} --use_paired_ends --min_contig 1000 ' \
                        '--vector_percent 1 --contigs_out {contigs}' \
                .format(fastqfiles=unmapped_reads,
                        threads=threads,
                        contigs=skesa_assembly_file)
            # Run the system call if the qualimap report does not exist
            if not os.path.isfile(skesa_assembly_file):
                out, err = run_subprocess(skesa_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the absolute path to the contigs (the file will exist, but may be empty)
            strain_skesa_output_fasta_dict[strain_name] = skesa_assembly_file
        return strain_skesa_output_fasta_dict

    @staticmethod
    def quast(
            strain_skesa_output_fasta_dict,
            strain_unmapped_reads_dict,
            strain_sorted_bam_dict,
            threads,
            logfile):
        """
        Run quast on the samples
        :param strain_skesa_output_fasta_dict: type DICT: Dictionary of strain name: absolute path to SKESA assembly
        :param strain_unmapped_reads_dict: type DICT: Dictionary of strain name: absolute path to unmapped FASTQ reads
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path to sorted BAM file
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to the logfile basename
        :return: quast_report_dict: Dictionary of strain name: absolute path to quast report
        """
        quast_report_dict = dict()
        for strain_name, assembly_file in strain_skesa_output_fasta_dict.items():
            output_dir = os.path.dirname(assembly_file)
            quast_report = os.path.join(output_dir, 'transposed_report.tsv')
            quast_report_dict[strain_name] = quast_report
            # Create the system call to quast. --debug is specified, as certain temporary files are either used
            # for downstream analyses (BAM file), or parsed (insert size estimation)
            cmd = 'quast --single {single} --ref-bam {bam} -t {threads} --k-mer-stats --circos ' \
                  '--conserved-genes-finding -o {outputdir} --debug {assembly}' \
                .format(single=strain_unmapped_reads_dict[strain_name],
                        bam=strain_sorted_bam_dict[strain_name],
                        threads=threads,
                        outputdir=os.path.dirname(assembly_file),
                        assembly=assembly_file
                        )
            # Run the quast system call if the final quast report doesn't already exist
            if not os.path.isfile(quast_report):
                out, err = run_subprocess(cmd)
                # Write the appropriate information to the logfile
                write_to_logfile(
                    out=f'{cmd}\n{out}',
                    err=err,
                    logfile=logfile)
        return quast_report_dict

    @staticmethod
    def parse_quast_report(
            quast_report_dict,
            summary_path):
        """
        Parse the quast reports for each assembly, and create a combined report in the supplied report path
        :param quast_report_dict: type DICT: Dictionary of strain name: absolute path to quast report
        :param summary_path: type STR: Absolute path to directory in which reports are to be created
        """
        # Initialise strings to store the header and body information from the quast reports
        header = str()
        body = str()
        for sample_name, quast_report in quast_report_dict.items():
            if os.path.isfile(quast_report):
                # Read in the report
                with open(quast_report, 'r') as report:
                    # Populate the header string if it doesn't already exist
                    if not header:
                        header = report.readline()
                    else:
                        _ = report.readline()
                    # Add the report data to the string
                    for line in report:
                        body += line
        # Create the report path as required
        make_path(summary_path)
        # Write the header and body strings to the combined report
        with open(os.path.join(summary_path, 'assembly_report.tsv'), 'w') as quast:
            quast.write(header)
            quast.write(body)

    @staticmethod
    def samtools_index(
            strain_sorted_bam_dict,
            strain_name_dict,
            threads,
            logfile):
        """
        Index the sorted BAM file with samtools index
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path to sorted BAM file
        :param strain_name_dict: type DICT: Dictionary of strain name: strain-specific working directory
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to the logfile basename
        """
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            # Extract the folder name from the dictionary
            strain_folder = strain_name_dict[strain_name]
            # Set the system call for the samtools index command
            index_cmd = 'samtools index -@ {threads} {sorted_bam}'.format(threads=threads,
                                                                          sorted_bam=sorted_bam)
            # Only run the command if the .bai index file does not exist
            if not os.path.isfile(sorted_bam + '.bai'):
                out, err = run_subprocess(index_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err')
                )

    @staticmethod
    def deepvariant_run_container(
            strain_sorted_bam_dict,
            strain_reference_abs_path_dict,
            strain_name_dict,
            gpu,
            container_platform,
            home,
            working_path,
            version,
            platform,
            threads,
            logfile):
        """
        Run all three deepvariant binaries (make examples, call variants, postprocess_variants) with a single command.
        Use a GPU-enabled image if required. Allow for either Docker or Singularity to be used
        :return:
        """
        strain_gvcf_tfrecords_dict = {}
        strain_vcf_dict = {}
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute path, and create the deepvariant working directory
            deepvariant_dir = os.path.join(strain_folder, 'deepvariant')
            os.makedirs(deepvariant_dir, exist_ok=True)
            gvcf_tfrecords = '{output}_gvcf'.format(output=os.path.join(deepvariant_dir, strain_name))
            strain_gvcf_tfrecords_dict[strain_name] = \
                '{gvcf_tfrecords}@{threads}.gz'.format(gvcf_tfrecords=gvcf_tfrecords,
                                                       threads=threads)
            ref_genome = strain_reference_abs_path_dict[strain_name]
            if container_platform == 'singularity':
                cmd = 'singularity run '
                if gpu:
                    cmd += '--nv '
                cmd += f'-B /usr/lib/locale/:/usr/lib/locale/ ' \
                       f'docker://google/deepvariant:{version} '
            else:

                # Create the string of volume to mount to the container. Add the working path if it has been provided
                volumes = f'{home}:{home}'.format(home=home) if not working_path \
                    else f'{home}:{home} -v {working_path}:{working_path}'
                cmd = f'docker run '
                if gpu:
                    cmd += f'--gpus 1 '
                cmd += f'-v {volumes} ' \
                       f'google/deepvariant:{version} '
            output_vcf = os.path.join(deepvariant_dir, f'{strain_name}.vcf')
            gvcf_file = os.path.join(deepvariant_dir, '{sn}.gvcf.gz'.format(sn=strain_name))
            log_dir = os.path.join(deepvariant_dir, 'logs')
            os.makedirs(log_dir, exist_ok=True)
            cmd += f'/opt/deepvariant/bin/run_deepvariant ' \
                   f'--model_type={platform} ' \
                   f'--ref={ref_genome} ' \
                   f'--reads={sorted_bam} ' \
                   f'--output_vcf={output_vcf} ' \
                   f'--output_gvcf={gvcf_file} ' \
                   f'--num_shards={threads} ' \
                   f'--logging_dir={log_dir}'
            # Update the dictionary with the path to the output file
            strain_vcf_dict[strain_name] = gvcf_file
            if not os.path.isfile(gvcf_file):
                out, err = run_subprocess(cmd)
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile)
        return strain_vcf_dict

    @staticmethod
    def deepvariant_make_examples(
            strain_sorted_bam_dict,
            strain_name_dict,
            strain_reference_abs_path_dict,
            vcf_path,
            home,
            threads,
            logfile,
            deepvariant_version,
            working_path=None):
        """
        Use make_examples from deepvariant to extract pileup images from sorted BAM files. Currently, deepvariant
        does not support Python 3, so it will be run in a Docker container. The make_examples command is multithreaded
        with parallel
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path to strain-specific
        sorted BAM file
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path to strain-specific working dir
        :param strain_reference_abs_path_dict: type DICT: Dictionary of strain name: absolute path to best reference
        genome
        :param vcf_path: type STR: Absolute path to folder in which all symlinks to .vcf files are to be created
        :param home: type STR: Absolute path to $HOME
        :param threads: type INT: Number of threads to use in the analyses
        :param logfile: type STR: Absolute path to logfile basename
        :param deepvariant_version: type STR: Version number of deepvariant docker image to use
        :param working_path: type STR: Absolute path to an additional volume to mount to docker container
        :return: strain_examples_dict: Dictionary of strain name: sorted list of deepvariant created example files
        :return: strain_variant_path: Dictionary of strain name: absolute path to deepvariant working dir
        :return: strain_gvcf_tfrecords_dict: Dictionary of strain name: absolute path to gVCF TFRecord file of
        Variant protocol buffers
        """
        # Initialise a dictionary to store the absolute paths of all the output example files
        strain_examples_dict = dict()
        strain_gvcf_tfrecords_dict = dict()
        strain_variant_path = dict()
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            # Extract the necessary variables
            strain_folder = strain_name_dict[strain_name]
            ref_genome = strain_reference_abs_path_dict[strain_name]
            # Set the absolute path, and create the deepvariant working directory
            deepvariant_dir = os.path.join(strain_folder, 'deepvariant')
            strain_variant_path[strain_name] = deepvariant_dir
            make_path(deepvariant_dir)
            # Create the absolute path to the base name for the outputs. This will be used by parallel to create
            # all the output files
            output_example = '{output}_tfrecord'.format(output=os.path.join(deepvariant_dir, strain_name))
            gvcf_tfrecords = '{output}_gvcf'.format(output=os.path.join(deepvariant_dir, strain_name))
            strain_gvcf_tfrecords_dict[strain_name] = \
                '{gvcf_tfrecords}@{threads}.gz'.format(gvcf_tfrecords=gvcf_tfrecords,
                                                       threads=threads)
            # Create the string of volume to mount to the container. Add the working path if it has been provided
            volumes = '{home}:{home}'.format(home=home) if not working_path \
                else '{home}:{home} -v {working_path}:{working_path}'.format(home=home,
                                                                             working_path=working_path)
            # Create the system call. This consists of multiple parts:
            # 1) seq: generates range of numbers from start to end e.g. 0 to threads - 1
            # 2) parallel: executes commands in parallel. Use the following arguments:
            # --halt 2: kill all jobs if any job fails, --joblog: logfile of completed jobs, --res: store output
            # in the supplied folder
            # 3) docker run: use the deepvariant docker image to process the samples with the following arguments:
            # --rm: remove the docker container when finished, -v: mount the $HOME directory as a volume in the
            # container, --task: enables multi-processing by splitting the input, and generating sharded output
            make_example_cmd = 'seq 0 {threads_minus_1} | ' \
                               'parallel --halt 2 --joblog {logfile} --res {deepvariant_dir} ' \
                               'docker run --rm -v {volumes} ' \
                               'google/deepvariant:{dvv} ' \
                               '/opt/deepvariant/bin/make_examples --mode calling --ref {ref} --reads {bam} ' \
                               '--examples {output_example}@{threads}.gz --gvcf {gvcf_tfrecords}@{threads}.gz ' \
                               '--task {{}}' \
                .format(threads_minus_1=threads - 1,
                        logfile=os.path.join(strain_folder, 'log.out'),
                        deepvariant_dir=deepvariant_dir,
                        dvv=deepvariant_version,
                        volumes=volumes,
                        ref=ref_genome,
                        bam=sorted_bam,
                        output_example=output_example,
                        threads=threads,
                        gvcf_tfrecords=gvcf_tfrecords)
            # Create a list of all the sharded output example files
            output_examples = glob(os.path.join(deepvariant_dir, '{strain_name}_tfrecord-*.gz'
                                                .format(strain_name=strain_name)))
            # Ensure that the final outputs don't already exist
            vcf_file_name = os.path.join(vcf_path, '{sn}.gvcf.gz'.format(sn=strain_name))
            if not os.path.isfile(vcf_file_name):
                # If there is only one thread, need to check to see if output_examples exists, and then check
                # the number of threads
                if not output_examples:
                    # If there are fewer output files than expected (when running multi-threaded) or only one thread,
                    # run the system call
                    if len(output_examples) < threads - 1 or threads == 1:
                        out, err = run_subprocess(make_example_cmd)
                        write_to_logfile(
                            out=out,
                            err=err,
                            logfile=logfile)
            # Populate the dictionary with the sorted list of all the sharded files
            strain_examples_dict[strain_name] = \
                sorted(glob(os.path.join(deepvariant_dir, '{strain_name}_tfrecord-*.gz'
                                         .format(strain_name=strain_name))))
        return strain_examples_dict, strain_variant_path, strain_gvcf_tfrecords_dict

    @staticmethod
    def deepvariant_call_variants(
            strain_variant_path_dict,
            strain_name_dict,
            vcf_path,
            home,
            threads,
            logfile,
            variant_caller,
            deepvariant_version,
            working_path=None):
        """
        Perform variant calling. Process deepvariant examples files with call_variant
        :param strain_variant_path_dict: type DICT: Dictionary of strain name: absolute path to deepvariant output dir
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path to strain-specific working dir
        :param vcf_path: type STR: Absolute path to folder in which all symlinks to .vcf files are to be created
        :param home: type STR: Absolute path to $HOME
        :param threads: type INT: Number of threads used in the analyses
        :param logfile: type STR: Absolute path to logfile basename
        :param variant_caller: type STR: Variant calling software to use. Choices are deepvariant, and deepvariant-gpu
        (has GPU support)
        :param deepvariant_version: type STR: Version number of deepvariant docker image to use
        :param working_path: type STR: Absolute path to an additional volume to mount to docker container
        :return: strain_call_variants_dict: Dictionary of strain name: absolute path to deepvariant call_variants
        outputs
        """
        # Initialise a dictionary to store the absolute path of the call_variants output file
        strain_call_variants_dict = dict()
        for strain_name in strain_variant_path_dict:
            # Extract the necessary variables from dictionaries
            deepvariant_dir = strain_variant_path_dict[strain_name]
            strain_folder = strain_name_dict[strain_name]
            # Set the absolute path of the call_variants output file
            call_variants_output = os.path.join(deepvariant_dir, '{sn}_call_variants_output_tfrecord.gz'
                                                .format(sn=strain_name))
            # Update the dictionary with the file path
            strain_call_variants_dict[strain_name] = call_variants_output
            output_example = '{output}_tfrecord'.format(output=os.path.join(deepvariant_dir, strain_name))
            # Create the string of volume to mount to the container. Add the working path if it has been provided
            volumes = '{home}:{home}'.format(home=home) if not working_path \
                else '{home}:{home} -v {working_path}:{working_path}'.format(home=home,
                                                                             working_path=working_path)
            # Set the absolute path to the deepvariant model to be used by call_variants
            model = '/opt/models/wgs/model.ckpt'
            # Set the system call.
            if variant_caller == 'deepvariant-gpu':
                #  Use nvidia-docker to run call_variants with GPU support.
                call_variants_cmd = 'nvidia-docker run --rm -v {volumes} ' \
                                    'google/deepvariant-gpu:{dvv} ' \
                                    '/opt/deepvariant/bin/call_variants --outfile {output_file} ' \
                                    '--examples {output_example}@{threads}.gz --checkpoint {model} ' \
                    .format(volumes=volumes,
                            dvv=deepvariant_version,
                            output_file=call_variants_output,
                            output_example=output_example,
                            threads=threads,
                            model=model)
            else:
                # Use docker to run call_variants.
                call_variants_cmd = 'docker run --rm -v {volumes} google/deepvariant:{dvv} ' \
                                    '/opt/deepvariant/bin/call_variants --outfile {output_file} ' \
                                    '--examples {output_example}@{threads}.gz --checkpoint {model}' \
                    .format(volumes=volumes,
                            dvv=deepvariant_version,
                            output_file=call_variants_output,
                            output_example=output_example,
                            threads=threads,
                            model=model)
            # Run the system call if the output file and the final outputs don't already exist
            vcf_file_name = os.path.join(vcf_path, '{sn}.gvcf.gz'.format(sn=strain_name))
            if not os.path.isfile(vcf_file_name) and not os.path.isfile(call_variants_output):
                out, err = run_subprocess(call_variants_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(
                    out=out,
                    err=err,
                    logfile=logfile,
                    samplelog=os.path.join(strain_folder, 'log.out'),
                    sampleerr=os.path.join(strain_folder, 'log.err')
                )
        return strain_call_variants_dict

    @staticmethod
    def deepvariant_postprocess_variants_multiprocessing(
            strain_call_variants_dict,
            strain_variant_path_dict,
            strain_name_dict,
            strain_reference_abs_path_dict,
            strain_gvcf_tfrecords_dict,
            vcf_path,
            home,
            logfile,
            threads,
            deepvariant_version,
            working_path=None):
        """
        Create .gvcf.gz outputs
        :param strain_call_variants_dict: type DICT: Dictionary of strain name: absolute path to deepvariant
        call_variants outputs
        :param strain_variant_path_dict: type DICT: Dictionary of strain name: absolute path to deepvariant output dir
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path to strain-specific working dir
        :param strain_reference_abs_path_dict: type DICT: Dictionary of strain name: absolute path to best reference
        genome
        :param strain_gvcf_tfrecords_dict: type DICT: strain_call_variants_dict: Dictionary of strain name:
        absolute path to deepvariant call_variants outputs
        :param vcf_path: type STR: Absolute path to folder in which all symlinks to .vcf files are to be created
        :param home: type STR: Absolute path to $HOME
        :param logfile: type STR: Absolute path to logfile basename
        :param threads: type INT: Number of concurrent processes to spawn
        :param deepvariant_version: type STR: Version number of deepvariant docker image to use
        :param working_path: type STR: Absolute path to an additional volume to mount to docker container
        :return: strain_vcf_dict: Dictionary of strain name: absolute path to deepvariant output VCF file
        """
        # Initialise a dictionary to store the absolute path of the .vcf.gz output files
        strain_vcf_dict = dict()
        # Create a multiprocessing pool. Limit the number of processes to the number of threads
        p = multiprocessing.Pool(processes=threads)
        # Create a list of all the strain names
        strain_list = [strain_name for strain_name in strain_call_variants_dict]
        # Determine the number of strains present in the analyses
        list_length = len(strain_list)
        # Use multiprocessing.Pool.starmap to process the samples in parallel
        # Supply the list of strains, as well as a list the length of the number of strains of each required variable
        for vcf_dict in p.starmap(
                VCFMethods.deepvariant_postprocess_variants,
                zip(strain_list,
                    [strain_call_variants_dict] * list_length,
                    [strain_variant_path_dict] * list_length,
                    [strain_name_dict] * list_length,
                    [strain_reference_abs_path_dict] * list_length,
                    [strain_gvcf_tfrecords_dict] * list_length,
                    [vcf_path] * list_length,
                    [home] * list_length,
                    [logfile] * list_length,
                    [deepvariant_version] * list_length,
                    [working_path] * list_length)):
            # Update the dictionaries
            strain_vcf_dict.update(vcf_dict)
        # Close and join the pool
        p.close()
        p.join()
        return strain_vcf_dict

    @staticmethod
    def deepvariant_postprocess_variants(
            strain_name, strain_call_variants_dict,
            strain_variant_path_dict,
            strain_name_dict,
            strain_reference_abs_path_dict,
            strain_gvcf_tfrecords_dict,
            vcf_path,
            home,
            logfile,
            deepvariant_version,
            working_path=None):
        """
        Run the postprocess_variants script in the deepvariant Docker images. Creates global VCF output files
        :param strain_name: type STR: Name of strain currently being processed
        :param strain_call_variants_dict: type DICT: Dictionary of strain name: absolute path to deepvariant
        call_variants outputs
        :param strain_variant_path_dict: type DICT: Dictionary of strain name: absolute path to deepvariant output dir
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path to strain-specific working dir
        :param strain_reference_abs_path_dict: type DICT: Dictionary of strain name: absolute path to best reference
        genome
        :param strain_gvcf_tfrecords_dict: type DICT: strain_call_variants_dict: Dictionary of strain name:
        absolute path to deepvariant call_variants outputs
        :param vcf_path: type STR: Absolute path to folder in which all symlinks to .vcf files are to be created
        :param home: type STR: Absolute path to $HOME
        :param logfile: type STR: Absolute path to logfile basename
        :param deepvariant_version: type STR: Version number of deepvariant docker image to use
        :param working_path: type STR: Absolute path to an additional volume to mount to docker container
        :return: strain_vcf_dict: Dictionary of strain name: absolute path to .gvcf.gz output file
        """
        # Initialise a dictionary to store the absolute path of the .vcf.gz output files
        strain_vcf_dict = dict()
        # for strain_name, call_variants_output in strain_call_variants_dict.items():
        call_variants_output = strain_call_variants_dict[strain_name]
        # Extract the required variables from the dictionaries
        strain_folder = strain_name_dict[strain_name]
        deepvariant_dir = strain_variant_path_dict[strain_name]
        ref_genome = strain_reference_abs_path_dict[strain_name]
        gvcf_records = strain_gvcf_tfrecords_dict[strain_name]
        # Set the absolute path to the output file
        vcf_file = os.path.join(deepvariant_dir, '{sn}.vcf.gz'.format(sn=strain_name))
        gvcf_file = os.path.join(deepvariant_dir, '{sn}.gvcf.gz'.format(sn=strain_name))
        # Update the dictionary with the path to the output file
        strain_vcf_dict[strain_name] = gvcf_file
        # Create the string of volume to mount to the container. Add the working path if it has been provided
        volumes = '{home}:{home}'.format(home=home) if not working_path \
            else '{home}:{home} -v {working_path}:{working_path}'.format(home=home,
                                                                         working_path=working_path)
        # Set the system call. Use docker to run postprocess_variants
        postprocess_variants_cmd = 'docker run --rm -v {volumes} google/deepvariant:{dvv}' \
                                   ' /opt/deepvariant/bin/postprocess_variants --ref {ref} ' \
                                   '--nonvariant_site_tfrecord_path {gvcf_records} ' \
                                   '--infile {call_variants_output} --outfile {vcf_file} ' \
                                   '--gvcf_outfile {gvcf_file}' \
            .format(volumes=volumes,
                    dvv=deepvariant_version,
                    ref=ref_genome,
                    gvcf_records=gvcf_records,
                    call_variants_output=call_variants_output,
                    vcf_file=vcf_file,
                    gvcf_file=gvcf_file)
        # Ensure that the final outputs don't already exist
        vcf_file_name = os.path.join(vcf_path, '{sn}.gvcf.gz'.format(sn=strain_name))
        # Run the system call if the output file and the final output file don't exist
        if not os.path.isfile(vcf_file_name) and not os.path.isfile(gvcf_file):
            out, err = run_subprocess(postprocess_variants_cmd)
            # Write STDOUT and STDERR to the logfile
            write_to_logfile(
                out=out,
                err=err,
                logfile=logfile,
                samplelog=os.path.join(strain_folder, 'log.out'),
                sampleerr=os.path.join(strain_folder, 'log.err')
            )
        return strain_vcf_dict

    @staticmethod
    def parse_gvcf(strain_vcf_dict):
        """
        Determine the number of SNPs called that pass quality filters
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to gVCF file
        :return: strain_num_high_quality_snps_dict: Dictionary of strain name: number of high quality SNPs present
        """
        # Initialise a dictionary to store the number of high quality SNPs
        strain_num_high_quality_snps_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            strain_num_high_quality_snps_dict[strain_name] = int()
            # Ensure that the gVCF file was created
            if os.path.isfile(vcf_file):
                # Use gzip to open the compressed gVCF file
                with gzip.open(vcf_file, 'r') as gvcf:
                    for line in gvcf:
                        # Convert the line to a string from bytes
                        line = line.decode()
                        # Skip the header section
                        if line.startswith('#CHROM'):
                            for subline in gvcf:
                                subline = subline.decode()
                                # Split the line based on the columns
                                ref_genome, pos, id_stat, ref, alt_string, qual, filter_stat, info_string, \
                                    format_stat, strain = subline.split('\t')
                                # Initialise a string to hold the clean 'alt_string'
                                alt = str()
                                # Matches will have the following alt_string format: <*>. For SNP calls, the alt_string
                                # will have the following format: G,<*>. While insertions will look like: TGCC,<*>.
                                # Replace the <*>, and split on the comma
                                alt_split = alt_string.replace('<*>', '').split(',')
                                # If alt_split has a length greater than one e.g. not a match, which will look like
                                # [''], while a SNP and an insertion will be ['G', ''] and ['TGCC', ''], respectively
                                if len(alt_split) > 1:
                                    # Iterate through the list, and ensure that only the 'G' or the 'TGCC' are examined
                                    # rather than the empty ''
                                    for sub_alt in alt_split:
                                        if sub_alt:
                                            # Set the alt string as the 'G' or the 'TGCC'
                                            alt = sub_alt
                                # Filter the lines to find the high quality SNPs
                                # They must have a 'PASS' in the filter column. As well both the reference and query
                                # base must be of length one (a SNP rather than an indel), and the quality must pass
                                # the threshold
                                if filter_stat == 'PASS' and len(ref) == 1 and len(alt) == 1 and float(qual) > 14:
                                    # Add the passing SNP to the dictionary
                                    strain_num_high_quality_snps_dict[strain_name] += 1
        return strain_num_high_quality_snps_dict

    @staticmethod
    def parse_vcf(strain_vcf_dict):
        """
        Parse the .vcf file. Filter positions with QUAL < 150, and count the number of high
        quality SNPs (QUAL >= 150, and reference length = 1 (not an indel)). Also store indels and zero coverage
        regions. Overwrite the original file
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to .vcf file
        :return: strain_num_high_quality_snps_dict: Dictionary of strain name: number of high quality SNPs
        """
        # Initialise dictionaries to store the number of high quality SNPs, and the absolute path to the filtered
        # .vcf file
        strain_num_high_quality_snps_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            try:
                # Set the absolute path of the filtered .vcf file (replace the .vcf extension with _filtered.vcf)
                filtered_vcf = vcf_file.replace('.gvcf', '_filtered.gvcf')
                # Initialise a count for the number of high quality SNPs
                strain_num_high_quality_snps_dict[strain_name] = 0
                # # Iterate through all the records in the .vcf file
                with open(vcf_file, 'r') as unfiltered:
                    with open(filtered_vcf, 'w') as filtered:
                        for line in unfiltered:
                            # Add the VCF file header information to the filtered file
                            if line.startswith('#'):
                                filtered.write(line)
                            else:
                                # Split the line based on the columns
                                ref_genome, pos, id_stat, ref, alt_string, qual, filter_stat, info_string, \
                                    format_stat, strain = line.rstrip().split('\t')
                                # Find the depth entry. e.g. DP=11
                                depth_group = re.search('(DP=[0-9]+)', info_string)
                                # Split the depth matching group on '=' and convert the depth to an int
                                depth = int(str(depth_group.group()).split('=')[1])
                                # Store the zero coverage entries
                                if depth == 0:
                                    filtered.write(line)
                                else:
                                    #
                                    if len(ref) == 1:
                                        # Only store SNPs with a quality score greater or equal to 150
                                        if float(qual) >= 150:
                                            # Populate the dictionaries with the number of high quality SNPs
                                            strain_num_high_quality_snps_dict[strain_name] += 1
                                            filtered.write(line)
                                    # Store all indels
                                    else:
                                        filtered.write(line)
                # Remove the giant unfiltered file
                os.remove(vcf_file)
                # Rename the filtered file with the original file name
                shutil.move(src=filtered_vcf,
                            dst=vcf_file)
            except FileNotFoundError:
                pass
        return strain_num_high_quality_snps_dict

    @staticmethod
    def copy_vcf_files(
            strain_vcf_dict,
            vcf_path):
        """
        Create a folder with copies of the .vcf files
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to .vcf files
        :param vcf_path: type STR: Absolute path to folder in which all .gvcf.gz files are to be copied
        :return:
        """
        make_path(vcf_path)
        for strain_name, vcf_file in strain_vcf_dict.items():
            # Extract the file name of the gVCF file
            vcf_file_name = os.path.basename(vcf_file)
            if os.path.isfile(vcf_file):
                shutil.copyfile(src=vcf_file,
                                dst=os.path.join(vcf_path, vcf_file_name))
