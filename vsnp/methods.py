#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, write_to_logfile
from biotools import mash
from glob import glob
import os
__author__ = 'adamkoziol'


class Methods(object):

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
        strain_dict = filer(filelist=fastq_files,
                            returndict=True)
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
    def file_link(strain_folder_dict, strain_name_dict):
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
                symlink_path = relative_symlink(src_file=fastq_file,
                                                output_dir=strain_folder,
                                                export_output=True)
                # Add the absolute path of the symlink to the dictionary
                try:
                    strain_fastq_dict[strain_name].append(symlink_path)
                except KeyError:
                    strain_fastq_dict[strain_name] = [symlink_path]
        return strain_fastq_dict

    @staticmethod
    def run_reformat_reads(strain_fastq_dict, strain_name_dict, logfile):
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
                reformat_cmd = 'reformat.sh in={fastq} qchist={qchist} lhist={lhist}'\
                    .format(fastq=fastq_file,
                            qchist=qual_histo,
                            lhist=length_histo)
                # Only run the analyses if the output file doesn't already exist
                if not os.path.isfile(length_histo):
                    out, err = run_subprocess(command=reformat_cmd)
                    # Write the stdout, and stderr to the main logfile, as well as to the strain-specific logs
                    write_to_logfile(out=out,
                                     err=err,
                                     logfile=logfile,
                                     samplelog=os.path.join(strain_folder, 'log.out'),
                                     sampleerr=os.path.join(strain_folder, 'log.err'))
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
            # variables are initialised outside of the histo loop
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
    def call_mash_sketch(strain_fastq_dict, strain_name_dict, logfile):
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
            # Set the absolute paths of the sketch file with and without the .msh extension (used for calling MASH)
            fastq_sketch_no_ext = os.path.join(strain_folder, '{sn}_sketch'.format(sn=strain_name))
            fastq_sketch = fastq_sketch_no_ext + '.msh'
            # Create the system call - cat together the FASTQ files, and pipe them into MASH
            # -p requests the number of desired threads, -m sets the minimum copies of each k-mer required to pass
            # noise filter for reads to 2 (ignores single copy kmers), - indicates that MASH should use stdin
            # -o is the absolute path to the sketch output file
            mash_sketch_command = 'cat {fastq} | mash sketch -m 2 - -o {output_file}'\
                .format(fastq=' '.join(fastq_files),
                        output_file=fastq_sketch_no_ext)
            # Only make the system call if the output sketch file doesn't already exist
            if not os.path.isfile(fastq_sketch):
                out, err = run_subprocess(command=mash_sketch_command)
                # Write the stdout, and stderr to the main logfile, as well as to the strain-specific logs
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the absolute path of the sketch file
            fastq_sketch_dict[strain_name] = fastq_sketch
        return fastq_sketch_dict

    @staticmethod
    def call_mash_dist(strain_fastq_dict, strain_name_dict, fastq_sketch_dict, ref_sketch_file, logfile):
        """
        Run a MASH dist of a pre-sketched set of FASTQ reads against the custom MASH sketch file of the reference
        genomes supported by vSNP
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
            out_tab = os.path.join(strain_folder, '{sn}_mash.tab'.format(sn=strain_name))
            # Create the system call: -p is the number of threads requested, -m sets the minimum copies of each k-mer
            # required to pass noise filter for reads to 2 (ignores single copy kmers). MASH outputs are piped to
            # the sort function, which sorts the data as follows: g: general numeric sort, K: keydef
            mash_dist_command = 'mash dist -m 2 {ref_sketch_file} {fastq_sketch} | sort -gk3 > {out}' \
                .format(ref_sketch_file=ref_sketch_file,
                        fastq_sketch=fastq_sketch_file,
                        out=out_tab)
            if not os.path.isfile(out_tab):
                out, err = run_subprocess(command=mash_dist_command)
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
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
    def mash_best_ref(mash_dist_dict, accession_species_dict):
        """
        Parse the MASH dist output table to determine the closest reference sequence, as well as the total
        number of matching hashes the strain and that reference genome share
        :param mash_dist_dict: type DICT: Dictionary of strain name: absolute path of MASH dist output table
        :param accession_species_dict: type DICT: Dictionary of reference accession: species code
        :return: strain_best_ref: Dictionary of strain name: closest MASH-calculated reference genome
        :return: strain_ref_matches: Dictionary of strain name: number of matching hashes between query and
        closest reference genome
        :return: strain_species: Dictionary of strain name: species code
        """
        # Initialise dictionaries to store the strain-specific closest reference genome, number of matching hashes
        # between read sets and the reference genome, as well as the species code
        strain_best_ref = dict()
        strain_ref_matches = dict()
        strain_species = dict()
        for strain_name, mash_dist_table in mash_dist_dict.items():
            with open(mash_dist_table, 'r') as mash_dist:
                # Extract all the data included on each line of the table outputs
                best_ref, query_id, mash_distance, p_value, matching_hashes = mash_dist.readline().rstrip().split('\t')
            # Split the total of matching hashes from the total number of hashes
            matching_hashes = int(matching_hashes.split('/')[0])
            # Populate the dictionaries appropriately
            strain_best_ref[strain_name] = best_ref
            strain_ref_matches[strain_name] = matching_hashes
            strain_species[strain_name] = accession_species_dict[best_ref]
        return strain_best_ref, strain_ref_matches, strain_species


