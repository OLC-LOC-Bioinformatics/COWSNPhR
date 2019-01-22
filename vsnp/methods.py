#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, write_to_logfile
from Bio import SeqIO
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
                file_size = os.path.getsize(fastq_file) / (1024*1024.0)
                strain_fastq_size_dict[strain_name].append(file_size)
        return strain_fastq_size_dict

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
            with open(mash_dist_table, 'r') as mash_dist:
                # Extract all the data included on each line of the table outputs
                best_ref, query_id, mash_distance, p_value, matching_hashes = mash_dist.readline().rstrip().split('\t')
            # Split the total of matching hashes from the total number of hashes
            matching_hashes = int(matching_hashes.split('/')[0])
            # Populate the dictionaries appropriately
            strain_best_ref_dict[strain_name] = best_ref
            strain_ref_matches_dict[strain_name] = matching_hashes
            strain_species_dict[strain_name] = accession_species_dict[best_ref]
        return strain_best_ref_dict, strain_ref_matches_dict, strain_species_dict

    @staticmethod
    def reference_folder(strain_best_ref_dict, dependency_path):
        """
        Create a dictionary of base strain name to the folder containing all the closest reference genome dependency
        files
        :param dependency_path: type STR: Absolute path to dependency folder
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: closest reference genome
        :return: reference_link_path_dict: Dictionary of strain name: relative path to symlinked reference genome
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
        return reference_link_path_dict

    @staticmethod
    def bowtie2_build(reference_link_path_dict, dependency_path, logfile):
        """
        Use bowtie2-build to index the reference genomes
        :param reference_link_path_dict: type DICT: Dictionary of base strain name: reference folder path
        :param dependency_path: type STR: Absolute path to dependency folder
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_bowtie2_index_dict: Dictionary of strain name: Absolute path to bowtie2 index
        :return: strain_reference_abs_path_dict: Dictionary of strain name: absolute path to reference file
        """
        # Initialise a dictionary to store the absolute path to the bowtie2 index and reference genome
        strain_bowtie2_index_dict = dict()
        strain_reference_abs_path_dict = dict()
        for strain_name, ref_link in reference_link_path_dict.items():
            # Set the absolute path, and strip off the file extension for use in the build call
            ref_abs_path = os.path.abspath(os.path.join(dependency_path, ref_link))
            base_name = os.path.splitext(ref_abs_path)[0]
            abs_ref_link = os.path.abspath(os.path.join(dependency_path, ref_link))
            build_cmd = 'bowtie2-build {ref_file} {base_name}'.format(ref_file=abs_ref_link,
                                                                      base_name=base_name)
            # Only run the system call if the index files haven't already been created
            if not os.path.isfile('{base_name}.1.bt2'.format(base_name=base_name)):
                out, err = run_subprocess(build_cmd)
                # Write the stdout and stderr to the log files
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile)
            # Populate the dictionaries
            strain_bowtie2_index_dict[strain_name] = base_name
            strain_reference_abs_path_dict[strain_name] = abs_ref_link
        return strain_bowtie2_index_dict, strain_reference_abs_path_dict

    @staticmethod
    def bowtie2_map(strain_fastq_dict, strain_name_dict, strain_bowtie2_index_dict, threads, logfile):
        """
        Create a sorted BAM file by mapping the strain-specific FASTQ reads against the closest reference genome with
        bowtie2, converting the SAM outputs from bowtie2 to BAM format with samtools view, and sorting the BAM file
        with samtools sort. The individual commands are piped together to prevent the creation of unnecessary
        intermediate files
        :param strain_fastq_dict: type DICT: Dictionary of strain name: list of FASTQ files
        :param strain_name_dict: type DICT: Dictionary of strain name: strain-specific working folder
        :param strain_bowtie2_index_dict: type DICT: Dictionary of strain name: Absolute path to bowtie2 index
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_sorted_bam_dict: Dictionary of strain name: absolute path to sorted BAM files
        """
        # Initialise a dictionary to store the absolute path of the sorted BAM files
        strain_sorted_bam_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            # Extract the required variables from the appropriate dictionaries
            strain_folder = strain_name_dict[strain_name]
            reference_index = strain_bowtie2_index_dict[strain_name]
            abs_ref_link = reference_index + '.fasta'
            # Set the absolute path of the sorted BAM file
            sorted_bam = os.path.join(strain_folder, '{sn}_sorted.bam'.format(sn=strain_name))
            # Compound mapping command: bowtie2 (with read groups enabled: --rg-id  and --rg flags)|
            # samtools view (-h: include headers, -b: out BAM, -T: target file)
            # samtools rmdup to remove duplicate reads
            # samtools sort
            map_cmd = 'bowtie2 --rg-id {sn} --rg SM:{sn} --rg PL:ILLUMINA --rg PI:250 -x {ref_index} ' \
                      '-U {fastq} -p {threads} | ' \
                      'samtools view -@ {threads} -h -bT {abs_ref_link} - | ' \
                      'samtools rmdup - -S - | ' \
                      'samtools sort - -@ {threads} -o {sorted_bam}'.format(sn=strain_name,
                                                                            ref_index=reference_index,
                                                                            fastq=','.join(fastq_files),
                                                                            threads=threads,
                                                                            abs_ref_link=abs_ref_link,
                                                                            sorted_bam=sorted_bam)
            # Only run the system call if the sorted BAM file doesn't already exist
            if not os.path.isfile(sorted_bam):
                out, err = run_subprocess(map_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the absolute path to the sorted BAM file
            strain_sorted_bam_dict[strain_name] = sorted_bam
        return strain_sorted_bam_dict

    @staticmethod
    def extract_unmapped_reads(strain_sorted_bam_dict, strain_name_dict, threads, logfile):
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
            unmapped_cmd = 'samtools bam2fq -@ {threads} -f4 {sorted_bam} | gzip > {unmapped_reads}'\
                .format(threads=threads,
                        sorted_bam=sorted_bam,
                        unmapped_reads=unmapped_reads)
            # Run the system call if the reads file does not exist
            if not os.path.isfile(unmapped_reads):
                out, err = run_subprocess(unmapped_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            strain_unmapped_reads_dict[strain_name] = unmapped_reads
        return strain_unmapped_reads_dict

    @staticmethod
    def assemble_unmapped_reads(strain_unmapped_reads_dict, strain_name_dict, threads, logfile):
        """
        Run SKESA to attempt to assembled any unmapped reads
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
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the absolute path to the contigs (the file will exist, but may be empty)
            strain_skesa_output_fasta_dict[strain_name] = skesa_assembly_file
        return strain_skesa_output_fasta_dict

    @staticmethod
    def assembly_stats(strain_skesa_output_fasta_dict):
        strain_unmapped_contigs_dict = dict()
        for strain_name, assembly_file in strain_skesa_output_fasta_dict.items():
            strain_unmapped_contigs_dict[strain_name] = 0
            if os.path.getsize(assembly_file) > 0:
                for _ in SeqIO.parse(assembly_file, 'fasta'):
                    strain_unmapped_contigs_dict[strain_name] += 1
        return strain_unmapped_contigs_dict

    @staticmethod
    def samtools_index(strain_sorted_bam_dict, strain_name_dict, threads, logfile):
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
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))

    @staticmethod
    def run_qualimap(strain_sorted_bam_dict, strain_name_dict, logfile):
        """
        Run qualimap on the sorted BAM files
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path of sorted BAM file
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path of strain-specific working folder
        :param logfile: type STR: Absolute path of logfile basename
        :return: strain_qualimap_report_dict: Dictionary of strain name: absolute path of qualimap report
        """
        # Initialise the dictionary to store the absolute path of the qualimap report
        strain_qualimap_report_dict = dict()
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            # Extract the absolute path of the working directory from the dictionary
            strain_folder = strain_name_dict[strain_name]
            # Set the path, and create a folder in which the qualimap results are to be stored
            qualimap_output_dir = os.path.join(strain_folder, 'qualimap')
            make_path(qualimap_output_dir)
            # Set the absolute path of the qualimap report to be parsed
            qualimap_report = os.path.join(qualimap_output_dir, 'genome_results.txt')
            # Create the qualimap system call
            qualimap_cmd = 'qualimap bamqc -bam {sorted_bam} -outdir {qualimap_out_dir}'\
                .format(sorted_bam=sorted_bam,
                        qualimap_out_dir=qualimap_output_dir)
            # Run the system call if the qualimap report does not exist
            if not os.path.isfile(qualimap_report):
                out, err = run_subprocess(qualimap_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the path to the report
            strain_qualimap_report_dict[strain_name] = qualimap_report
        return strain_qualimap_report_dict

    @staticmethod
    def parse_qualimap(strain_qualimap_report_dict):
        """
        Create a dictionary of the key: value pairs in the qualimap report
        :param strain_qualimap_report_dict: type DICT: Dictionary of strain name: absolute path of qualimap report
        :return: strain_qualimap_outputs_dict: Dictionary of strain name: dictionary of key: value from qualimap report
        """
        # Initialise a dictionary to store the qualimap outputs in dictionary format
        strain_qualimap_outputs_dict = dict()
        for strain_name, qualimap_report in strain_qualimap_report_dict.items():
            # Initialise a dictionary to hold the qualimap results - will be overwritten for each strain
            qualimap_dict = dict()
            with open(qualimap_report, 'r') as report:
                for line in report:
                    # Sanitise the keys and values using self.analyze
                    key, value = Methods.analyze(line)
                    # If the keys and values exist, enter them into the dictionary
                    if (key, value) != (None, None):
                        # Populate the dictionary with the sanitised key: value pair. Strip of the 'X' from the depth
                        # values
                        qualimap_dict[key] = value.rstrip('X')
            # Populate the dictionary with the sanitised outputs
            strain_qualimap_outputs_dict[strain_name] = qualimap_dict
        return strain_qualimap_outputs_dict

    @staticmethod
    def analyze(line):
        """
        Parse lines in qualimap reports. Split lines into key, value pairs, and sanitise lines, so that these pairs
        can be added to a dictionary
        :param line: type STR: Current line from a qualimap report
        :return: Sanitised key: value pair
        """
        # Split on ' = '
        if ' = ' in line:
            key, value = line.split(' = ')
            # Replace occurrences of: "number of ", "'", and " " with empty strings
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Remove commas
            value = value.replace(",", "").replace(" ", "").rstrip()
        # If '=' is absent, we are not interested in this line. Set the key and value to None
        else:
            key, value = None, None
        return key, value

    @staticmethod
    def reference_regions(strain_reference_abs_path_dict, logfile, size=100000):
        """
        Create a regions file to be used by freebayes-parallel. This file is the base pair position of the reference
        genome in regions of a defined size (in this case 100000). The file contains the contig name: bp position range
        NC_002945.4:0-100000
        NC_002945.4:100000-200000
        :param strain_reference_abs_path_dict: type DICT: Dictionary of strain name: absolute path of reference genome
        :param logfile: type STR: Absolute path to logfile basename
        :param size: type INT: Size of regions to create from reference genome
        :return: strain_ref_regions_dict: Dictionary of strain name: regions file
        """
        # Initialise a dictionary to store the absolute path of the regions file
        strain_ref_regions_dict = dict()
        for strain_name, ref_genome in strain_reference_abs_path_dict.items():
            # Set the absolute path of the regions file
            regions_file = ref_genome + '.regions'
            # Create the system call to fasta_generate_regions.py (part of the freebayes package)
            regions_cmd = 'fasta_generate_regions.py {ref_genome}.fai {size} > {regions_file}'\
                .format(ref_genome=ref_genome,
                        size=size,
                        regions_file=regions_file)
            # Run the system call if the regions file does not exist
            if not os.path.isfile(regions_file):
                out, err = run_subprocess(regions_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile)
            # Populate the dictionary with the path to the regions file
            strain_ref_regions_dict[strain_name] = regions_file
        return strain_ref_regions_dict

    @staticmethod
    def freebayes(strain_sorted_bam_dict, strain_name_dict, strain_reference_abs_path_dict, strain_ref_regions_dict,
                  threads, logfile):
        """
        Run freebayes-parallel on each sample
        :param strain_sorted_bam_dict: type DICT: Dictionary of strain name: absolute path of sorted BAM file
        :param strain_name_dict: type DICT: Dictionary of strain name: absolute path of strain-specific working folder
        :param strain_reference_abs_path_dict: type DICT: Dictionary of strain name: absolute path of reference genome
        :param strain_ref_regions_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        regions file
        :param threads: type INT: Number of threads to request for the analyses
        :param logfile: type STR: Absolute path to logfile basename
        :return: strain_vcf_dict: Dictionary of strain name: freebayes-created .vcf file
        """
        # Initialise a dictionary to store the absolute path to the .vcf files
        strain_vcf_dict = dict()
        for strain_name, sorted_bam in strain_sorted_bam_dict.items():
            # Extract the paths to the strain-specific working directory, the reference genome, and the regions file
            strain_folder = strain_name_dict[strain_name]
            ref_genome = strain_reference_abs_path_dict[strain_name]
            ref_regions_file = strain_ref_regions_dict[strain_name]
            # Set the absolute path to, and create the freebayes working directory
            freebayes_out_dir = os.path.join(strain_folder, 'freebayes')
            make_path(freebayes_out_dir)
            # Set the name of the output .vcf file
            freebayes_out_vcf = os.path.join(freebayes_out_dir, '{sn}.vcf'.format(sn=strain_name))
            # Create the system call to freebayes-parallel
            # Use the regions file to allow for parallelism
            freebayes_cmd = 'freebayes-parallel {ref_regions} {threads} -E -1 -e 1 -u --strict-vcf ' \
                            '-f {ref_genome} {sorted_bam} > {out_vcf}'\
                .format(ref_regions=ref_regions_file,
                        ref_genome=ref_genome,
                        threads=threads,
                        sorted_bam=sorted_bam,
                        out_vcf=freebayes_out_vcf)
            # Run the system call if the .vcf file does not exist
            if not os.path.isfile(freebayes_out_vcf):
                out, err = run_subprocess(freebayes_cmd)
                # Write STDOUT and STDERR to the logfile
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            # Populate the dictionary with the path to the .vcf file
            strain_vcf_dict[strain_name] = freebayes_out_vcf
        return strain_vcf_dict
