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
    def call_mash_sketch(strain_fastq_dict, strain_name_dict, threads, logfile):
        """

        :param strain_fastq_dict:
        :param strain_name_dict:
        :param threads:
        :param logfile:
        :return:
        """
        fastq_sketch_dict = dict()
        for strain_name, fastq_files in strain_fastq_dict.items():
            strain_folder = strain_name_dict[strain_name]
            fastq_sketch_no_ext = os.path.join(strain_folder, '{sn}_sketch'.format(sn=strain_name))
            fastq_sketch = fastq_sketch_no_ext + '.msh'
            file_list = os.path.join(strain_folder, 'fastq_files.txt')
            # Create a text file containing the FASTQ files to be used in the MASH analysis
            # with open(file_list, 'w') as fastq_list:
            #     fastq_list.write('\n'.join(fastq_files))
            mash_sketch_command = 'cat {fastq} | mash sketch -p {threads} -m 2 -r - -o {output_file}'\
                .format(fastq=' '.join(fastq_files),
                        threads=threads,
                        file_list=file_list,
                        output_file=fastq_sketch_no_ext)
            print(mash_sketch_command)
            if not os.path.isfile(fastq_sketch):
                out, err = run_subprocess(command=mash_sketch_command)
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=logfile,
                                 samplelog=os.path.join(strain_folder, 'log.out'),
                                 sampleerr=os.path.join(strain_folder, 'log.err'))
            #
            fastq_sketch_dict[strain_name] = fastq_sketch
        return fastq_sketch_dict

    @staticmethod
    def call_mash_dist(strain_fastq_dict, strain_name_dict, fastq_sketch_dict, ref_sketch_file, threads, logfile):
        strain_mash_outputs = dict()
        for strain_name in strain_fastq_dict:
            strain_folder = strain_name_dict[strain_name]
            fastq_sketch_file = fastq_sketch_dict[strain_name]
            out_tab = os.path.join(strain_folder, '{sn}_mash.tab'.format(sn=strain_name))
            mash_dist_command = 'mash dist -p {threads} -m 2 {ref_sketch_file} {fastq_sketch} | sort -gk2 > {out}' \
                .format(threads=threads,
                        ref_sketch_file=ref_sketch_file,
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