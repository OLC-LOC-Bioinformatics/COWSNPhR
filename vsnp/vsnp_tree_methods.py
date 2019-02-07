#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, write_to_logfile
from Bio import SeqIO
from glob import glob
import xlsxwriter
import shutil
import vcf
import os

__author__ = 'adamkoziol'


class VSNPTreeMethods(object):

    @staticmethod
    def file_list(path):
        """
        Create a list of all VCF files present in the supplied path. Accepts .vcf extension only
        :param path: type STR: absolute path of folder containing VCF files
        :return vcf_files: sorted list of all VCF files present in the supplied path
        """
        # Use glob to find the acceptable extensions of VCF files in the supplied path
        vcf_files = glob(os.path.join(path, '*.vcf'))
        # Sort the list of VCF files
        vcf_files = sorted(vcf_files)
        # Ensure that there are actually files present in the path
        assert vcf_files, 'Cannot find VCF files in the supplied path: {path}'.format(path=path)
        return vcf_files

    @staticmethod
    def strain_list(vcf_files):
        """
        Parse a list of absolute paths of VCF files to yield the name of the strain.
        e.g. /path/to/03-1057.vcf will return a dictionary of /path/to/03-1057: path/to/03-1057.vcf
        :param vcf_files: type LIST: List of absolute of paths to vcf_files
        :return strain_dict: dictionary of the absolute paths to the base strain name: VCF files
        """
        # Initialise a dictionary to store the base name of the VCF file: absolute path to file
        strain_dict = dict()
        for vcf_file in vcf_files:
            # Split off the extension
            base_name = os.path.splitext(vcf_file)[0]
            # Remove _filtered and _zc from the base name
            base_name = base_name if not base_name.endswith('_filtered') else base_name.split('_filtered')[0]
            base_name = base_name if not base_name.endswith('_zc') else base_name.split('_zc')[0]
            # Populate the dictionary with the base name: absolute path to the file
            strain_dict[base_name] = vcf_file
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
        Create folders for each strain. Create relative symlinks to the original VCF files from within the folder
        :param strain_folder_dict: type DICT: Dictionary of strain folder path: VCF files
        :param strain_name_dict: type DICT: Dictionary of base strain name: strain folder path
        :return: strain_vcf_dict: Dictionary of strain name: list of absolute path(s) of VCF file
        """
        # Initialise a dictionary to store strain name: absolute path to VCF file links
        strain_vcf_dict = dict()
        for strain_name, strain_folder in strain_name_dict.items():
            # Create the strain folder path if required
            make_path(strain_folder)
            # Use the strain_folder value from the strain_name_dict as the key to extract the list of VCF files
            # associated with each strain
            vcf_file = strain_folder_dict[strain_folder]
            # Create relative symlinks between the original VCF files and the strain folder
            symlink_path = relative_symlink(src_file=vcf_file,
                                            output_dir=strain_folder,
                                            export_output=True)
            # Add the absolute path of the symlink to the dictionary
            try:
                strain_vcf_dict[strain_name].append(symlink_path)
            except KeyError:
                strain_vcf_dict[strain_name] = [symlink_path]
        return strain_vcf_dict
