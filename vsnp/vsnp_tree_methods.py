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
        :return: strain_vcf_dict: Dictionary of strain name: absolute path of VCF file
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
            strain_vcf_dict[strain_name] = symlink_path
        return strain_vcf_dict

    @staticmethod
    def parse_accession_species(ref_species_file):
        """
        Parse the reference genus accession: species code .csv file included in the mash dependencies path
        :param ref_species_file: type STR: Absolute path to file containing reference file name: species code
        :return: accession_species_dict: Dictionary of reference accession: species code
        """
        # Initialise a dictionary to store the species code
        accession_species_dict = dict()
        with open(ref_species_file, 'r') as species_file:
            for line in species_file:
                # Extract the accession and the species code pair from the line
                accession, species = line.rstrip().split(',')
                # Populate the dictionary with the accession: species pair
                accession_species_dict[accession] = species
        return accession_species_dict

    @staticmethod
    def load_vcf(strain_vcf_dict):
        """
        Using vcf.Reader(), load the VCF files. Store the Reader objects, as well as the extracted reference sequence,
        and its associated species code in dictionaries
        :param strain_vcf_dict: type DICT: Dictionary of strain name: list of absolute path to VCF file
        :return: strain_vcf_object_dict: Dictionary of strain name: VCF Reader object
        :return: strain_best_ref_dict: Dictionary of strain name: extracted reference sequence
        """
        # Initialise dictionaries to store the Reader objects, and the name of the reference sequence
        strain_vcf_object_dict = dict()
        strain_best_ref_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            # Create the PyVCF 'Reader' object from the .vcf file
            vcf_reader = vcf.Reader(open(vcf_file), 'r')
            # Iterate through all the records in the .vcf file
            for record in vcf_reader:
                # Extract the reference sequence stored as the .CHROM attribute
                strain_best_ref_dict[strain_name] = record.CHROM
                # Don't need to iterate through all the records once the reference sequence is known
                break
            # Add the Reader object to the dictionary
            strain_vcf_object_dict[strain_name] = vcf_reader
        return strain_vcf_object_dict, strain_best_ref_dict

    @staticmethod
    def determine_ref_species(strain_best_ref_dict, accession_species_dict):
        """
        Query the dictionary of reference file: species codes with the extracted strain-specific reference genome
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference sequence
        :param accession_species_dict: type DICT: Dictionary of reference accession: species code
        :return: strain_species_dict: Dictionary of strain name: species code
        :return: strain_best_ref_dict: Updated dictionary of strain name: best reference file
        """
        # Initialise a dictionary to store the extracted species code
        strain_species_dict = dict()
        for strain_name, best_ref in strain_best_ref_dict.items():
            for ref_file, species_code in accession_species_dict.items():
                # Remove any '.' from the best ref file
                # The match should work as follows: best_ref = NC_002945.4, ref_file = NC_002945v4. Strip off the
                # .4 from best_ref, and check to see if that string is present in ref_file
                if best_ref.split('.')[0] in ref_file:
                    # Populate the dictionary with the extracted species code
                    strain_species_dict[strain_name] = species_code
                    # Overwrite the strain_best_ref_dict for each strain with ref_file; this will ensure compatibility
                    # with downstream analyses
                    strain_best_ref_dict[strain_name] = ref_file
        return strain_species_dict, strain_best_ref_dict

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
        return reference_link_path_dict, reference_link_dict
