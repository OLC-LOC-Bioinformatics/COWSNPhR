#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, write_to_logfile
from Bio import SeqIO
from glob import glob
import xlsxwriter
import shutil
import pandas
import xlrd
import vcf
import os

__author__ = 'adamkoziol'


class VSNPTreeMethods(object):

    @staticmethod
    def file_list(file_path):
        """
        Create a list of all VCF files present in the supplied path. Accepts .vcf extension only
        :param file_path: type STR: absolute path of folder containing VCF files
        :return vcf_files: sorted list of all VCF files present in the supplied path
        """
        # Use glob to find the acceptable extensions of VCF files in the supplied path
        vcf_files = glob(os.path.join(file_path, '*.vcf'))
        # Sort the list of VCF files
        vcf_files = sorted(vcf_files)
        # Ensure that there are actually files present in the path
        assert vcf_files, 'Cannot find VCF files in the supplied path: {path}'.format(path=file_path)
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
        :return: strain_best_ref_dict: Dictionary of strain name: extracted reference genome name
        """
        # Initialise dictionaries to store the Reader objects, and the name of the reference sequence
        strain_vcf_object_dict = dict()
        strain_best_ref_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            # Create the PyVCF 'Reader' object from the .vcf file
            vcf_reader = vcf.Reader(open(vcf_file), 'r')
            # Add the Reader object to the dictionary
            strain_vcf_object_dict[strain_name] = list()
            # Iterate through all the records in the .vcf file
            for record in vcf_reader:
                # Extract the reference sequence stored as the .CHROM attribute
                strain_best_ref_dict[strain_name] = record.CHROM
                # Add the record to the dictionary
                strain_vcf_object_dict[strain_name].append(record)
        return strain_vcf_object_dict, strain_best_ref_dict

    @staticmethod
    def determine_ref_species(strain_best_ref_dict, accession_species_dict):
        """
        Query the dictionary of reference file: species codes with the extracted strain-specific reference genome
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param accession_species_dict: type DICT: Dictionary of reference accession: species code
        :return: strain_species_dict: Dictionary of strain name: species code
        :return: strain_best_ref_dict: Updated dictionary of strain name: best reference file
        """
        # Initialise a dictionary to store the extracted species code
        strain_species_dict = dict()
        strain_best_ref_fasta_dict = dict()
        for strain_name, best_ref in strain_best_ref_dict.items():
            for ref_file, species_code in accession_species_dict.items():
                # Remove any '.' from the best ref file
                # The match should work as follows: best_ref = NC_002945.4, ref_file = NC_002945v4. Strip off the
                # .4 from best_ref, and check to see if that string is present in ref_file
                if best_ref.split('.')[0] in ref_file:
                    # Populate the dictionary with the extracted species code
                    strain_species_dict[strain_name] = species_code
                    # Add the path to the reference file to the dictionary
                    strain_best_ref_fasta_dict[strain_name] = ref_file
        return strain_species_dict, strain_best_ref_fasta_dict

    @staticmethod
    def reference_folder(strain_best_ref_fasta_dict, dependency_path):
        """
        Create a dictionary of base strain name to the folder containing all the closest reference genome dependency
        files
        :param dependency_path: type STR: Absolute path to dependency folder
        :param strain_best_ref_fasta_dict: type DICT: Dictionary of strain name: path to closest reference genome FASTA
        file
        :return: reference_link_path_dict: Dictionary of strain name: relative path to reference genome dependency
        folder
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
        for strain_name, best_ref in strain_best_ref_fasta_dict.items():
            reference_link_path_dict[strain_name] = os.path.dirname(reference_link_dict[best_ref])
        return reference_link_path_dict, reference_link_dict

    @staticmethod
    def extract_defining_snps(reference_link_path_dict, strain_species_dict, dependency_path):
        """
        Load the Excel files containing species-specific groupings of defining SNPs
        :param reference_link_path_dict: type DICT: Dictionary of strain name: relative path to reference genome
        dependency folder
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :param dependency_path: type STR: Absolute path to dependencies
        :return: defining_snp_dict: Dictionary of species code: dictionary of grouping: reference genome: defining SNP
        """
        # Initialise a dictionary to store the species-specific groups of defining SNPs
        defining_snp_dict = dict()
        for strain_name, best_ref_path in reference_link_path_dict.items():
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            # Set the name of the Excel file containing the groups and their defining SNPs
            defining_snp_xlsx = os.path.join(dependency_path, best_ref_path, 'DefiningSNPsGroupDesignations.xlsx')
            if os.path.isfile(defining_snp_xlsx):
                # Only populate the dictionary once
                if species not in defining_snp_dict:
                    # Get the dictionary read to accept the species-specific group: SNP pair
                    defining_snp_dict[species] = dict()
                    # Use Pandas to read in the Excel file, and convert the Pandas object to a dictionary
                    snp_dict = pandas.read_excel(defining_snp_xlsx).to_dict()
                    # Iterate through the 'Grouping' column to extract the grouping name
                    for i, grouping in snp_dict['Grouping'].items():
                        # Use the iterator to extract the matching value from the 'Absolute position' column
                        try:
                            reference, position = str(snp_dict['Absolute position'][i]).split('-')
                            defining_snp_dict[species][grouping] = {reference: position}
                        # Ignore nan entries
                        except ValueError:
                            pass
                    # TB best reference af2122 has a second group: defining SNP column pair (parsed as Unnamed: 3 and
                    # Unnamed: 4 by Pandas. Add these pairs, too
                    try:
                        for i, grouping in snp_dict['Unnamed: 3'].items():
                            # Most columns are empty, and can be ignored
                            try:
                                reference, position = str(snp_dict['Unnamed: 4'][i]).split('-')
                                defining_snp_dict[species][grouping] = {reference: position}
                            except ValueError:
                                pass
                    except KeyError:
                        pass
        return defining_snp_dict

    @staticmethod
    def load_snp_positions(strain_vcf_object_dict, strain_species_dict, strain_best_ref_dict):
        """
        Parse the VCF files, and extract all the strain-specific SNP locations and reference sequence
        :param strain_vcf_object_dict: type DICT: Dictionary of strain name: VCF Reader object
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :return: ref_snp_positions: Dictionary of reference name: absolute position: reference base call
        :return: strain_snp_positions: Dictionary of strain name: all strain-specific SNP positions
        :return: strain_snp_sequence: Dictionary of strain name: SNP position: strain-specific sequence
        """
        # Initialise a dictionary to store all the SNP locations, and the reference sequence for each reference
        # genome
        ref_snp_positions = dict()
        strain_snp_positions = dict()
        strain_snp_sequence = dict()
        for strain_name, vcf_object in strain_vcf_object_dict.items():
            strain_snp_positions[strain_name] = list()
            strain_snp_sequence[strain_name] = dict()
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            # Extract the species code from the dictionary
            best_ref = strain_best_ref_dict[strain_name]
            # Iterate through all the records
            for record in vcf_object:
                # Initialise the reference genome key as required
                if best_ref not in ref_snp_positions:
                    ref_snp_positions[best_ref] = dict()
                try:
                    # Flu VCF files must be treated slightly differently
                    if species == 'flu':
                        # Use both AC=1 and AC=2 as valid position
                        if str(record.ALT[0]) != "None" and len(record.REF) == 1 and record.QUAL >= 150:
                            # Add the position, and the reference base to the dictionary
                            ref_snp_positions[best_ref][record.POS] = record.REF
                            strain_snp_positions[strain_name].append(record.POS)
                            strain_snp_sequence[strain_name][record.POS] = record.ALT[0]
                    else:
                        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 2 and len(
                                record.REF) == 1 and record.QUAL >= 150:
                            ref_snp_positions[best_ref][record.POS] = record.REF
                            strain_snp_positions[strain_name].append(record.POS)
                            strain_snp_sequence[strain_name][record.POS] = record.ALT[0]
                except KeyError:
                    pass
        return ref_snp_positions, strain_snp_positions, strain_snp_sequence

    @staticmethod
    def determine_groups(strain_snp_positions, defining_snp_dict):
        """
        Determine which defining SNPs are present in strains
        :param strain_snp_positions: type DICT: Dictionary of strain name: all strain-specific SNP positions
        :param defining_snp_dict: type DICT: Dictionary of species code: dictionary of grouping: reference genome:
        defining SNP
        :return: strain_groups: Dictionary of strain name: list of group(s) for which the strain contains the defining
        SNP
        """
        # Initialise a dictionary to store the list of groups to which the strain belongs
        strain_groups = dict()
        for strain_name, snp_positions in strain_snp_positions.items():
            strain_groups[strain_name] = list()
            # Unpack the dictionary
            for species, nested_dict in defining_snp_dict.items():
                for group, ref_snp_dict in nested_dict.items():
                    for ref, snp in ref_snp_dict.items():
                        # Inverted positions have a trailing '!', and therefore cannot be typecast to int
                        try:
                            position = int(snp)
                        except ValueError:
                            # Inverted
                            position = int(snp.rstrip('!'))
                        # If the position is in the list of strain-specific SNP positions, add the group to the list
                        if position in snp_positions:
                            strain_groups[strain_name].append(group)
        return strain_groups

    @staticmethod
    def load_filter_file(reference_link_path_dict, strain_best_ref_dict, dependency_path):
        """
        Load the Excel files containing curated lists of locations or ranges of locations in the reference genome
        that must be filtered prior to performing phylogenetic analyses
        :param reference_link_path_dict: type DICT: Dictionary of strain name: relative path to reference genome
        dependency folder
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param dependency_path: type STR: Absolute path to dependencies
        :return: filter_dict: Dictionary of reference file: group name: locations to filter
        """
        # Initialise a dictionary to store the locations to filter
        filter_dict = dict()
        for strain_name, best_ref_path in reference_link_path_dict.items():
            # Extract the name of the best reference from the dictionary e.g. NC_017250.1
            best_ref = strain_best_ref_dict[strain_name]
            # Set the name of the Excel file storing the regions to filter
            filter_file = os.path.join(dependency_path, best_ref_path, 'Filtered_Regions.xlsx')
            if os.path.isfile(filter_file):
                # Only load the file once per reference file
                if best_ref not in filter_dict:
                    # Open the file using xlrd
                    wb = xlrd.open_workbook(filter_file)
                    # Create a list of all the sheet names
                    sheets = wb.sheet_names()
                    # Iterate through all the sheets
                    for sheet in sheets:
                        # Initialise the dictionary with the sheet name - this will be the same as the best reference
                        # file name
                        filter_dict[sheet] = dict()
                        # Load each worksheet
                        ws = wb.sheet_by_name(sheet)
                        # Iterate through all the columns in the worksheet
                        for col_num in range(ws.ncols):
                            # Extract the group name from the header e.g. Bsuis1-All
                            group_name = ws.col_values(col_num)[0]
                            # Initialise the group name as a list in the nested dictionary
                            filter_dict[sheet][group_name] = list()
                            # Create a list of all the non-header entries in the column
                            entries = ws.col_values(col_num)[1:]
                            # Remove blank cells and typecast the entries to strings
                            entries = [str(entry) for entry in entries if entry]
                            for value in entries:
                                # Certain filtered SNPs are actually ranges
                                if "-" not in value:
                                    # Convert the value to an integer via a float
                                    value = int(float(value))
                                    # Add the position to the dictionary
                                    filter_dict[sheet][group_name].append(value)
                                elif "-" in value:
                                    # Split the range on '-' e.g. '524691-524833' becomes ['524691', '524833']
                                    values = value.split("-")
                                    # Iterate through all the position in the range of values
                                    for position in range(int(values[0]), int(values[1]) + 1):
                                        # Add each position in the range to the dictionary
                                        filter_dict[sheet][group_name].append(position)
        return filter_dict

    @staticmethod
    def filter_positions(strain_snp_positions, strain_groups, strain_best_ref_dict, filter_dict, strain_snp_sequence):
        """
        Use the filtered positions to filter SNPs found in the strains
        :param strain_snp_positions: type DICT: Dictionary of strain name: all strain-specific SNP positions
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param filter_dict: type DICT: Dictionary of reference file: group name: locations to filter
        :param strain_snp_sequence: type DICT: Dictionary of strain name: SNP position: strain-specific sequence
        :return: strain_filtered_sequences: Dictionary of strain name: SNP pos: SNP sequence
        """
        # Initialise a dictionary to store the location and sequence of SNPs that pass filter
        strain_filtered_sequences = dict()
        for strain_name, snp_positions in strain_snp_positions.items():
            strain_filtered_sequences[strain_name] = dict()
            # Extract the necessary variables from dictionaries
            groups = strain_groups[strain_name]
            best_ref = strain_best_ref_dict[strain_name]
            for group, positions in filter_dict[best_ref].items():
                # All strains of a particular species fall within the 'All' category
                if 'All' in group or group in groups:
                    # Initialise the dictionary with the group
                    strain_filtered_sequences[strain_name][group] = dict()
                    for snp_pos in snp_positions:
                        # Use the filter positions to remove unwanted positions
                        if snp_pos not in positions:
                            # Populate the dictionary with the position and the extracted SNP sequence from
                            # the dictionary
                            strain_filtered_sequences[strain_name][group][snp_pos] = \
                                strain_snp_sequence[strain_name][snp_pos]
        return strain_filtered_sequences

    @staticmethod
    def load_genbank_file(reference_link_path_dict, strain_best_ref_dict, dependency_path):
        """
        Use SeqIO to parse the best reference genome GenBank file for annotating SNP locations
        :param reference_link_path_dict: type DICT: Dictionary of strain name: relative path to reference genome
        dependency folder
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param dependency_path: type STR: Absolute path to dependencies
        :return: best_ref_gbk_dict: Dictionary of best ref: SeqIO parsed GenBank file-sourced records from closest
        reference genome
        """
        # Initialise a dictionary to store the SeqIO parsed GenBank files
        best_ref_gbk_dict = dict()
        for strain_name, best_ref_path in reference_link_path_dict.items():
            # Extract the species code from the dictionary
            best_ref = strain_best_ref_dict[strain_name]
            gbk_file = glob(os.path.join(
                dependency_path, best_ref_path, '{br}*.gbk'.format(br=os.path.splitext(best_ref)[0])))[0]
            # Only parse the file if it has not already been parsed
            if best_ref not in best_ref_gbk_dict:
                # Use SeqIO to first parse the GenBank file, and convert the parsed object to a dictionary
                gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
                # Add the GenBank dictionary to the best reference-specific dictionary
                best_ref_gbk_dict[best_ref] = gbk_dict
        return best_ref_gbk_dict

    # @staticmethod
    # def species_specific_working_dir(strain_species_dict, reference_link_path_dict, path):
    #     print('')
    #     strain_specific_dir = dict()
    #     for strain_name, species_code in strain_species_dict.items():
    #         best_ref_path = reference_link_path_dict[strain_name]
    #         print(strain_name, species_code, best_ref_path, path)
    #     return strain_specific_dir
