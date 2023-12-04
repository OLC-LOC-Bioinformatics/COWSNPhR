#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import make_path, run_subprocess, write_to_logfile
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing
from ete3 import Tree
from glob import glob
import xlsxwriter
import shutil
import pandas
import gzip
import math
import xlrd
import os

__author__ = 'adamkoziol'


class TreeMethods(object):

    @staticmethod
    def file_list(file_path):
        """
        Create a list of all gVCF files present in the supplied path. Accepts .gvcf.gz extension only
        :param file_path: type STR: absolute path of folder containing VCF files
        :return vcf_files: sorted list of all VCF files present in the supplied path
        """
        # Use glob to find the acceptable extensions of VCF files in the supplied path
        vcf_files = glob(os.path.join(file_path, '*.gvcf*'))
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
        :return strain_name_dict: Dictionary of strain name: absolute path to VCF file
        """
        # Initialise a dictionary to store the base name of the VCF file: absolute path to file
        strain_name_dict = dict()
        for vcf_file in vcf_files:
            # Split off the extension
            base_name = os.path.splitext(vcf_file)[0]
            base_name = base_name if not base_name.endswith('.gvcf') else base_name.split('.gvcf')[0]
            # Extract the base name from the absolute path plus base name
            strain_name = os.path.basename(base_name)
            strain_name_dict[strain_name] = vcf_file
        return strain_name_dict

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
    def load_gvcf(strain_vcf_dict, threads, qual_cutoff=20):
        """
        Create a multiprocessing pool to parse gVCF files concurrently
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to gVCF file
        :param qual_cutoff: type INT: Quality cutoff value to use. Default is 20
        :param threads: type INT: Number of processes to run concurrently
        :return:
        """
        # Initialise dictionaries to store the parsed gVCF outputs and the closest reference genome
        strain_parsed_vcf_dict = dict()
        strain_best_ref_dict = dict()
        strain_best_ref_set_dict = dict()
        # Create a multiprocessing pool. Limit the number of processes to the number of threads
        p = multiprocessing.Pool(processes=threads)
        # Create a list of all the strain names
        strain_list = [strain_name for strain_name in strain_vcf_dict]
        # Determine the number of strains present in the analyses
        list_length = len(strain_list)
        # Use multiprocessing.Pool.starmap to process the samples in parallel
        # Supply the list of strains, as well as a list the length of the number of strains of each required variable
        for parsed_vcf, strain_best_ref, strain_best_ref_set in p.starmap(TreeMethods.load_gvcf_multiprocessing,
                                                                          zip(strain_list,
                                                                              [strain_vcf_dict] * list_length,
                                                                              [qual_cutoff] * list_length)):
            # Update the dictionaries
            strain_parsed_vcf_dict.update(parsed_vcf)
            strain_best_ref_dict.update(strain_best_ref)
            strain_best_ref_set_dict.update(strain_best_ref_set)
        # Close and join the pool
        p.close()
        p.join()
        return strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

    @staticmethod
    def load_gvcf_multiprocessing(strain_name, strain_vcf_dict, qual_cutoff):
        """
        Load the gVCF files into a dictionary
        :param strain_name: type STR: Name of strain being processed
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to gVCF file
        :param qual_cutoff: type INT: Quality cutoff value to use.
        :return: parsed_vcf_dict: Dictionary of strain name: key: value pairs CHROM': ref_genome, 'REF': ref base,
            'ALT': alt base, 'QUAL': quality score, 'LENGTH': length of feature, 'FILTER': deepvariant filter call,
            'STATS': dictionary of format data
        :return: strain_best_ref_dict: Dictionary of strain name: reference genome parsed from gVCF file. Note that
            this will select only a single 'best reference genome' even if there are multiple contigs in the file
            against which this strain was reference mapped
        :return: strain_best_ref_set_dict: Dictionary of strain name: all reference genomes parsed from gVCF file
        """
        '''
        gVCF header information
        ##fileformat=VCFv4.2
        ##FILTER=<ID=PASS,Description="All filters passed">
        ##FILTER=<ID=RefCall,Description="Genotyping model thinks this site is reference.">
        ##FILTER=<ID=LowQual,Description="Confidence in this variant being real is below calling threshold.">
        ##INFO=<ID=END,Number=1,Type=Integer,Description="End position (for use with symbolic alleles)">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
        ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block.">
        ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
        ##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Variant allele fractions.">
        ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest 
            integer">
        ##contig=<ID=NC_002945.4,length=4349904>
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	<STRAIN NAME>
        '''
        # Initialise dictionaries to store the parsed gVCF outputs and the closest reference genome
        strain_parsed_vcf_dict = dict()
        strain_best_ref_dict = dict()
        strain_best_ref_set_dict = dict()
        vcf_file = strain_vcf_dict[strain_name]
        strain_parsed_vcf_dict[strain_name] = dict()
        # Use gzip to open the compressed file
        with gzip.open(vcf_file, 'r') as gvcf:
            for line in gvcf:
                # Convert the line to string from bytes
                line = line.decode()
                # Skip all the headers
                if line.startswith('#CHROM'):
                    for subline in gvcf:
                        # Convert to string
                        subline = subline.decode()
                        # Split the line on tabs. The components correspond to the #CHROM comment above
                        ref_genome, pos, id_stat, ref, alt_string, qual, filter_stat, info_string, format_stat, \
                            strain = subline.split('\t')
                        # The 'Format' entry consists of several components: GT:GQ:DP:AD:VAF:PL for SNP positions,
                        # and GT:GQ:MIN_DP:PL for all other entries (see quoted information above)
                        # Perform a dictionary comprehension to associate each format component with its
                        # corresponding 'strain' component e.g. FORMAT: GT:GQ:DP:AD:VAF:PL
                        # 'STRAIN' 1/1:54:18:0,18,0:1,0:60,55,0,990,990,990 yields'GT': '1/1', 'GQ': 54, 'DP': '18',
                        # 'AD': 18:0,18,0, 'VAF': 1,0, 'PL': 60,55,0,990,990,990
                        format_dict = {value: strain.split(':')[i].rstrip()
                                       for i, value in enumerate(format_stat.split(':'))}
                        # Initialise the dictionary as required
                        if strain_name not in strain_best_ref_dict:
                            strain_best_ref_dict[strain_name] = ref_genome
                            strain_best_ref_set_dict[strain_name] = {ref_genome}
                        else:
                            strain_best_ref_set_dict[strain_name].add(ref_genome)
                        # The 'info' string will look like this for non-SNP calls: END=1056 (and just a . for SNPs)
                        # Strip off the 'END='
                        if info_string.startswith('END='):
                            info = info_string.split('END=')[1]
                        else:
                            info = pos
                        # Initialise a string to store the sanitised 'ALT" call
                        alt = str()
                        # For SNP calls, the alt_string will look like this: G,<*>, or A,G,<*>, while matches are
                        # simply <*>. Replace the <*> with the reference call, and create a list by splitting on
                        # commas
                        alt_split = alt_string.replace('<*>', ref).split(',')
                        # Initialise the length of the list to 1
                        alt_length = 1
                        # Check if the length of the list is greater than 1 i.e. a SNP call
                        if len(alt_split) > 1:
                            # Iterate through the list
                            for sub_alt in alt_split:
                                # Add each allele to the alt string e.g. initial G,<*> -> G, T -> GT, and A,G,<*> ->
                                # A, G, C -> AGC
                                alt += sub_alt
                                # If there is an insertion, e.g. CGAGACCG,<*>, set alt_length to the length of the
                                # insertion
                                if len(sub_alt) > alt_length:
                                    alt_length = len(sub_alt)
                        # Typecast pos to be an integer
                        pos = int(pos)
                        if ref_genome not in strain_parsed_vcf_dict[strain_name]:
                            strain_parsed_vcf_dict[strain_name][ref_genome] = dict()
                        # SNPs must have a deepvariant filter of 'PASS', be of length one, and have a quality
                        # score above the threshold
                        if filter_stat == 'PASS' and len(ref) == 1 and alt_length == 1 and float(qual) > qual_cutoff:
                            # Populate the dictionary with the required key: value pairs
                            strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                                'CHROM': ref_genome,
                                'REF': ref,
                                'ALT': alt,
                                'QUAL': qual,
                                'LENGTH': 1,
                                'FILTER': filter_stat,
                                'STATS': format_dict
                            }
                        # Insertions must still have a deepvariant filter of 'PASS', but must have a length
                        # greater than one
                        elif filter_stat == 'PASS' and alt_length > 1:
                            strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                                'CHROM': ref_genome,
                                'REF': ref,
                                'ALT': alt,
                                'QUAL': qual,
                                'LENGTH': alt_length,
                                'FILTER': 'INSERTION',
                                'STATS': format_dict
                            }
                        # If the position in the 'info' field does not match pos, and the minimum depth of a gVCF
                        # block is 0, this is considered a deletion
                        elif int(info) != pos and format_dict['MIN_DP'] == 0:
                            # Subtract the starting position (pos) from the final position (info)
                            length = int(info) - pos
                            # Iterate through the range of the deletion, and populate the dictionary for each
                            # position encompassed by this range (add +1 due to needing to include the final
                            # position in the dictionary)
                            for i in range(int(pos), int(info) + 1):
                                strain_parsed_vcf_dict[strain_name][ref_genome][i] = {
                                    'CHROM': ref_genome,
                                    'REF': ref,
                                    'ALT': alt,
                                    'QUAL': qual,
                                    'LENGTH': length,
                                    'FILTER': 'DELETION',
                                    'STATS': format_dict
                                }
        return strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

    @staticmethod
    def load_vcf(strain_vcf_dict, min_depth=10):
        """
        Using vcf.Reader(), load the VCF files. Store the Reader objects, as well as the extracted reference sequence,
        and its associated species code in dictionaries
        :param strain_vcf_dict: type DICT: Dictionary of strain name: list of absolute path to VCF file
        :param min_depth: type INT: Integer of the minimum mapping depth at a site in order for it to be considered
        in the analysis
        :return: strain_vcf_object_dict: Dictionary of strain name: VCF Reader object
        :return: strain_best_ref_dict: Dictionary of strain name: extracted reference genome name
        """
        # Initialise dictionaries to store the parsed gVCF outputs and the closest reference genome
        strain_parsed_vcf_dict = dict()
        strain_best_ref_dict = dict()
        strain_best_ref_set_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            strain_parsed_vcf_dict[strain_name] = dict()
            if vcf_file.endswith('.gz'):
                # Use gzip to open the compressed gVCF file
                filtered = gzip.open(vcf_file, 'r')
            else:
                filtered = open(vcf_file, 'r')
            for line in filtered:
                if vcf_file.endswith('.gz'):
                    # Convert the line to a string from bytes
                    line = line.decode()
                # Add the VCF file header information to the filtered file
                if line.startswith('#'):
                    pass
                else:
                    # Split the line based on the columns
                    ref_genome, pos, id_stat, ref, alt_string, qual, filter_stat, info_string, \
                        format_stat, strain_info = line.rstrip().split('\t')
                    # Initialise the ref_genome key
                    if ref_genome not in strain_parsed_vcf_dict[strain_name]:
                        strain_parsed_vcf_dict[strain_name][ref_genome] = dict()
                    # Initialise the dictionary as required
                    if strain_name not in strain_best_ref_dict:
                        strain_best_ref_dict[strain_name] = ref_genome
                        strain_best_ref_set_dict[strain_name] = {ref_genome}
                    else:
                        strain_best_ref_set_dict[strain_name].add(ref_genome)
                    # # Find the depth entry. e.g. DP=11
                    # depth_group = re.search('(DP=[0-9]+)', info_string)
                    # # Split the depth matching group on '=' and convert the depth to an int
                    # depth = int(str(depth_group.group()).split('=')[1])
                    # Typecast pos to int
                    pos = int(pos)
                    # The 'Format' entry consists of several components: GT:DP:AD:RO:QR:AO:QA:GL for SNP positions,
                    # andGQ:DP:MIN_DP:QR:RO:QA:AO for all other entries (see quoted information above)
                    # Perform a dictionary comprehension to associate each format component with its
                    # corresponding 'strain' component
                    format_dict = dict()
                    for i, category in enumerate(format_stat.split(':')):
                        #
                        value = strain_info.split(':')[i]
                        format_dict[category] = value
                        if category == 'MIN_DP' and value == '0':
                            # The block of deleted sequence will stretch from the current position until the 'END='
                            # position e.g.
                            # Contig_1_138.744  136699  . A  <*> 0 . END=136738 GT:GQ:MIN_DP:PL 0/0:1:0:0,0,0
                            # the deletion is from 136699 - 136738
                            # Split the 'END=136738' string on the '='
                            end_str, end_pos = info_string.split('=')
                            # Typecast the string to an integer
                            end_pos = int(end_pos)
                            # If the end position is not the same as the current position, the deletion is larger than
                            # one bp, iterate through the range of missing positions, and add them all to the dictionary
                            if int(end_pos) > pos:
                                for del_pos in range(pos, end_pos + 1):
                                    # Store the zero coverage entries
                                    strain_parsed_vcf_dict[strain_name][ref_genome][del_pos] = {
                                        'CHROM': ref_genome,
                                        'REF': ref,
                                        'ALT': alt_string,
                                        'QUAL': qual,
                                        'LENGTH': 1,
                                        'FILTER': 'DELETION',
                                        'STATS': format_dict
                                    }
                            else:
                                # Store the zero coverage entries for the one bp deletions
                                strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                                    'CHROM': ref_genome,
                                    'REF': ref,
                                    'ALT': alt_string,
                                    'QUAL': qual,
                                    'LENGTH': 1,
                                    'FILTER': 'DELETION',
                                    'STATS': format_dict
                                }
                    try:
                        depth = float(format_dict['DP'].split(',')[0])
                    except KeyError:
                        depth = 0
                    if len(ref) == 1:
                        # Populate the dictionary with the appropriate filter information
                        if filter_stat == 'PASS':
                            split_alt = alt_string.split(',')[0]
                            if depth >= min_depth:
                                # Populate the dictionaries with the number of high quality SNPs
                                strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                                    'CHROM': ref_genome,
                                    'REF': ref,
                                    'ALT': alt_string,
                                    'QUAL': qual,
                                    'LENGTH': len(split_alt),
                                    'FILTER': 'PASS',
                                    'STATS': format_dict
                                }
                        else:
                            if pos not in strain_parsed_vcf_dict[strain_name][ref_genome]:
                                # Populate the dictionaries with the regions that match
                                strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                                    'CHROM': ref_genome,
                                    'REF': ref,
                                    'ALT': alt_string,
                                    'QUAL': qual,
                                    'LENGTH': 1,
                                    'FILTER': 'MATCH',
                                    'STATS': format_dict
                                }
                    # Store all indels
                    else:
                        strain_parsed_vcf_dict[strain_name][ref_genome][pos] = {
                            'CHROM': ref_genome,
                            'REF': ref,
                            'ALT': alt_string,
                            'QUAL': qual,
                            'LENGTH': len(alt_string),
                            'FILTER': 'INSERTION',
                            'STATS': format_dict
                        }
            filtered.close()
        return strain_parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

    @staticmethod
    def summarise_gvcf_outputs(strain_parsed_vcf_dict):
        """
        Count the number of locations that PASS filter (SNP call), are considered INSERTIONS, or DELETIONS
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: reference position: parsed VCF dictionary
        for that reference position
        :return: pass_dict: Dictionary of strain name: number of locations with 'PASS' filter
        :return: insertion_dict: Dictionary of strain name: number of locations with 'INSERTION' filter
        :return: deletion_dict: Dictionary of strain name: number of locations with 'DELETION' filter
        """
        # Initialise dictionaries
        pass_dict = dict()
        insertion_dict = dict()
        deletion_dict = dict()
        for strain_name, ref_dict in strain_parsed_vcf_dict.items():
            # Add the strain name key to the dictionaries
            pass_dict[strain_name] = 0
            insertion_dict[strain_name] = 0
            deletion_dict[strain_name] = 0
            for ref_chrom, pos_dict in ref_dict.items():
                for pos, filter_dict in pos_dict.items():
                    # Extract the value from the 'FILTER' key from the dictionary
                    filter_value = filter_dict['FILTER']
                    # Increment the appropriate dictionary if the filter matches
                    if filter_value == 'PASS':
                        pass_dict[strain_name] += 1
                    # As the dictionary is based on the reference position, insertions will be considered a single base
                    # Use the 'LENGTH filter to add the total insertion length to the dictionary
                    elif filter_value == 'INSERTION':
                        insertion_dict[strain_name] += filter_dict['LENGTH']
                    elif filter_value == 'DELETION':
                        deletion_dict[strain_name] += 1
        return pass_dict, insertion_dict, deletion_dict

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
        :return reference_link_dict: Dictionary of reference file name: absolute path containing reference file
        :return reference_strain_dict: Dictionary of strain name: absolute path to best reference file
        """
        # Initialise dictionaries
        reference_link_dict = dict()
        reference_link_path_dict = dict()
        reference_strain_dict = dict()
        # Read in the .csv file with reference file name: relative symbolic link information
        with open(os.path.join(dependency_path, 'reference_links.csv'), 'r') as reference_paths:
            for line in reference_paths:
                # Extract the link information
                reference, linked_file = line.rstrip().split(',')
                reference_link_dict[reference] = os.path.join(dependency_path, linked_file)
        # Use the strain-specific best reference genome name to extract the relative symlink information
        for strain_name, best_ref in strain_best_ref_fasta_dict.items():
            reference_link_path_dict[strain_name] = os.path.join(os.path.dirname(reference_link_dict[best_ref]))
            reference_strain_dict[strain_name] = os.path.join(dependency_path, reference_link_dict[best_ref])
        return reference_link_path_dict, reference_link_dict, reference_strain_dict

    @staticmethod
    def consolidate_group_ref_genomes(reference_link_dict, strain_best_ref_dict):
        """
        Brucella is mapped against a FASTA file with multiple contigs (NC_017251-NC_017250.fasta), after parsing the
        gVCF file, the best_ref will be one of: NC_017250.1 or NC_017251.1. Link the best_ref to the combined file
        :param reference_link_dict: type DICT: Dictionary of strain name:
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted best reference genome from
        gVCF file
        :return: strain_consolidated_ref_dict
        """
        # Initialise a dictionary to store the consolidated best_ref name
        strain_consolidated_ref_dict = dict()
        for strain_name, best_ref in strain_best_ref_dict.items():
            for reference in reference_link_dict:
                # Find the best_ref extracted from the gVCF file e.g. NC_017251.1 in the name of all the reference files
                # e.g. NC_017251-NC_017250.fasta
                if best_ref.split('.')[0] in reference:
                    consolidated_ref = reference.split('.')[0]
                    # Set the name of the combined reference name e.g. NC_017251-NC_017250
                    strain_consolidated_ref_dict[strain_name] = consolidated_ref
        return strain_consolidated_ref_dict

    @staticmethod
    def extract_defining_snps(reference_strain_dict, strain_species_dict):
        """
        Load the Excel files containing species-specific groupings of defining SNPs
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :return: defining_snp_dict: Dictionary of species code: dictionary of grouping: reference genome: defining SNP
        """
        # Initialise a dictionary to store the species-specific groups of defining SNPs
        defining_snp_dict = dict()
        for strain_name, best_ref_path in reference_strain_dict.items():
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            # Set the name of the Excel file containing the groups and their defining SNPs
            defining_snp_xlsx = os.path.join(best_ref_path, 'DefiningSNPsGroupDesignations.xlsx')
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
    def load_gvcf_snp_positions(strain_parsed_vcf_dict, strain_consolidated_ref_dict):
        """
        Parse the gVCF files, and extract all the query and reference genome-specific SNP locations as well as the
        reference sequence
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param strain_consolidated_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :return: consolidated_ref_snp_positions: Dictionary of reference name: reference chromosome: pos: ref sequence
        :return: ref_snp_positions: Dictionary of reference chromosome name: absolute position: reference base call
        :return: strain_snp_positions: Dictionary of strain name: reference chromosome: list of strain-specific
        SNP positions
        """
        # Initialise dictionaries to store the SNP positions
        consolidated_ref_snp_positions = dict()
        strain_snp_positions = dict()
        ref_snp_positions = dict()
        for strain_name, ref_dict in strain_parsed_vcf_dict.items():
            best_ref = strain_consolidated_ref_dict[strain_name]
            strain_snp_positions[strain_name] = dict()
            # Initialise the reference genome key as required
            if best_ref not in consolidated_ref_snp_positions:
                consolidated_ref_snp_positions[best_ref] = dict()
            # Iterate through all the positions
            for ref_chrom, vcf_dict in ref_dict.items():
                for pos, pos_dict in vcf_dict.items():
                    chrom = pos_dict['CHROM']
                    if chrom not in ref_snp_positions:
                        ref_snp_positions[chrom] = dict()
                    if chrom not in strain_snp_positions[strain_name]:
                        strain_snp_positions[strain_name][chrom] = list()
                    if chrom not in consolidated_ref_snp_positions[best_ref]:
                        consolidated_ref_snp_positions[best_ref][chrom] = dict()
                    consolidated_ref_snp_positions[best_ref][chrom][pos] = pos_dict['REF']
                    # Only consider locations that are called 'PASS' in the dictionary
                    if pos_dict['FILTER'] == 'PASS':
                        # Populate the dictionary with the position and the reference sequence at that position
                        ref_snp_positions[chrom][pos] = pos_dict['REF']
                        strain_snp_positions[strain_name][chrom].append(pos)
        return consolidated_ref_snp_positions, strain_snp_positions, ref_snp_positions

    @staticmethod
    def group_strains(strain_snp_positions):
        """
        Find all the SNP positions present in the strains in each group
        :param strain_snp_positions: type DICT: Dictionary of strain name: ref chromosome: positions
        :return: group_positions_set: Dictionary of species: group: reference chromosome name: positions
        return: strain_groups: Dictionary of strain name: ['group']
        return strain_species_dict: Dictionary of strain name: 'species'
        """
        # Initialise a dictionary to store group-specific SNP positions for each reference chromosome
        group_positions_set = dict()
        strain_groups = dict()
        strain_species_dict = dict()
        # In order to re-use code, the species and group must be provided. Use generic values
        species = 'species'
        group = 'group'
        for strain_name in strain_snp_positions:
            strain_groups[strain_name] = [group]
            strain_species_dict[strain_name] = species
            for ref_chrom, pos_list in sorted(strain_snp_positions[strain_name].items()):
                # Initialise the species key in the dictionary if required
                if species not in group_positions_set:
                    group_positions_set[species] = dict()
                # Initialise the group key in the dictionary if required
                if group not in group_positions_set[species]:
                    group_positions_set[species][group] = dict()
                if ref_chrom not in group_positions_set[species][group]:
                    group_positions_set[species][group][ref_chrom] = set()
                # Add the group-specific positions to the set
                for pos in pos_list:
                    group_positions_set[species][group][ref_chrom].add(pos)
        return group_positions_set, strain_groups, strain_species_dict

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
        for strain_name, snp_dict in strain_snp_positions.items():
            strain_groups[strain_name] = ['All']
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
                        for ref_chrom, snp_positions in snp_dict.items():
                            # If the position is in the list of strain-specific SNP positions, add the group to the list
                            if position in snp_positions:
                                strain_groups[strain_name].append(group)
        return strain_groups

    @staticmethod
    def determine_group_snp_positions(strain_snp_positions, strain_groups, strain_species_dict):
        """
        Find all the group-specific SNP positions
        :param strain_snp_positions: type DICT: Dictionary of strain name: reference chromosome: list of strain-specific
        SNP positions
        :param strain_groups: type DICT:  Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :return: group_positions_set: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        """
        # Initialise a dictionary to store all of the species: group-specific SNP positions
        group_positions_set = dict()
        for strain_name, groups in strain_groups.items():
            for ref_chrom, pos_list in strain_snp_positions[strain_name].items():
                # Extract the species code from the dictionary
                species = strain_species_dict[strain_name]
                # Initialise the species key in the dictionary if required
                if species not in group_positions_set:
                    group_positions_set[species] = dict()
                for group in groups:
                    # Initialise the group key in the dictionary if required
                    if group not in group_positions_set[species]:
                        group_positions_set[species][group] = dict()
                    if ref_chrom not in group_positions_set[species][group]:
                        group_positions_set[species][group][ref_chrom] = set()
                    # Add the group-specific positions to the set
                    for pos in pos_list:
                        group_positions_set[species][group][ref_chrom].add(pos)
        return group_positions_set

    @staticmethod
    def density_filter_snps(group_positions_set, window_size=1000, threshold=2):
        """
        Remove any SNPs from regions of a defined window size with two or more SNPs
        :param group_positions_set: type DICT: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        :param window_size: type INT: Window size to use when filtering SNPs. Default is 1000
        :param threshold: type INT: Number of SNPs present within the window to trigger filtering
        :return: filtered_group_positions: Dictionary of species code: group name: reference chromosome: set of SNP
        unfiltered SNP positions
        """
        # Divide the window size by two to yield the range to use when filtering putative recombinant SNPs
        bp_range = int(window_size / 2)
        # Initialise a dictionary to store the SNPs that pass filter
        filtered_group_positions = dict()
        for species, group_dict in group_positions_set.items():
            # Initialise the species key
            filtered_group_positions[species] = dict()
            for group, ref_dict in group_dict.items():
                # Initialise the group key
                filtered_group_positions[species][group] = dict()
                for ref_chrom, pos_list in ref_dict.items():
                    if ref_chrom not in filtered_group_positions[species][group]:
                        filtered_group_positions[species][group][ref_chrom] = set()
                    for pos in pos_list:
                        # Count of SNPs in range of this SNP
                        add_pos = 0
                        # Simple window of 1001 bp (pos - 500 to pos + 500)
                        for i in range(pos - bp_range, pos + bp_range):
                            # If there is another SNP in this range, increment add_pos
                            if i in pos_list and i != pos:
                                add_pos += 1
                        # If the position passes the filter, add it to the dictionary
                        if add_pos <= threshold:
                            filtered_group_positions[species][group][ref_chrom].add(pos)
        return filtered_group_positions

    @staticmethod
    def mask_ref_genome(reference_strain_dict, logfile):
        """
        Use nucmer to mask repetitive and low complexity regions in the reference genome.
        Uses logic from: https://bioinformatics.stackexchange.com/a/2374
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        dependency folder type STR: Name and absolute path of reference genome
        :param logfile: type STR: Absolute path to logfile
        :return: coords_dict: Dictionary of reference genome name: absolute path to nucmer-created coords file
        """
        # Initialise the dictionary
        coords_dict = dict()
        for strain_name, ref_genome in reference_strain_dict.items():
            # Remove the file extension from the path and filename of the reference file
            ref_name = os.path.splitext(ref_genome)[0]
            # Create the nucmer masking command
            # --maxmatch: use all anchor matches regardless of their uniqueness
            # --nosimplify: simplify alignments by removing shadowed clusters.
            #   Turn this option off if aligning a sequence to itself to look for repeats (default --simplify)
            # --coords: automatically generate the original NUCmer1.1 coords output file using the 'show-coords' program
            # -p output file prefix
            nucmer_cmd = 'nucmer --maxmatch --nosimplify --coords -p {ref_name}_mummer_maskfile {ref_file} {ref_file}'\
                .format(ref_name=ref_name,
                        ref_file=ref_genome)
            # The desired output file is the .coords file
            coordsfile = ref_genome.replace('.fasta', '_mummer_maskfile.coords')
            # If the delta file does not exist, run the nucmer command
            if not os.path.isfile(coordsfile):
                out, err = run_subprocess(command=nucmer_cmd)
                # Write the stdout and stderr to the logfile
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=nucmer_cmd,
                                                           out=out),
                                 err=err,
                                 logfile=logfile)
            # Add the path to the coords files
            coords_dict[ref_name] = coordsfile
        # Return the dictionary
        return coords_dict

    @staticmethod
    def determine_coordinates(strain_groups, coords_dict, cutoff=90):
        """
        Parse the coords file generated by nucmer, and store the coordinates of the repeat and low-complexity regions
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param coords_dict: type DICT: Dictionary of name of reference genome: name and absolute path to the nucmer-
        generated coords file
        :param cutoff: type INT: Percentage cutoff to use for identity values
        :return: mask_pos_dict: Dictionary of species: group: reference chromosome: set of masked positions
        """
        # Initialise the dictionary to store the masked positions
        mask_pos_dict = dict()
        for species, group_list in strain_groups.items():
            if species not in mask_pos_dict:
                mask_pos_dict[species] = dict()
            for group in group_list:
                if group not in mask_pos_dict[species]:
                    mask_pos_dict[species][group] = dict()
                # Unpack the reference name: name and path of coords file dictionary
                for ref_name, coordsfile in coords_dict.items():
                    # Set the name and path of the output .bed file by replacing the .coords file extension with .bed
                    bedfile = coordsfile.replace('.coords', '.bed')
                    with open(bedfile, 'w') as bed:
                        with open(coordsfile, 'r') as coords:
                            for line in coords:
                                # Ignore the header information
                                if line.startswith('='):
                                    for data in coords:
                                        # Split the data line on whitespace - filter the list for non-empty values that
                                        # are not '|'
                                        # [S1]  [E1]| [S2] [E2]| [LEN 1] [LEN 2] |  [% IDY]  | [TAGS]
                                        # =============================================================================
                                        #   1   340|  1   340  |   340      340  |   100.00  | Contig_10_403 Contig_10_
                                        s1, e1, s2, e2, len1, len2, pid, contig1, contig2 \
                                            = [value for value in data.rstrip().split() if value and value != '|']
                                        # Ensure that the match isn't to itself, and that the percent identity of the
                                        # match is above the desired cutoff value
                                        if s1 != s2 and e1 != e2 and float(pid) > cutoff:
                                            # Write the name of the chromosome, as well as the start and end positions
                                            # of the repeat regions to the .bed file
                                            bed.write('{chrom}\t{chrom_start}\t{chrom_end}\n'
                                                      .format(chrom=contig1,
                                                              chrom_start=s1,
                                                              chrom_end=e1))
                                            if contig1 not in mask_pos_dict[species][group]:
                                                mask_pos_dict[species][group][contig1] = set()
                                            # Add the range of the masked region to the set
                                            for i in range(int(s1), int(e1) + 1):
                                                mask_pos_dict[species][group][contig1].add(i)
        return mask_pos_dict

    @staticmethod
    def load_supplied_mask(strain_groups, maskfile):
        """
        Load the coordinates from the user-supplied mask file
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param maskfile: type STR: Name and absolute path to the user-supplied maskfile
        :return: supplied_mask_pos_dict: Dictionary of species: group: reference chromosome: set of masked positions
        """
        # Initialise the dictionary
        supplied_mask_pos_dict = dict()
        for species, group_list in strain_groups.items():
            if species not in supplied_mask_pos_dict:
                supplied_mask_pos_dict['species'] = dict()
            for group in group_list:
                if group not in supplied_mask_pos_dict['species']:
                    supplied_mask_pos_dict['species'][group] = dict()
                # Open the mask file
                with open(maskfile, 'r') as mask:
                    for line in mask:
                        # Split the line on tabs
                        # e.g. Contig_13_99.0332	49320	49360
                        chrom, chrom_start, chrom_end = line.rstrip().split('\t')
                        # Initialise the reference chromosome in the dictionary if necessary
                        if chrom not in supplied_mask_pos_dict['species'][group]:
                            supplied_mask_pos_dict['species'][group][chrom] = set()
                        # Add the range of the masked region to the set
                        for i in range(int(chrom_start), int(chrom_end) + 1):
                            supplied_mask_pos_dict['species'][group][chrom].add(i)
        return supplied_mask_pos_dict

    @staticmethod
    def filter_masked_snp_positions(group_positions_set, filtered_group_positions, mask_pos_dict,
                                    supplied_mask_pos_dict):
        """
        Determine which SNP positions are present in the masked area of the reference genome, and remove them from the
        analysis
        :param group_positions_set: type DICT: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        :param filtered_group_positions: type DICT: Dictionary of species: group: reference chromosome: set of
        density-filtered SNP positions
        :param mask_pos_dict: type DICT: Dictionary of species: group: reference chromosome: set of
        nucmer-calculated masked positions
        :param supplied_mask_pos_dict: type DICT: Dictionary of species: group: reference chromosome: set of user-
        supplied masked positions
        :return: filtered_masked_group_positions: Dictionary of species: group: reference chromosome: set of SNP
        positions that pass filter
        :return: filter_reasons: Dictionary of species: group: reference chromosome: position: list of reasons
        position was excluded
        """
        filtered_masked_group_positions = dict()
        filter_reasons = dict()
        for species, group_dict in group_positions_set.items():
            # Initialise the species key
            filtered_masked_group_positions[species] = dict()
            filter_reasons[species] = dict()
            for group, ref_dict in group_dict.items():
                # Initialise the group key
                filtered_masked_group_positions[species][group] = dict()
                filter_reasons[species][group] = dict()
                # Iterate through all the reference chromosomes in the dictionary
                for ref_chrom, pos_list in ref_dict.items():
                    # Initialise the reference chromosome in the dictionaries
                    if ref_chrom not in filtered_masked_group_positions[species][group]:
                        filtered_masked_group_positions[species][group][ref_chrom] = set()
                    if ref_chrom not in filter_reasons[species][group]:
                        filter_reasons[species][group][ref_chrom] = dict()
                    for pos in sorted(pos_list):
                        # Initialise the list of reasons that a position is filtered
                        reasons = list()
                        # Density filtering check
                        if pos not in filtered_group_positions[species][group][ref_chrom]:
                            reasons.append('density')
                        # Nucmer-calculated mask
                        try:
                            calculated_mask = mask_pos_dict[species][group][ref_chrom]
                            # Check to see if this position is present in the set of calculated masked positions
                            if pos in calculated_mask:
                                reasons.append('calculated mask')
                        except KeyError:
                            pass
                        # User-supplied mask
                        try:
                            user_mask = supplied_mask_pos_dict[species][group][ref_chrom]
                            # Check to see if this position is present in the set of user-supplied masked positions
                            if pos in user_mask:
                                reasons.append('user mask')
                        except KeyError:
                            pass
                        # If the position was filtered, add the list of reasons to the dictionary
                        if reasons:
                            filter_reasons[species][group][ref_chrom][pos] = reasons
                        # If the position didn't get filtered, add it to the set of unfiltered positions
                        else:
                            filtered_masked_group_positions[species][group][ref_chrom].add(pos)
        return filtered_masked_group_positions, filter_reasons

    @staticmethod
    def load_snp_sequence(strain_parsed_vcf_dict, strain_consolidated_ref_dict, group_positions_set, strain_groups,
                          strain_species_dict, consolidated_ref_snp_positions, iupac):
        """
        Parse the gVCF-derived dictionaries to determine the strain-specific sequence at every SNP position for every
        group
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param strain_consolidated_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param group_positions_set: type DICT: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :param consolidated_ref_snp_positions: type DICT: Dictionary of reference name: absolute position: reference
        base call
        :param iupac: type DICT: Dictionary of degenerate code: nucleotides included in group
        :return: group_strain_snp_sequence: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :return: species_group_best_ref: Dictionary of species code: group name; best ref
        """
        # Initialise a dictionary to store all the SNP locations, and the reference sequence for each reference
        # genome
        group_strain_snp_sequence = dict()
        species_group_best_ref = dict()
        for strain_name, ref_dict in strain_parsed_vcf_dict.items():
            # Extract the set of groups to which the strain belongs
            groups = strain_groups[strain_name]
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            best_ref = strain_consolidated_ref_dict[strain_name]
            # Initialise the dictionary with the species key if required
            if species not in group_strain_snp_sequence:
                group_strain_snp_sequence[species] = dict()
                species_group_best_ref[species] = dict()
            for group in groups:
                # Add the group and strain name keys to the dictionary
                if group not in group_strain_snp_sequence[species]:
                    group_strain_snp_sequence[species][group] = dict()
                if strain_name not in group_strain_snp_sequence[species][group]:
                    group_strain_snp_sequence[species][group][strain_name] = dict()
                if group not in species_group_best_ref[species]:
                    species_group_best_ref[species][group] = best_ref
                for ref_chrom, position_set in group_positions_set[species][group].items():
                    # The reference sequence only needs to be determined once. This variable will check to see if it
                    # already has been added to the dictionary
                    write_ref = False
                    if best_ref not in group_strain_snp_sequence[species][group]:
                        group_strain_snp_sequence[species][group][best_ref] = dict()
                    if ref_chrom not in group_strain_snp_sequence[species][group][best_ref]:
                        group_strain_snp_sequence[species][group][best_ref][ref_chrom] = dict()
                        write_ref = True
                    for pos in sorted(position_set):
                        ref_pos = consolidated_ref_snp_positions[best_ref][ref_chrom][pos]
                        # Include the reference position if necessary
                        if write_ref:
                            group_strain_snp_sequence[species][group][best_ref][ref_chrom][pos] = ref_pos
                        # gVCF blocks will compress stretches of normal matches. I haven't added these regions to the
                        # dictionary, so these positions will yield KeyErrors
                        try:
                            # Extract the gVCF dictionary from the position-specific dictionary
                            pos_dict = ref_dict[ref_chrom][pos]
                            if ref_chrom not in group_strain_snp_sequence[species][group][strain_name]:
                                group_strain_snp_sequence[species][group][strain_name][ref_chrom] = dict()
                            # Deletions are recorded as a '-'
                            if pos_dict['FILTER'] == 'DELETION':
                                group_strain_snp_sequence[species][group][strain_name][ref_chrom][pos] = '-'
                            elif pos_dict['FILTER'] == 'INSERTION':
                                group_strain_snp_sequence[species][group][strain_name][ref_chrom][
                                    pos] = pos_dict['REF']
                            else:
                                if pos_dict['LENGTH'] == 1:
                                    # Determine if the allele frequency is a mixed population
                                    try:
                                        # Create a list consisting of the reference position for now
                                        mixed_components = [ref_pos]
                                        # Split the 'ALT' value on commas e.g. 'T', '<*>'
                                        other_components = list()
                                        for component in pos_dict['ALT'].split(','):
                                            # Ignore the '<*>'
                                            if component != '<*>':
                                                # If the component is a single base, add it to the list - otherwise,
                                                # it represents an insertion
                                                if len(component) == 1:
                                                    other_components.append(component)
                                        if other_components:
                                            depth = float(pos_dict['STATS']['DP'])
                                            # The allele depth is composed of three values e.g. 0,17,0
                                            # The first is the number of reference-matching bases, the second is the
                                            # number of bases matching the alternate call, and the third is 'other'?
                                            allele_depth = float(pos_dict['STATS']['AD'].split(',')[1])
                                            allele_freq = allele_depth / depth
                                            # If the alternate allele constitutes less than 75%, but over 50% of the
                                            # total depth, determine what the degenerate base call is
                                            if allele_freq < 0.75:
                                                if allele_freq >= 0.50:
                                                    for component in other_components:
                                                        mixed_components.append(component)
                                                    # Determine the IUPAC code for any multi-allelic sites
                                                    for code, components in iupac.items():
                                                        if sorted(mixed_components) == sorted(components):
                                                            group_strain_snp_sequence[species][group][strain_name][
                                                                ref_chrom][pos] = code
                                                # Use the 'REF' allele
                                                else:
                                                    group_strain_snp_sequence[species][group][strain_name][ref_chrom][
                                                        pos] = pos_dict['REF']
                                            # Otherwise use the alt allele
                                            else:
                                                group_strain_snp_sequence[species][group][strain_name][ref_chrom][
                                                    pos] = pos_dict['ALT'][0]
                                    # The FreeBayes stats dictionary doesn't have the VAF. Assign the alternate allele
                                    # as in the original vSNP.
                                    except KeyError:
                                        # IF "AC" (alternate called alleles) is 1, find the IUPAC code of the ref + alt
                                        # allele combination e.g. 13-1950 pos 714775: ref: G, alt: A, call: R
                                        if pos_dict['STATS']['AC'] == '1':
                                            for code, components in iupac.items():
                                                if sorted([pos_dict['REF'], pos_dict['ALT']]) == sorted(components):
                                                    group_strain_snp_sequence[species][group][strain_name][
                                                        ref_chrom][pos] = code
                                        # Otherwise use the alt allele
                                        else:
                                            group_strain_snp_sequence[species][group][strain_name][ref_chrom][pos] = \
                                                pos_dict['ALT'][0]
                        # If the entry isn't in the dictionary, it is because it matches the reference sequence or
                        # because there is a deletion
                        except KeyError:
                            # Create a copy of the dictionary with the strain-specific reference chromosome information
                            # derived from the gVCF file
                            data = strain_parsed_vcf_dict[strain_name][ref_chrom]
                            # Find the closest position below the current position in the dictionary
                            closest_pos = min(data.keys(), key=lambda k: abs(k - pos))
                            # If the position is a DELETION, store a -
                            if data[closest_pos]['FILTER'] == 'DELETION':
                                group_strain_snp_sequence[species][group][strain_name][ref_chrom][pos] = '-'
                            # Otherwise, the position should match the reference genome sequence
                            else:
                                try:
                                    if ref_chrom not in group_strain_snp_sequence[species][group][strain_name]:
                                        group_strain_snp_sequence[species][group][strain_name][ref_chrom] = dict()
                                    group_strain_snp_sequence[species][group][strain_name][ref_chrom][pos] \
                                        = consolidated_ref_snp_positions[best_ref][ref_chrom][pos]
                                except KeyError:
                                    pass
        return group_strain_snp_sequence, species_group_best_ref

    @staticmethod
    def find_identical_calls(group_strain_snp_sequence):
        """
        Remove any positions that have all identical SNP calls
        :param group_strain_snp_sequence: type DICT: Dictionary of species: group: strain name: reference chromosome:
        position: sequence
        :return: ident_group_positions: Dictionary of species: group: reference chromosome: set of identical positions
        """
        # Initialise dictionary to store the number of strains with a particular base for positions of interest
        snp_count_dict = dict()
        # Dictionary to store the number of strains present in each grouping
        species_group_strain_dict = dict()
        # Dictionary to store the identical SNP positions
        ident_group_positions = dict()
        for species, group_dict in group_strain_snp_sequence.items():
            snp_count_dict[species] = dict()
            species_group_strain_dict[species] = dict()
            for group, strain_dict in group_dict.items():
                snp_count_dict[species][group] = dict()
                # Count the number of strains in the dictionary - if the sequence is the same for all of them, it
                # will be excluded
                species_group_strain_dict[species][group] = len(strain_dict)
                for strain_name, ref_dict in strain_dict.items():
                    # Iterate through all the reference chromosomes
                    for ref_chrom, pos_dict in ref_dict.items():
                        if ref_chrom not in snp_count_dict[species][group]:
                            snp_count_dict[species][group][ref_chrom] = dict()
                        for pos, seq in pos_dict.items():
                            if pos not in snp_count_dict[species][group][ref_chrom]:
                                snp_count_dict[species][group][ref_chrom][pos] = dict()
                            # If the base at this position hasn't been encountered, set the count to 1
                            if seq not in snp_count_dict[species][group][ref_chrom][pos]:
                                snp_count_dict[species][group][ref_chrom][pos][seq] = 1
                            # Otherwise, increment the count of this sequence
                            else:
                                snp_count_dict[species][group][ref_chrom][pos][seq] += 1
        # Find all the positions that are not all the same
        for species, group_dict in snp_count_dict.items():
            # Initialise the species key
            ident_group_positions[species] = dict()
            for group, ref_dict in group_dict.items():
                # Initialise the group key
                ident_group_positions[species][group] = dict()
                for ref_chrom, pos_dict in ref_dict.items():
                    # Initialise the reference chromosome key as required
                    if ref_chrom not in ident_group_positions[species][group]:
                        ident_group_positions[species][group][ref_chrom] = set()
                    for pos, seq_dict in pos_dict.items():
                        for current_seq, count in seq_dict.items():
                            # If the count is equal than the number of strains, add the position to the set
                            # of identical positions
                            if count >= species_group_strain_dict[species][group]:
                                ident_group_positions[species][group][ref_chrom].add(pos)
        return ident_group_positions

    @staticmethod
    def create_multifasta(group_strain_snp_sequence, fasta_path, group_positions_set, strain_parsed_vcf_dict,
                          species_group_best_ref, reference_strain_dict, ident_group_positions, nested=True):
        """
        Create a multiple sequence alignment in FASTA format for each group from all the SNP positions for the group
        :param group_strain_snp_sequence: type DICT: Dictionary of species: group: strain name: reference chromosome:
        position: sequence
        :param fasta_path: type STR: Absolute path of folder in which alignments are to be created
        :param group_positions_set: type DICT: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param ident_group_positions: type DICT: Dictionary of species: group: reference chromosome: set of identical
        positions
        :param nested: type BOOL: Boolean on whether the multi-FASTA files should be created in the normal directory
        structure, or within the fasta_path
        :return: group_fasta_dict: Dictionary of species code: group name: FASTA file created for the group
        :return: group_folders: Set of absolute paths to folders for each group
        :return: species_folders: Set of absolute path to folders for each species
        """
        # Initialise variables to return
        group_fasta_dict = dict()
        group_folders = set()
        species_folders = set()
        # Clear out the fasta_path to ensure that no previously processed FASTA files are present, as new outputs will
        # be appended to the old outputs
        try:
            shutil.rmtree(fasta_path)
        except FileNotFoundError:
            pass
        for species, group_dict in group_strain_snp_sequence.items():
            # Initialise the species key
            group_fasta_dict[species] = dict()
            # Add the absolute path of the species-specific folder to the set of all species folders
            species_folders.add(os.path.join(fasta_path, species))
            for group, strain_dict in group_dict.items():
                # Set the output_dir as appropriate based on whether nesting is requested
                if nested:
                    output_dir = os.path.join(fasta_path, species, group)
                else:
                    output_dir = fasta_path
                make_path(output_dir)
                # Add the group-specific folder to the set of all group folders
                group_folders.add(output_dir)
                best_ref = species_group_best_ref[species][group]
                # Use SeqIO to parse all the records in the reference FASTA file
                ref_records = SeqIO.to_dict(SeqIO.parse(reference_strain_dict[best_ref], 'fasta'))
                for strain_name, chrom_dict in strain_dict.items():
                    # Create a string to store the strain-specific sequence
                    strain_group_seq = str()
                    # Iterate through all the chromosomes in the reference genome
                    for ref_chrom, position_set in group_positions_set[species][group].items():
                        for pos in sorted(position_set):
                            # Don't add the sequence if it is identical for all strains
                            if pos not in ident_group_positions[species][group][ref_chrom]:
                                # Adjust sequence
                                ref_seq = str(ref_records[ref_chrom].seq)[pos - 1]
                                try:
                                    sequence = chrom_dict[ref_chrom][pos]
                                    if len(sequence) > 1:
                                        strain_group_seq += ref_seq
                                    else:
                                        # Append each base to the string
                                        strain_group_seq += sequence
                                except KeyError:
                                    try:
                                        sequence_dict = strain_parsed_vcf_dict[strain_name][ref_chrom][pos]
                                        if sequence_dict['FILTER'] == 'DELETION':
                                            strain_group_seq += '-'
                                        else:
                                            strain_group_seq += ref_seq
                                    except KeyError:
                                        strain_group_seq += ref_seq
                    # Create a SeqRecord from the sequence string. Use the strain name as the id
                    record = SeqRecord(Seq(strain_group_seq),
                                       id=strain_name,
                                       description='')
                    # Set the name of the FASTA alignment file
                    group_fasta = os.path.join(output_dir, 'alignment.fasta')
                    # Use SeqIO to append the FASTA sequence to the alignment file
                    with open(group_fasta, 'a+') as fasta:
                        SeqIO.write(record, fasta, 'fasta')
                    # Add the alignment file to the set of all alignment files
                    group_fasta_dict[species][group] = group_fasta
        return group_folders, species_folders, group_fasta_dict

    @staticmethod
    def snp_summary(group_strain_snp_sequence, species_group_best_ref, reference_strain_dict, group_positions_set,
                    filter_reasons, strain_parsed_vcf_dict, filtered_group_positions, mask_pos_dict,
                    supplied_mask_pos_dict, ident_group_positions, summary_path):
        """
        Create two summary tables. snv_summary.tsv has details on every SNV position extracted from the global VCF
        files. contig_summary_tsv has details for every reference contig
        :param group_strain_snp_sequence: type DICT: Dictionary of species: group: strain name: reference chromosome:
        position: sequence
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param group_positions_set: type DICT: Dictionary of species code: group name: reference chromosome: set of
        group-specific SNP positions
        :param filter_reasons: type DICT: Dictionary of species: group: reference chromosome: position: list of reasons
        position was excluded
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param filtered_group_positions: type DICT: Dictionary of species: group: reference chromosome: set of
        density-filtered SNP positions
        :param mask_pos_dict: type DICT: Dictionary of species: group: reference chromosome: set of
        nucmer-calculated masked positions
        :param supplied_mask_pos_dict: type DICT: Dictionary of species: group: reference chromosome: set of user-
        supplied masked positions
        :param ident_group_positions: type DICT: Dictionary of species: group: reference chromosome: set of identical
        positions
        :param summary_path: type STR: Absolute path to folder in which summary reports are to be created
        """
        # Create the summary path as required
        make_path(summary_path)
        with open(os.path.join(summary_path, 'snv_summary.tsv'), 'w') as snp_summary:
            with open(os.path.join(summary_path, 'contig_summary.tsv'), 'w') as pos_summary:
                for species, group_dict in group_strain_snp_sequence.items():
                    for group, strain_dict in group_dict.items():
                        # Extract the name of the reference genome from the species_group_best_ref dictionary using
                        # the species code and the group name
                        best_ref = species_group_best_ref[species][group]
                        # Use SeqIO to parse all the records in the reference FASTA file
                        ref_records = SeqIO.to_dict(SeqIO.parse(reference_strain_dict[best_ref], 'fasta'))
                        # Create the header strings
                        snp_summary_header = 'Contig\tPos\tStatus\tReason\t'
                        pos_summary_header = 'Contig\tTotalLength\tTotalInvalid\tTotalValid\tTotalValidInCore\t' \
                                             'PercentValidInCore\tPercentTotalValidInCore\n'
                        # Initialise variables
                        snp_summary_body = str()
                        pos_summary_body = str()
                        summary_length = 0
                        summary_invalid = 0
                        summary_valid = 0
                        summary_valid_in_core = 0
                        strain_names = ['{best_ref}(ref)'.format(best_ref=best_ref)]
                        for ref_chrom, position_set in group_positions_set[species][group].items():
                            total_length = len(str(ref_records[ref_chrom].seq))
                            total_invalid = 0
                            total_valid = 0
                            total_valid_in_core = 0
                            for pos, ref_seq in enumerate(str(ref_records[ref_chrom].seq)):
                                # Adjust sequence to account for 0-based indexing
                                ref_seq = str(ref_records[ref_chrom].seq)[pos - 1]
                                # Initialise the valid and core variables
                                valid = True
                                core = True
                                # Determine if the position is present in the set of all SNV positions of the sample
                                if pos in position_set:
                                    snp_summary_body += '{ref_chrom}\t'.format(ref_chrom=ref_chrom)
                                    # If the position is present in the filter_reasons dictionary, it is neither a core,
                                    # nor a valid position
                                    try:
                                        # Add the reasons that the position is invalid to the validity string
                                        validity = 'invalid\t{reason}'\
                                            .format(reason=';'.join(filter_reasons[species][group][ref_chrom][pos]))
                                        valid = False
                                        core = False
                                    except KeyError:
                                        # If the position is not in the dictionary of identical postions, it is both
                                        # valid and core
                                        if pos not in ident_group_positions[species][group][ref_chrom]:
                                            validity = 'valid\t'
                                        # Otherwise, the SNV allele likely didn't have a high enough fraction, so the
                                        # position is invalid, but core
                                        else:
                                            validity = 'invalid\tmajority reference call'
                                            valid = False
                                            core = True
                                    sequence_string = str()
                                    for strain_name, chrom_dict in sorted(strain_dict.items()):
                                        # Don't process the reference strain
                                        # if strain_name != best_ref:
                                        # Add the strain name to the list of all strains
                                        if strain_name not in strain_names:
                                            strain_names.append(strain_name)
                                        try:
                                            sequence = chrom_dict[ref_chrom][pos]
                                            # Append each base to the string
                                            sequence_string += '{seq}\t'.format(seq=sequence)
                                        except KeyError:
                                            # If the position isn't in chrom_dict, it is either because it is
                                            # identical to the reference, or it was deleted
                                            try:
                                                # If it was deleted, add a -, otherwise, use the reference base
                                                sequence_dict = strain_parsed_vcf_dict[strain_name][ref_chrom][pos]
                                                if sequence_dict['FILTER'] == 'DELETION':
                                                    sequence_string += '-\t'
                                                else:
                                                    sequence_string += '{seq}\t'\
                                                        .format(seq=ref_seq)
                                            except KeyError:
                                                sequence_string += '{seq}\t' \
                                                    .format(seq=ref_seq)
                                    snp_summary_body += '{pos}\t{validity}\t{ref_seq}\t{seq_string}\n'\
                                        .format(pos=pos,
                                                validity=validity,
                                                ref_seq=ref_seq,
                                                seq_string=sequence_string)
                                # If the position isn't in the set of sample SNVs, check its status
                                else:
                                    for strain_name, chrom_dict in strain_dict.items():
                                        # If the position is missing, then it is not a core position
                                        try:
                                            sequence_dict = strain_parsed_vcf_dict[strain_name][ref_chrom][pos]
                                            if sequence_dict['FILTER'] == 'DELETION':
                                                core = False
                                        except KeyError:
                                            pass
                                    # Density filtering check. A density-filtered position is neither valid or core
                                    if pos in filtered_group_positions[species][group][ref_chrom]:
                                        valid = False
                                        core = False
                                    # Calculated mask. A masked position is neither valid nor core
                                    try:
                                        calculated_mask = mask_pos_dict[species][group][ref_chrom]
                                        if pos in calculated_mask:
                                            valid = False
                                            core = False
                                    except KeyError:
                                        pass
                                    # User-supplied mask
                                    try:
                                        user_mask = supplied_mask_pos_dict[species][group][ref_chrom]
                                        if pos in user_mask:
                                            valid = False
                                            core = False
                                    except KeyError:
                                        pass
                                # Determine if the position is either valid or core
                                if not valid:
                                    total_invalid += 1
                                else:
                                    total_valid += 1
                                    if core:
                                        total_valid_in_core += 1
                            # Calculate the necessary values
                            percent_valid_in_core = '{:.2f}'.format(total_valid_in_core/total_valid*100)
                            percent_total_valid_in_core = '{:.2f}'.format(total_valid_in_core/total_length*100)
                            summary_length += total_length
                            summary_invalid += total_invalid
                            summary_valid += total_valid
                            summary_valid_in_core += total_valid_in_core
                            pos_summary_body += '{ref_chrom}\t{total_length}\t{total_invalid}\t{total_valid}\t' \
                                                '{total_valid_in_core}\t{percent_valid_in_core}\t' \
                                                '{percent_total_valid_in_core}\n'\
                                .format(ref_chrom=ref_chrom,
                                        total_length=total_length,
                                        total_invalid=total_invalid,
                                        total_valid=total_valid,
                                        total_valid_in_core=total_valid_in_core,
                                        percent_valid_in_core=percent_valid_in_core,
                                        percent_total_valid_in_core=percent_total_valid_in_core)
                    summary_percent_valid_in_core = '{:.2f}'.format(summary_valid_in_core / summary_valid * 100)
                    summary_percent_total_valid_in_core = '{:.2f}'.format(summary_valid_in_core / summary_length * 100)
                    pos_summary_body += 'summary\t{summary_length}\t{summary_invalid}\t{summary_valid}\t' \
                                        '{summary_valid_in_core}\t{summary_percent_valid_in_core}\t' \
                                        '{summary_percent_total_valid_in_core}\n' \
                        .format(summary_length=summary_length,
                                summary_invalid=summary_invalid,
                                summary_valid=summary_valid,
                                summary_valid_in_core=summary_valid_in_core,
                                summary_percent_valid_in_core=summary_percent_valid_in_core,
                                summary_percent_total_valid_in_core=summary_percent_total_valid_in_core)
                    snp_summary_header += '{strain_name}\n'.format(strain_name='\t'.join(strain_names))
                    snp_summary.write(snp_summary_header)
                    snp_summary.write(snp_summary_body)
                    pos_summary.write(pos_summary_header)
                    pos_summary.write(pos_summary_body)

    @staticmethod
    def run_fasttree(group_fasta_dict, strain_consolidated_ref_dict, strain_groups, logfile):
        """
        Create maximum-likelihood tree using FastTree
        :param group_fasta_dict: type DICT: Dictionary of species code: group name: FASTA file created for the group
        :param strain_consolidated_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param logfile: type STR: Absolute path to logfile basename
        :return: species_group_trees: Dictionary of species code: group name: dictionary of tree type: absolute path
        to FastTree output tree
        """
        # Initialise a dictionary to store the absolute paths of the output trees
        species_group_trees = dict()
        for species, group_dict in group_fasta_dict.items():
            # Initialise the species key in the dictionary if necessary
            if species not in species_group_trees:
                species_group_trees[species] = dict()
            for strain_name, best_ref in strain_consolidated_ref_dict.items():
                for group, fasta_file in group_dict.items():
                    if group in strain_groups[strain_name]:
                        # Initialise the group key in the dictionary if necessary
                        if group not in species_group_trees[species]:
                            species_group_trees[species][group] = dict()
                        # Set the path of the working dir
                        output_dir = os.path.dirname(fasta_file)
                        best_tree = os.path.join(output_dir,
                                                 'best_tree.tre'
                                                 .format(species=species,
                                                         group=group))
                        species_group_trees[species][group]['best_tree'] = best_tree
                        # Create a system call to FastTree
                        fasttree_cmd = 'FastTree -gamma  -nt {fasta_file} > {tree_file}' \
                            .format(fasta_file=fasta_file,
                                    tree_file=best_tree)
                        # Run the system call if the output best tree doesn't already exist
                        if not os.path.isfile(best_tree):
                            out, err = run_subprocess(command=fasttree_cmd)
                            # Write the stdout and stderr to the main logfiles
                            write_to_logfile(out=out,
                                             err=err,
                                             logfile=logfile)
        return species_group_trees

    @staticmethod
    def parse_tree_order(species_group_trees):
        """
        Extract the order of the strains from the phylogenetic trees using ete3
        :param species_group_trees: type DICT: Dictionary of Dictionary of species code: group name: dictionary of
        tree type: absolute path to output tree
        :return: species_group_order_dict: Dictionary of species code: group name: list of ordered strains
        """
        # Initialise the dictionary to store the ordered list of strain
        species_group_order_dict = dict()
        for species, group_dict in species_group_trees.items():
            species_group_order_dict[species] = dict()
            for group, options_dict in group_dict.items():
                species_group_order_dict[species][group] = list()
                for tree_type, tree_file in options_dict.items():
                    # Only extract the order from the best trees
                    if tree_type == 'best_tree':
                        # Load the tree file with ete3
                        tree = Tree(tree_file)
                        # Traverse the tree using the 'postorder' strategy: 1) Traverse the left subtree,
                        # 2) Traverse the right subtree, 3) Visit the root
                        for node in tree.traverse('postorder'):
                            # When the node.name is present, append the node name to the list of names
                            if node.name:
                                species_group_order_dict[species][group].append(node.name)
        return species_group_order_dict

    @staticmethod
    def copy_trees(species_group_trees, tree_path):
        """
        Copy the output trees to a common folder
        :param species_group_trees: type DICT: Dictionary of Dictionary of species code: group name: dictionary of
        tree type: absolute path to output tree
        :param tree_path: type STR: Absolute path to folder into which tree files are to be copied
        """
        # Create the tree_path folder
        make_path(tree_path)
        for species, group_dict in species_group_trees.items():
            for group, options_dict in group_dict.items():
                for tree_type, tree_file in options_dict.items():
                    # Extract the file name from the path
                    tree_name = os.path.basename(tree_file)
                    # Set the name of the destination file
                    destination_file = os.path.join(tree_path, tree_name)
                    # Copy the file to the destination folder
                    if not os.path.isfile(destination_file):
                        shutil.copyfile(src=tree_file,
                                        dst=destination_file)

    @staticmethod
    def prokka(reference_strain_dict, logfile):
        """
        Run prokka to annotate the reference genome
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param logfile: type STR: Absolute path to the logfile
        """
        # Prepare command
        for strain_name, ref_fasta in reference_strain_dict.items():
            ref_path = os.path.dirname(ref_fasta)
            ref_strain = os.path.splitext(os.path.basename(ref_fasta))[0]
            gbk_file = ref_fasta.replace('.fasta', '.gbk')
            cmd = 'prokka --force --outdir {output_folder} --prefix {prefix} {ref_file}' \
                .format(output_folder=ref_path,
                        prefix=ref_strain,
                        ref_file=ref_fasta)
            if not os.path.isfile(gbk_file):
                out, err = run_subprocess(command=cmd)
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=cmd,
                                                           out=out),
                                 err=err,
                                 logfile=logfile)

    @staticmethod
    def load_filter_file(reference_strain_dict, strain_best_ref_dict):
        """
        Load the Excel files containing curated lists of locations or ranges of locations in the reference genome
        that must be filtered prior to performing phylogenetic analyses
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param strain_best_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :return: filter_dict: Dictionary of reference file: group name: locations to filter
        """
        # Initialise a dictionary to store the locations to filter
        filter_dict = dict()
        for strain_name, best_ref_path in reference_strain_dict.items():
            # Extract the name of the best reference from the dictionary e.g. NC_017250.1
            best_ref = strain_best_ref_dict[strain_name]
            # Set the name of the Excel file storing the regions to filter
            filter_file = os.path.join(best_ref_path, 'Filtered_Regions.xlsx')
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
    def load_genbank_file_multiple(reference_strain_dict, strain_best_ref_set_dict):
        """
        Use SeqIO to parse the best reference genome GenBank file for annotating SNP locations
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param strain_best_ref_set_dict: type DICT: Dictionary of strain name: set of strain-specific reference genomes
        :return: full_best_ref_gbk_dict: Dictionary of best ref: ref position: SeqIO parsed GenBank file-sourced
        records from the closest reference genome for that position
        """
        # Initialise a dictionary to store the SeqIO parsed GenBank files
        best_ref_gbk_dict = dict()
        for strain_name, best_ref_path in reference_strain_dict.items():
            # Extract the species code from the dictionary
            best_ref_set = strain_best_ref_set_dict[strain_name]
            for best_ref in best_ref_set:
                gbk_file = glob(os.path.join(best_ref_path, '{br}*.gbk'.format(br=os.path.splitext(best_ref)[0])))[0]
                # Only parse the file if it has not already been parsed
                if best_ref not in best_ref_gbk_dict:
                    # Use SeqIO to first parse the GenBank file, and convert the parsed object to a dictionary
                    gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
                    # Add the GenBank dictionary to the best reference-specific dictionary
                    best_ref_gbk_dict[best_ref] = gbk_dict
        # Initialise a dictionary to store all the positions, and associated information
        full_best_ref_gbk_dict = dict()
        for best_ref, gbk_dict in best_ref_gbk_dict.items():
            full_best_ref_gbk_dict[best_ref] = dict()
            for ref_name, record in gbk_dict.items():
                for feature in record.features:
                    # Ignore the full record
                    if feature.type != 'source':
                        # Iterate through the full length of the feature, and add each position to the dictionary
                        for i in range(int(feature.location.start), int(feature.location.end) + 1):
                            full_best_ref_gbk_dict[best_ref][i] = feature
        return full_best_ref_gbk_dict

    @staticmethod
    def load_genbank_file_single(reference_strain_dict):
        """
        Use SeqIO to parse the best reference genome GenBank file for annotating SNP locations
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :return: full_best_ref_gbk_dict: Dictionary of best ref: ref position: SeqIO parsed GenBank file-sourced
        records from the closest reference genome for that position
        """
        # Initialise a dictionary to store the SeqIO parsed GenBank files
        best_ref_gbk_dict = dict()
        full_best_ref_gbk_dict = dict()
        for strain_name, best_ref_file in reference_strain_dict.items():
            best_ref = os.path.splitext(best_ref_file)[0]
            gbk_file = best_ref_file.replace('.fasta', '.gbk')
            # Only parse the file if it has not already been parsed
            if best_ref not in best_ref_gbk_dict:
                # Use SeqIO to first parse the GenBank file, and convert the parsed object to a dictionary
                gbk_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, "genbank"))
                # Add the GenBank dictionary to the best reference-specific dictionary
                best_ref_gbk_dict[best_ref] = gbk_dict
        for best_ref, gbk_dict in best_ref_gbk_dict.items():
            full_best_ref_gbk_dict[best_ref] = dict()
            for ref_name, record in gbk_dict.items():
                if record.name not in full_best_ref_gbk_dict:
                    full_best_ref_gbk_dict[record.name] = dict()
                for feature in record.features:
                    # Ignore the full record
                    if feature.type != 'source':
                        # Iterate through the full length of the feature, and add each position to the dictionary
                        for i in range(int(feature.location.start), int(feature.location.end) + 1):
                            full_best_ref_gbk_dict[record.name][i] = feature
        return full_best_ref_gbk_dict

    @staticmethod
    def annotate_snps(group_strain_snp_sequence, full_best_ref_gbk_dict, strain_best_ref_set_dict, ref_snp_positions):
        """
        Use GenBank records to annotate each SNP with 'gene', 'locus', and 'product' details
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :param full_best_ref_gbk_dict: type DICT: Dictionary of best ref: ref position: SeqIO parsed GenBank
        file-sourced records from the closest reference genome for that position
        :param strain_best_ref_set_dict: type DICT: Dictionary of strain name: set of strain-specific reference genomes
        :param ref_snp_positions: type DICT: Dictionary of reference chromosome name: absolute position: reference base
         call
        :return: species_group_annotated_snps_dict: Dictionary of species code: group name: reference chromosome:
        reference position: annotation dictionary
        """
        # Initialise a dictionary to store the annotations for the group-specific SNPs
        species_group_annotated_snps_dict = dict()
        for species, group_dict in group_strain_snp_sequence.items():
            # Initialise the key in the dictionary if necessary
            if species not in species_group_annotated_snps_dict:
                species_group_annotated_snps_dict[species] = dict()
            for group, strain_dict in group_dict.items():
                if group not in species_group_annotated_snps_dict[species]:
                    species_group_annotated_snps_dict[species][group] = dict()
                for strain_name, ref_dict in strain_dict.items():
                    for ref_chrom, pos_dict in ref_dict.items():
                        if strain_name in strain_best_ref_set_dict:
                            # As the reference only needs to be added to the dictionary once, ensure that it is not
                            # present before continuing
                            if ref_chrom not in species_group_annotated_snps_dict[species][group]:
                                species_group_annotated_snps_dict[species][group][ref_chrom] = dict()
                            # Unpack the dictionary using the ref_chrom as the key
                            gbk_pos_dict = full_best_ref_gbk_dict[ref_chrom]
                            for pos in pos_dict:
                                # Ensure that the position is in the specific chromosome being considered e.g.
                                # 'NC_017250.1' vs 'NC_017251.1'
                                if pos in ref_snp_positions[ref_chrom]:
                                    species_group_annotated_snps_dict[species][group][ref_chrom][pos] = dict()
                                    # Non-coding regions will not be present in the dictionary
                                    try:
                                        # Extract the SeqIO-parsed GenBank feature from the GenBank record
                                        # dictionary
                                        feature = gbk_pos_dict[pos]
                                        # Populate the 'locus', 'gene', and 'product' key: value pairs. Extract
                                        # the values from the feature.qualifiers OrderedDict list
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['locus'] = \
                                            feature.qualifiers['locus_tag'][0]
                                        # Not all features have the 'gene' key. Add 'None' if this is the case
                                        try:
                                            gene = feature.qualifiers['gene'][0]
                                        except KeyError:
                                            gene = 'None'
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['gene'] \
                                            = gene
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['product'] \
                                            = feature.qualifiers['product'][0]
                                        # Extract the location of the coding sequence e.g. [48999:49497](+)
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['location'] \
                                            = feature.location
                                    # Populate negative key: value pairs if the position is not in the dictionary
                                    except KeyError:
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['locus'] \
                                            = 'None'
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['gene'] \
                                            = 'None'
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['product'] \
                                            = 'None'
                                        species_group_annotated_snps_dict[species][group][ref_chrom][pos]['location'] \
                                            = 'None'
        return species_group_annotated_snps_dict

    @staticmethod
    def determine_snp_number(group_strain_snp_sequence, species_group_best_ref):
        """
        Determine the number of strains that have a SNP at each group-specific position
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :return: species_group_snp_num_dict: Dictionary of species code: group name: reference chromosome:
        position: number of strains that have a SNP at that position
        """
        # Initialise a dictionary to store the number of strains that have a SNP for each group-specific position
        species_group_snp_num_dict = dict()
        for species, group_dict in group_strain_snp_sequence.items():
            # Set the species key
            species_group_snp_num_dict[species] = dict()
            for group, strain_dict in group_dict.items():
                # Set the group key
                species_group_snp_num_dict[species][group] = dict()
                # Extract the name of the reference genome from the species_group_best_ref dictionary using the species
                # code and the group name
                best_ref = species_group_best_ref[species][group]
                best_ref_dict = strain_dict[best_ref]
                for strain_name, ref_dict in sorted(strain_dict.items()):
                    # The reference genome should not be considered when counting SNP positions
                    # if strain_name != best_ref:
                    for ref_chrom, sequence_dict in ref_dict.items():
                        if ref_chrom not in species_group_snp_num_dict[species][group]:
                            species_group_snp_num_dict[species][group][ref_chrom] = dict()
                        for pos, sequence in sequence_dict.items():
                            # If the base at the current position is not the same as the reference genome
                            # (strain_dict[best_ref][ref_chrom][pos]), the position is a SNP
                            if sequence != best_ref_dict[ref_chrom][pos]:
                                # Initialise the number of SNPs at a position to 1
                                if pos not in species_group_snp_num_dict[species][group][ref_chrom]:
                                    species_group_snp_num_dict[species][group][ref_chrom][pos] = 1
                                # Otherwise increment the number of SNPs observed at the position
                                else:
                                    species_group_snp_num_dict[species][group][ref_chrom][pos] += 1
        return species_group_snp_num_dict

    @staticmethod
    def determine_aa_sequence(group_strain_snp_sequence, species_group_best_ref, strain_parsed_vcf_dict,
                              species_group_annotated_snps_dict, reference_strain_dict, species_group_snp_num_dict,
                              iupac):
        """
        Determine the amino acid residue at the SNP position
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param species_group_annotated_snps_dict: type DICT: Dictionary of species code: group name: reference
        :param reference_strain_dict: type DICT: Dictionary of strain name: absolute path to reference genome
        :param species_group_snp_num_dict: type DICT: Dictionary of species code: group name: reference chromosome:
        position: number of strains that have a SNP at that position
        :param iupac: type DICT: Dictionary of degenerate code: nucleotides included in group
        :return translated_snp_residue_dict: Dictionary of species: group: ref_chrom: pos: {pos-specific sequence info}
        """
        # Initialise a dictionary to store the annotations for the group-specific SNPs
        translated_snp_residue_dict = dict()
        ref_translated_snp_residue_dict = dict()
        for species, group_dict in group_strain_snp_sequence.items():
            # Initialise the key in the dictionary if necessary
            if species not in translated_snp_residue_dict:
                translated_snp_residue_dict[species] = dict()
                ref_translated_snp_residue_dict[species] = dict()
            for group, strain_dict in group_dict.items():
                # Extract the name of the reference genome from the species_group_best_ref dictionary using the species
                # code and the group name
                best_ref = species_group_best_ref[species][group]
                # Use SeqIO to parse all the records in the reference FASTA file
                ref_records = SeqIO.to_dict(SeqIO.parse(reference_strain_dict[best_ref], 'fasta'))
                # Initialise the group key in the dictionary as required
                if group not in translated_snp_residue_dict[species]:
                    translated_snp_residue_dict[species][group] = dict()
                    ref_translated_snp_residue_dict[species][group] = dict()
                for strain_name, ref_dict in strain_dict.items():
                    # if strain_name != best_ref:
                    # Add the strain_name to the nested dictionary
                    if strain_name not in translated_snp_residue_dict[species][group]:
                        translated_snp_residue_dict[species][group][strain_name] = dict()
                    for ref_chrom in ref_dict:
                        # Initialise the ref_chrom key in the dictionary as required
                        if ref_chrom not in translated_snp_residue_dict[species][group][strain_name]:
                            translated_snp_residue_dict[species][group][strain_name][ref_chrom] = dict()
                            ref_translated_snp_residue_dict[species][group][ref_chrom] = dict()
                            # Iterate through all the SNP positions found in the group
                            for pos in species_group_snp_num_dict[species][group][ref_chrom]:
                                # Extract the SNP sequence from the ref_dict
                                try:
                                    snp_seq = ref_dict[ref_chrom][pos]
                                # If the strain does not have an entry at that position, it matches the reference
                                # sequence
                                except KeyError:
                                    snp_seq = strain_dict[best_ref][ref_chrom][pos]
                                # Extract the location of the coding sequence from the reference chromosome
                                # containing the SNP position
                                location = \
                                    species_group_annotated_snps_dict[species][group][ref_chrom][pos]['location']
                                # rRNA products are not coding
                                product = \
                                    species_group_annotated_snps_dict[species][group][ref_chrom][pos]['product']
                                # 'None' locations and rRNA products do not correspond to coding regions.
                                # Deletions do not have a corresponding amino acid sequence
                                if location != 'None' and 'ribosomal RNA' not in product and snp_seq != '-':
                                    # Extract the coding sequence of reference genome containing the SNP position
                                    ref_nt_seq = ref_records[ref_chrom].seq[location.start:location.end]
                                    # If the SNP is a degenerate base, parse the 'ALT' entry from the gVCF file
                                    if snp_seq in iupac:
                                        # The 'ALT' entry contains the alternate base, and the reference: C,<*>
                                        alt_dict = strain_parsed_vcf_dict[strain_name][ref_chrom][pos]['ALT']
                                        # Split the entry on the comma, and update the snp_seq variable with
                                        # the alternate base from C,<*>
                                        snp_seq = alt_dict.split(',')[0]
                                    # Create a variable to store the raw SNP sequence
                                    alt_seq = snp_seq
                                    # If the coding sequence containing the SNP position is on the positive strand,
                                    # the location of the SNP is calculated by subtracting the defined start pos
                                    # of the coding sequence (+1 due to 0-based indexing) from the global SNP pos
                                    if location.strand == 1:
                                        # e.g. 49340 (pos) - 48999 (location.start) + 1 = 342 (snp_loc)
                                        '''
                                        CDS 49000..49497
                                        /locus_tag="PELMLHNL_01252"
                                        '''
                                        # + 1
                                        snp_loc = pos - location.start
                                    # If the coding sequence containing the SNP position is on the negative strand,
                                    # the location of the SNP is calculated by subtracting the global SNP position
                                    # from the end of the defined end position of the coding sequence (don't need
                                    # to worry about 0-based indexing)
                                    else:
                                        # The reference sequence is treated as the reverse complement
                                        ref_nt_seq = ref_nt_seq.reverse_complement()
                                        # e.g. 3216 (location.end) - 2118 (pos) = 1098 (snp_loc)
                                        '''
                                        CDS complement(1285..3216)
                                        /gene="ccmF_1"
                                        '''
                                        snp_loc = location.end - pos
                                        # Convert the snp sequence to be the reverse complement
                                        snp_seq = Seq(snp_seq)
                                        snp_seq = str(snp_seq.reverse_complement())
                                    # Strings are not mutable, so convert the coding sequence to a list in
                                    # anticipation of modifying the sequence at the SNP position
                                    snp_nt_seq = list(ref_nt_seq)
                                    # Update the list at the SNP index with the new sequence: e.g. at position 342,
                                    # the base 'T' is substituted with a 'G'
                                    # Convert the list to be a Biopython Seq object
                                    snp_nt_seq = Seq(''.join(snp_nt_seq))
                                    # Determine the position of the amino acid residue within the translated
                                    # sequence by dividing the location in the nucleic acid sequence by three and
                                    # rounding down e.g. 342 (snp_loc) / 3 = 114
                                    aa_loc = math.floor(snp_loc / 3)
                                    # Populate the dictionary with all the necessary data
                                    translated_snp_residue_dict[species][group][strain_name][ref_chrom][pos] = {
                                        'snp_nt_seq_raw': snp_seq,
                                        'snp_nt_seq_alt': alt_seq,
                                        'snp_nt_seq_cds': snp_nt_seq[snp_loc],
                                        'snp_aa_seq_cds': str(snp_nt_seq.translate())[aa_loc],
                                        'ref_nt_seq_raw': ref_records[ref_chrom].seq[location.start:location.end]
                                        [snp_loc],
                                        'ref_nt_seq_cds': ref_nt_seq[snp_loc],
                                        'ref_aa_seq_cds': str(ref_nt_seq.translate())[aa_loc],
                                        'cds_strand': location.strand
                                    }
                                    ref_translated_snp_residue_dict[species][group][ref_chrom][pos] = {
                                        'ref_nt_seq_raw': ref_records[ref_chrom].seq[location.start:location.end]
                                        [snp_loc],
                                        'ref_nt_seq_cds': ref_nt_seq[snp_loc],
                                        'ref_aa_seq_cds': str(ref_nt_seq.translate())[aa_loc],
                                        'cds_strand': location.strand
                                    }
        return translated_snp_residue_dict, ref_translated_snp_residue_dict

    @staticmethod
    def create_snp_matrix(species_group_best_ref, group_strain_snp_sequence, matrix_path):
        """
        Create a matrix of the pairwise SNPs between each strain
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :param matrix_path: type STR: Absolute path to folder in which matrix files are to be created
        """
        # Create the matrix_path if required
        make_path(matrix_path)
        # Initialise a dictionary to store the pairwise SNPs between all the strains
        snp_matrix = dict()
        # Initialise a set to track all SNV positions
        positions = dict()
        for species, group_dict in group_strain_snp_sequence.items():
            for group, strain_dict in group_dict.items():
                # Extract the name of the reference strain
                consolidated_ref = species_group_best_ref[species][group]
                # Iterate through the dictionary twice - get a strain name to be used to compare against the other
                # strains
                for compare_strain in group_strain_snp_sequence[species][group]:
                    # Initialise the SNP dictionary with the target name
                    snp_matrix[compare_strain] = dict()
                    if compare_strain not in positions:
                        positions[compare_strain] = dict()
                    # The target strain will not have any SNPs against itself, but this needs to be recorded for
                    # creating the square matrix file
                    snp_matrix[compare_strain][compare_strain] = 0
                    # Iterate through the dictionary again to get the query strain
                    for strain_name, ref_dict in group_strain_snp_sequence[species][group].items():
                        # Query strain must be different from the target strain
                        if strain_name != compare_strain:
                            # Initialise the keys in the snp_matrix dictionary as required
                            if strain_name not in snp_matrix[compare_strain]:
                                snp_matrix[compare_strain].update({strain_name: 0})
                            if strain_name not in positions[compare_strain]:
                                positions[compare_strain][strain_name] = dict()
                            # Iterate through the reference chromosomes, and position dictionaries nested in the
                            # reference dictionary
                            for ref_chrom, pos_dict in ref_dict.items():
                                if ref_chrom not in positions[compare_strain][strain_name]:
                                    positions[compare_strain][strain_name][ref_chrom] = dict()
                                # Extract the position and sequence at that position on the reference chromosome
                                for pos, seq in pos_dict.items():
                                    # Add the position to the dictionary
                                    positions[compare_strain][strain_name][ref_chrom][pos] = seq
                                    # Determine the sequence of the target strain at this position
                                    try:
                                        # If this position is in the dictionary, extract the sequence from the
                                        # dictionary
                                        compare_seq = \
                                            group_strain_snp_sequence[species][group][compare_strain][ref_chrom][pos]
                                    except KeyError:
                                        # If this position is not in the dictionary, then it is because it matches the
                                        # reference sequence, so that should be used
                                        compare_seq = \
                                            group_strain_snp_sequence[species][group][consolidated_ref][ref_chrom][pos]
                                    # If the target and query sequences do not match, increment the snp_matrix counter
                                    if seq != compare_seq:
                                        snp_matrix[compare_strain][strain_name] += 1
                # Determine which positions are missing from certain strain combinations
                for compare_strain, strain_name_dict in positions.items():
                    for strain_name, ref_chrom_dict in strain_name_dict.items():
                        for ref_chrom, pos_dict in ref_chrom_dict.items():
                            for pos, seq in pos_dict.items():
                                # Ignore this position if it is in the dictionary of positions for the "other" strain
                                if pos in positions[strain_name][compare_strain][ref_chrom]:
                                    continue
                                # Extract the reference sequence for this position
                                ref_seq = group_strain_snp_sequence[species][group][consolidated_ref][ref_chrom][pos]
                                # Only increment the count if this position is different than the reference sequence
                                if seq != ref_seq:
                                    snp_matrix[strain_name][compare_strain] += 1
                # Write the snp_matrix dictionary to a .csv file
                with open(os.path.join(matrix_path, 'snv_matrix.tsv'), 'w') as matrix_file:
                    # Initialise variables to store the header, and body strings
                    header = 'Strain\t'
                    body = str()
                    # A count that is incremented for each query strain, so that a list can be accessed by the
                    # current index
                    count = 0
                    # Iterate through all the strains in the snp_matrix dictionary
                    for ref, compare_dict in sorted(snp_matrix.items()):
                        # The reference strain is special, so add '(ref)' to the strain name
                        if ref == consolidated_ref:
                            ref += ' (ref)'
                        # Update the header with the strain name
                        header += '{ref}\t'.format(ref=ref)
                        # Boolean on whether the name of the current query strain has been printed yet
                        print_strain = True
                        # Iterate through all the strains compared to this target strain
                        for query, num_snps in sorted(compare_dict.items()):
                            # Only include the strain name if it hasn't already been added to the strain column
                            if print_strain:
                                # Extract the strain name from the dictionary based on its sorted position
                                strain = sorted(compare_dict.items())[count][0]
                                # Add the '(ref)' string when processing the reference strain
                                if strain == consolidated_ref:
                                    strain += ' (ref)'
                                body += '{query}\t'.format(query=strain)
                                # Set the boolean to False, so that the strain name isn't added again in this row
                                print_strain = False
                                # Increment the count for the next sample name
                                count += 1
                            # Add the number of SNPs between the query and target genomes
                            body += '{num_snps}\t'.format(num_snps=num_snps)
                        # Add a newline character to separate each query strain
                        body += '\n'
                    header += '\n'
                    matrix_file.write(header)
                    matrix_file.write(body)

    @staticmethod
    def rank_snps(species_group_snp_num_dict):
        """
        Sort the group-specific SNP positions based on how prevalent they are in the group
        :param species_group_snp_num_dict: type DICT: Dictionary of species code: group name: reference chromosome:
        position: number of strains that have a SNP at that position
        :return: species_group_snp_rank: Dictionary of species code: group name: num strains with SNP:
        reference chromosome: SNP position
        :return: species_group_num_snps: Dictionary of species code: group name: total number of group-specific SNP
        positions
        """
        # Initialise dictionaries to store ranked SNP positions, and the total number of group-specific SNPs
        species_group_snp_rank = dict()
        species_group_num_snps = dict()
        # species_group_snp_num_dict[species][group][ref_chrom][pos]
        for species, group_dict in species_group_snp_num_dict.items():
            # Initialise the species key
            species_group_snp_rank[species] = dict()
            species_group_num_snps[species] = dict()
            for group, ref_dict in group_dict.items():
                # Initialise the group key
                species_group_snp_rank[species][group] = dict()
                species_group_num_snps[species][group] = int()
                for ref_chrom, pos_dict in ref_dict.items():
                    for pos, num_snps in pos_dict.items():
                        # Initialise the number of SNPs in the dictionary
                        if num_snps not in species_group_snp_rank[species][group]:
                            species_group_snp_rank[species][group][num_snps] = dict()
                        # Finally, initialise the reference chromosome in the dictionary
                        if ref_chrom not in species_group_snp_rank[species][group][num_snps]:
                            species_group_snp_rank[species][group][num_snps][ref_chrom] = list()
                        # Add the position to the dictionary under the current number of SNPs and reference chromosome
                        species_group_snp_rank[species][group][num_snps][ref_chrom].append(pos)
                        species_group_num_snps[species][group] += 1
        return species_group_snp_rank, species_group_num_snps

    @staticmethod
    def sort_snps(species_group_order_dict, species_group_snp_rank, species_group_best_ref,
                  group_strain_snp_sequence):
        """
        Sort the group-specific SNP positions on two criteria 1) Number of strains with a SNP at that position,
        2) Based on phylogenetic tree topology (a SNP that is only present in certain strains will only be added to
        the list once a strain in which it is present is assessed)
        :param species_group_order_dict: type DICT: Dictionary of species code: group name: list of ordered strains
        :param species_group_snp_rank: type DICT: Dictionary of species code: group name: num strains with SNP:
        reference chromosome: SNP position
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :return: species_group_sorted_snps: Dictionary of species code: group name: reference chromosome: ordered
        list of SNP positions
        """
        # Initialise a dictionary to store the group-specific SNP order
        species_group_sorted_snps = dict()
        for species, group_dict in species_group_order_dict.items():
            # Initialise the species key
            species_group_sorted_snps[species] = dict()
            for group, ordered_strain_list in group_dict.items():
                # Initialise the group key
                species_group_sorted_snps[species][group] = dict()
                # Extract the name of the reference genome from the species_group_best_ref genome
                best_ref = species_group_best_ref[species][group]
                # Extract the number of group-specific SNPs from the reverse-sorted dictionary (more SNPs first)
                for num_snps, ref_dict in sorted(species_group_snp_rank[species][group].items(), reverse=True):
                    if num_snps not in species_group_sorted_snps[species][group]:
                        species_group_sorted_snps[species][group][num_snps] = dict()
                    for strain_name in ordered_strain_list:
                        # Don't need to look at the reference genome when finding SNPs
                        # if strain_name != best_ref:
                        for ref_chrom, pos_list in ref_dict.items():
                            if ref_chrom not in species_group_sorted_snps[species][group][num_snps]:
                                species_group_sorted_snps[species][group][num_snps][ref_chrom] = list()
                            for pos in pos_list:
                                try:
                                    sequence = \
                                        group_strain_snp_sequence[species][group][strain_name][ref_chrom][pos]
                                    ref_seq = group_strain_snp_sequence[species][group][best_ref][ref_chrom][pos]
                                    if sequence != ref_seq and pos not in \
                                            species_group_sorted_snps[species][group][num_snps][ref_chrom]:
                                        # Add the position to the list
                                        species_group_sorted_snps[species][group][num_snps][ref_chrom].append(pos)
                                # If the position is not present in group_strain_snp_sequence, it is because
                                # the strain does not have a SNP at that position
                                except KeyError:
                                    pass

        for species, group_dict in species_group_sorted_snps.items():
            for group, num_dict in group_dict.items():
                for num_snps, chrom_dict in num_dict.items():
                    for ref_chrom in chrom_dict:
                        chrom_dict[ref_chrom] = sorted(chrom_dict[ref_chrom])
        return species_group_sorted_snps

    @staticmethod
    def create_summary_table(species_group_sorted_snps, species_group_order_dict, species_group_best_ref,
                             group_strain_snp_sequence, species_group_annotated_snps_dict, translated_snp_residue_dict,
                             ref_translated_snp_residue_dict, species_group_num_snps, summary_path, molecule):
        """
        Create an Excel table that summarises the sorted SNP positions, and adds the annotations
        :param species_group_sorted_snps: type DICT: Dictionary of species code: group name: reference chromosome:
        ordered list of SNP positions
        :param species_group_order_dict: type DICT: Dictionary of species code: group name: list of ordered strains
        :param species_group_best_ref: type DICT: Dictionary of species code: group name: best ref
        :param group_strain_snp_sequence: type DICT: Dictionary of species code: group name: strain name:
        reference chromosome: position: strain-specific sequence
        :param species_group_annotated_snps_dict: type DICT: Dictionary of species code: group name: reference
        chromosome: reference position: annotation dictionary
        :param translated_snp_residue_dict: type DICT: Dictionary of species: group: strain_name: ref_chrom: pos:
        {pos-specific sequence info}
        :param ref_translated_snp_residue_dict: type DICT: Dictionary of species: group: ref_chrom: pos:
        {pos-specific sequence info}
        :param species_group_num_snps: type DICT: Dictionary of species code: group name: total number of
        group-specific SNP positions
        :param summary_path: type STR: Absolute path to folder in which summary reports are to be created
        :param molecule: type STR: String of whether the desired outputs are nucleotide (nt) or amino acid residue (aa)
        """
        for species, group_dict in species_group_order_dict.items():
            for group, ordered_strain_list in group_dict.items():
                # Extract the name of the reference genome from the species_group_best_ref genome
                consolidated_ref = species_group_best_ref[species][group]
                total_snps = species_group_num_snps[species][group]
                # Initialise a variable to store the current column for the report; each reference chromosome will
                # be added to the report, and cannot overwrite the previous results
                current_col = 0
                # Set the name of the summary table
                summary_table = os.path.join(summary_path, '{molecule}_snv_sorted_table.xlsx'
                                             .format(molecule=molecule))
                # Create an xlsxwriter workbook object
                wb = xlsxwriter.Workbook(summary_table)
                # Create a worksheet in the workbook
                ws = wb.add_worksheet()
                # Create all the necessary formats for the workbook
                header, courier, bold_courier, top_bold_courier, annotation, format_dict, ambiguous_format = \
                    TreeMethods.format_workbook(wb=wb,
                                                molecule=molecule)
                # Initialise a variable to store the current row under consideration
                row = 0
                # Merge all the cells in the first row from the second column until the column corresponding to
                # the length of the total number of group-specific SNP positions
                ws.merge_range(first_row=row,
                               first_col=1,
                               last_row=0,
                               last_col=total_snps,
                               data='SNV Position',
                               cell_format=top_bold_courier)
                # Adjust the width of the columns from the 2nd until the column corresponding to the total number
                # of SNPs to 2
                ws.set_column(first_col=1,
                              last_col=total_snps,
                              width=2)
                # Write the 'Strain' header
                ws.write_string(row=row,
                                col=0,
                                string='Strain',
                                cell_format=top_bold_courier)
                # Determine the height of the cells for the annotation results
                annotation_height = 0
                for num_snps, chrom_dict in species_group_sorted_snps[species][group].items():
                    for ref_chrom, snp_order in chrom_dict.items():
                        # Determine the height to use based on the group-specific SNP with the longest annotation
                        # (multiplied by 5 as determined by trial and error)
                        try:
                            annotation_list = list()
                            for pos in snp_order:
                                annotation_dict = species_group_annotated_snps_dict[species][group][ref_chrom][pos]
                                annotation_string = '{product};{gene};{locus}'.format(
                                    product=annotation_dict['product'],
                                    gene=annotation_dict['gene'],
                                    locus=annotation_dict['locus'])
                                annotation_list.append(annotation_string)
                            current_height = 5 * max([len(annotation) for annotation in annotation_list])
                            # Determine the height to use for the header. Each cell consists of the ref chromosome name
                            # and the SNP pos (e.g. NC_002945.4_1057) rotated 270 degrees, so the cell has a height
                            # equal to the length of the header multiplied by 6 (as determined by trial and error)
                            max_snp_length = max([len(str(snp)) for snp in snp_order])
                            snp_height = 6 * (max_snp_length + len(ref_chrom))
                            # Adjust the row height
                            ws.set_row(row=1,
                                       height=snp_height)
                        except ValueError:
                            current_height = annotation_height
                        # As there are multiple reference chromosomes to consider, ensure that the length of the longest
                        # annotation of all chromosomes is being used
                        annotation_height = current_height if current_height > annotation_height else annotation_height
                for num_snps, chrom_dict in species_group_sorted_snps[species][group].items():
                    for ref_chrom, snp_order in chrom_dict.items():
                        row = 1
                        # Extract the strain_dict (dictionary of strain name: pos: pos sequence) from the dictionary
                        strain_dict = group_strain_snp_sequence[species][group]
                        # Extract consolidated reference genome pos: sequence dictionary
                        ref_dict = strain_dict[consolidated_ref][ref_chrom]
                        # Set the width of the first column to be the longest of the following items: 1) the length
                        # of the longest strain name, 2) the length of the consolidated reference, 3) length of the word
                        # 'Annotation'; the longest hardcoded string in the column
                        strain_length = max([len(consolidated_ref + ' (ref)'), len(max(strain_dict)), 10])
                        ws.set_column(first_col=0,
                                      last_col=0,
                                      width=strain_length)

                        # Write the header consisting of the reference chromosome + '_' + SNP position
                        # (e.g. NC_002945.4_1057) for every ordered SNP
                        for i, entry in enumerate(snp_order):
                            ws.write_string(row=row,
                                            col=current_col + i + 1,
                                            string='{ref_chrom}_{entry}'.format(ref_chrom=ref_chrom,
                                                                                entry=entry),
                                            cell_format=header)
                        # Increment the row for the reference genome sequences
                        row += 1
                        # Add the name of the consolidated reference sequence to the 'Strain' column
                        ws.write_string(row=row,
                                        col=0,
                                        string=consolidated_ref + ' (ref)',
                                        cell_format=courier)
                        # Iterate through the sorted SNPs, and write the reference sequence for each position
                        for i, pos in enumerate(snp_order):
                            if molecule == 'nt':
                                ws.write_string(row=row,
                                                col=current_col + i + 1,
                                                string=ref_dict[pos],
                                                cell_format=bold_courier)
                            else:
                                # Determine the sequence at the current position, and the format to use to write the
                                # sequence in the report
                                sequence, base_format = TreeMethods.format_sequence(
                                    pos=pos,
                                    sequence_dict=ref_translated_snp_residue_dict[species][group][ref_chrom],
                                    ref_dict=ref_dict,
                                    bold_courier=bold_courier,
                                    format_dict=format_dict,
                                    ambiguous_format=ambiguous_format,
                                    molecule=molecule,
                                    ref=True)
                                ws.write_string(row=row,
                                                col=current_col + i + 1,
                                                string=sequence,
                                                cell_format=base_format)
                        # Increment the row
                        row += 1
                        # Freeze the panes, so that the row containing the reference sequence is at the bottom of the
                        # frozen pane, and the column with the strain names is always present
                        ws.freeze_panes(row=row,
                                        col=1)
                        # Add the strain-specific data to the table
                        for strain_name in ordered_strain_list:
                            # Continue to unpack the dictionary to obtain sequence_dict: a dictionary of pos: pos seq
                            try:
                                sequence_dict = strain_dict[strain_name][ref_chrom]
                            except KeyError:
                                sequence_dict = strain_dict[consolidated_ref][ref_chrom]
                            # Don't need to look at the reference genome when finding SNPs
                            # if strain_name != consolidated_ref:
                            # Write the strain name in the 'Strain' column
                            ws.write_string(row=row,
                                            col=0,
                                            string=strain_name,
                                            cell_format=courier)
                            for i, pos in enumerate(snp_order):
                                if molecule == 'nt':
                                    # Determine the sequence of the nucleotide at the current position, and any
                                    # special formatting required
                                    sequence, base_format = TreeMethods.format_sequence(
                                        pos=pos,
                                        sequence_dict=sequence_dict,
                                        ref_dict=ref_dict,
                                        bold_courier=bold_courier,
                                        format_dict=format_dict,
                                        ambiguous_format=ambiguous_format,
                                        molecule=molecule)
                                # For amino acid reports, use the translated_snp_residue_dict instead
                                else:
                                    sequence_dict = \
                                        translated_snp_residue_dict[species][group][strain_name][ref_chrom]
                                    sequence, base_format = TreeMethods.format_sequence(
                                        pos=pos,
                                        sequence_dict=sequence_dict,
                                        ref_dict=ref_dict,
                                        bold_courier=bold_courier,
                                        format_dict=format_dict,
                                        ambiguous_format=ambiguous_format,
                                        molecule=molecule)
                                # Write the sequence in the appropriate format
                                ws.write_string(row=row,
                                                col=current_col + i + 1,
                                                string=sequence,
                                                cell_format=base_format)
                            # Increment the row for each sequence
                            row += 1

                        # Add the string 'Annotation' to the 'Strain' column
                        ws.write_string(row=row,
                                        col=0,
                                        string='Annotation',
                                        cell_format=top_bold_courier)
                        # Set the row height to the previously calculated height
                        ws.set_row(row=row,
                                   height=annotation_height)
                        # Extract the 'product' annotations for the SNP from the species_group_annotated_snps_dict
                        for i, pos in enumerate(snp_order):
                            annotation_dict = species_group_annotated_snps_dict[species][group][ref_chrom][pos]
                            annotation_string = '{product};{gene};{locus}'.format(product=annotation_dict['product'],
                                                                                  gene=annotation_dict['gene'],
                                                                                  locus=annotation_dict['locus'])
                            ws.write_string(row=row,
                                            col=current_col + i + 1,
                                            string=annotation_string,
                                            cell_format=annotation)
                        current_col += len(snp_order)
                        # Increment the row
                        row += 1
                        # Set the final row to a height of 1
                        ws.set_row(row=row,
                                   height=1)
                # Close the workbook
                wb.close()

    @staticmethod
    def format_workbook(wb, molecule, font_size=8):
        """
        Create the required formats for the summary table
        :param wb: xlsxwriter workbook object
        :param molecule: type STR: Creating nucleotide or amino acid outputs
        :param font_size: type INT: Font size to use in creating the report
        :return: workbook formats
        """
        # Add a bold format for header cells. Using a monotype font size 8, rotated 270 degrees
        header = wb.add_format({'bold': True,
                                'font_name': 'Courier New',
                                'font_size': font_size,
                                'rotation': '90'})
        # Format for data cells. Monotype, size 8, top vertically justified
        courier = wb.add_format({'font_name': 'Courier New',
                                 'font_size': font_size,
                                 'align': 'top'})
        # Bold courier format
        bold_courier = wb.add_format({'font_name': 'Courier New',
                                      'font_size': font_size,
                                      'bold': True})
        # Bold courier top aligned
        top_bold_courier = wb.add_format({'font_name': 'Courier New',
                                          'font_size': font_size,
                                          'bold': True,
                                          'align': 'top'})
        # Create a format for the annotations: bold, blue text, rotated -90 degrees
        annotation = wb.add_format({'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size,
                                    'font_color': '#0A028C',
                                    'rotation': '-90',
                                    'align': 'top'})
        if molecule == 'nt':
            # Create a dictionary to store formats for the difference nucleotides
            format_dict = {
                'A': wb.add_format({'bg_color': '#58FA82',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'C': wb.add_format({'bg_color': '#0000FF',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'G': wb.add_format({'bg_color': '#F7FE2E',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'T': wb.add_format({'bg_color': '#FF0000',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'N': wb.add_format({'bg_color': '#E2CFDD',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size})
            }
            # Create a format to handle IUPAC degenerate bases
            ambiguous_format = wb.add_format({'font_color': '#C70039',
                                              'bg_color': '#E2CFDD',
                                              'font_name': 'Courier New',
                                              'font_size': font_size,
                                              'bold': True})
        else:
            # Create a dictionary to store formats for the difference nucleotides
            format_dict = {
                'A': wb.add_format({'bg_color': '#C8C8C8',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'R': wb.add_format({'bg_color': '#145AFF',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'N': wb.add_format({'bg_color': '#00DCDC',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'D': wb.add_format({'bg_color': '#E60A0A',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'C': wb.add_format({'bg_color': '#E6E600',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'E': wb.add_format({'bg_color': '#E60A0A',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'Q': wb.add_format({'bg_color': '#00DCDC',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'G': wb.add_format({'bg_color': '#EBEBEB',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'H': wb.add_format({'bg_color': '#8282D2',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'I': wb.add_format({'bg_color': '#0F820F',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'L': wb.add_format({'bg_color': '#0F820F',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'K': wb.add_format({'bg_color': '#145AFF',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'M': wb.add_format({'bg_color': '#E6E600',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'F': wb.add_format({'bg_color': '#3232AA',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'P': wb.add_format({'bg_color': '#DC9682',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'S': wb.add_format({'bg_color': '#FA9600',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'T': wb.add_format({'bg_color': '#FA9600',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'W': wb.add_format({'bg_color': '#B45AB4',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'Y': wb.add_format({'bg_color': '#3232AA',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
                'V': wb.add_format({'bg_color': '#0F820F',
                                    'bold': True,
                                    'font_name': 'Courier New',
                                    'font_size': font_size}),
            }
            # Create a format to handle non-standard bases
            ambiguous_format = wb.add_format({'font_color': '#C70039',
                                              'bg_color': '#BEA06E',
                                              'font_name': 'Courier New',
                                              'font_size': font_size,
                                              'bold': True})
        return header, courier, bold_courier, top_bold_courier, annotation, format_dict, ambiguous_format

    @staticmethod
    def format_sequence(pos, sequence_dict, ref_dict, bold_courier, format_dict, ambiguous_format, molecule, ref=False):
        """
        Add the necessary formatting to the current position
        :param pos: type INT: Current position in the sequence of the reference chromosome
        :param sequence_dict: type DICT: Dictionary of strain-specific SNP data. Accessed with pos
        :param ref_dict: type DICT: Dictionary of the reference genome-specific location data. Accessed with pos
        :param bold_courier: XLSX writer formatting block for bold courier
        :param format_dict: type DICT: Dictionary of nucleotide/amino acid-specific XLSX writer formats
        :param ambiguous_format: XLSX writer formatting block for sequences that do not exist in format_dict
        :param molecule: type STR: Creating nucleotide or amino acid outputs
        :param ref: type BOOL: Create outputs for a reference genome rather than a query strain
        :return: sequence: The nucleotide base or amino acid residue at the current position of the strain under
        consideration
        return: base_format: The location-specific formatting to apply to the cell containing the sequence
        """
        # Extract the position-specific sequence from the dictionary
        if molecule == 'nt':
            try:
                sequence = sequence_dict[pos]
            # A missing key indicates that the strain-specific sequence matches the reference sequence
            except KeyError:
                sequence = ref_dict[pos]
        else:
            # The reference dictionary does not include the 'strain_name' key, and missing keys are due to the fact
            # that the SNP position falls in a non-coding region
            if ref:
                try:
                    sequence = sequence_dict[pos]['ref_aa_seq_cds']
                except KeyError:
                    sequence = 'NC'
            # Extract the strain-specific amino acid sequence at the SNP position
            else:
                try:
                    sequence = sequence_dict[pos]['snp_aa_seq_cds']
                # If the key is missing from the dictionary it is either because of a deletion or that the SNP position
                # falls in a non-coding region
                except KeyError:
                    sequence = '-'
        # Determine the format to use for the cell based on the sequence
        # If the sequence matches the reference sequence, it uses the standard black text on
        # white background format
        base_format = bold_courier
        if molecule == 'nt':
            if sequence != ref_dict[pos]:
                # If the sequence is one of A, C, G, T, or N, extract the appropriate format
                # from the dictionary
                try:
                    base_format = format_dict[sequence]
                # If the sequence is a degenerate base (e.g. M), or missing (-), use the
                # 'ambiguous' format
                except KeyError:
                    base_format = ambiguous_format
        else:
            try:
                base_format = format_dict[sequence]
            # If the sequence is a degenerate base (e.g. M), or missing (-), use the
            # 'ambiguous' format
            except KeyError:
                base_format = ambiguous_format
        return sequence, base_format
