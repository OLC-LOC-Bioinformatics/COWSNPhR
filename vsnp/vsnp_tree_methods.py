#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import filer, make_path, relative_symlink, run_subprocess, write_to_logfile
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
import multiprocessing
from glob import glob
import xlsxwriter
import shutil
import pandas
import gzip
import xlrd
import vcf
import os

__author__ = 'adamkoziol'


class VSNPTreeMethods(object):

    @staticmethod
    def file_list(file_path):
        """
        Create a list of all gVCF files present in the supplied path. Accepts .gvcf.gz extension only
        :param file_path: type STR: absolute path of folder containing VCF files
        :return vcf_files: sorted list of all VCF files present in the supplied path
        """
        # Use glob to find the acceptable extensions of VCF files in the supplied path
        vcf_files = glob(os.path.join(file_path, '*.gvcf.gz'))
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
    def load_vcf_library(strain_vcf_dict):
        """
        Use PyVCF to load the gVCF files into dictionaries. This seems very slow, so I am not using this now. Maybe I'm
        doing something wrong?
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to gVCF file
        :return: parsed_vcf_dict: Dictionary of strain name: key: value pairs CHROM': ref_genome, 'REF': ref base,
            'ALT': alt base, 'QUAL': quality score, 'LENGTH': length of feature, 'FILTER': deepvariant filter call,
            'STATS': dictionary of format data
        :return: strain_best_ref_dict: Dictionary of strain name: reference genome parsed from gVCF file. Note that
            this will select only a single 'best reference genome' even if there are multiple contigs in the file
            against which this strain was reference mapped
        :return: strain_best_ref_set_dict: Dictionary of strain name: all reference genomes parsed from gVCF file
        """
        # Initialise dictionaries to store the parsed gVCF outputs and the closest reference genome
        parsed_vcf_dict = dict()
        strain_best_ref_dict = dict()
        strain_best_ref_set_dict = dict()
        for strain_name, vcf_file in strain_vcf_dict.items():
            parsed_vcf_dict[strain_name] = dict()
            # Create the PyVCF 'Reader' object from the .gvcf file
            vcf_reader = vcf.Reader(filename=vcf_file)
            # Iterate through all the records in the .vcf file
            for record in vcf_reader:
                # print(record.POS)
                # Initialise the dictionary as required
                if strain_name not in strain_best_ref_dict:
                    strain_best_ref_dict[strain_name] = record.CHROM
                    strain_best_ref_set_dict[strain_name] = {record.CHROM}
                else:
                    strain_best_ref_set_dict[strain_name].add(record.CHROM)
                # Initialise the length of the list to 1
                alt_length = 1
                # Substitute out <*> for the reference call in the record.ALT list
                record.ALT = [str(entry).replace('<*>', record.REF) for entry in record.ALT]
                # Check if the length of the list is greater than 1 i.e. a SNP call
                if len(record.ALT) > 1:
                    for entry in record.ALT:
                        if len(entry) > 1:
                            alt_length = len(entry)
                # SNPs must have a deepvariant filter of 'PASS', be of length one, and have a quality
                # score above the threshold
                if record.FILTER is not None and record.FILTER != ['RefCall'] and len(record.REF) == 1 \
                        and alt_length == 1 and record.QUAL > 30:
                    record.FILTER = 'PASS'
                    parsed_vcf_dict[strain_name][record.POS] = record
                # Insertions must still have a deepvariant filter of 'PASS', but must have a length
                # greater than one
                elif record.FILTER is not None and record.FILTER != ['RefCall'] and len(record.REF) == 1:
                    record.FILTER = 'INSERTION'
                    parsed_vcf_dict[strain_name][record.POS] = record
                # If the position in the 'info' field does not match pos, and the minimum depth of a gVCF
                # block is 0, this is considered a deletion
                elif record.FILTER is None and record.INFO != record.POS and record.samples[0]['MIN_DP'] == '0':
                    # Iterate through the range of the deletion, and populate the dictionary for each
                    # position encompassed by this range (add +1 due to needing to include the final
                    # position in the dictionary)
                    record.FILTER = 'DELETION'
                    for i in range(record.POS, record.INFO + 1):
                        parsed_vcf_dict[strain_name][i] = record
        return parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

    @staticmethod
    def load_vcf(strain_vcf_dict, threads):
        """
        Create a multiprocessing pool to parse gVCF files concurrently
        :param strain_vcf_dict:
        :param threads:
        :return:
        """
        # Initialise dictionaries to store the parsed gVCF outputs and the closest reference genome
        parsed_vcf_dict = dict()
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
        for parsed_vcf, strain_best_ref, strain_best_ref_set in p.starmap(VSNPTreeMethods.load_vcf_multiprocessing,
                                                                          zip(strain_list,
                                                                              [strain_vcf_dict] * list_length)):
            # Update the dictionaries
            parsed_vcf_dict.update(parsed_vcf)
            strain_best_ref_dict.update(strain_best_ref)
            strain_best_ref_set_dict.update(strain_best_ref_set)
        # Close and join the pool
        p.close()
        p.join()
        return parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

    @staticmethod
    def load_vcf_multiprocessing(strain_name, strain_vcf_dict):
        """
        Load the gVCF files into a dictionary
        :param strain_name: type STR: Name of strain being processed
        :param strain_vcf_dict: type DICT: Dictionary of strain name: absolute path to gVCF file
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
        parsed_vcf_dict = dict()
        strain_best_ref_dict = dict()
        strain_best_ref_set_dict = dict()
        vcf_file = strain_vcf_dict[strain_name]
        parsed_vcf_dict[strain_name] = dict()
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
                        # SNPs must have a deepvariant filter of 'PASS', be of length one, and have a quality
                        # score above the threshold
                        if filter_stat == 'PASS' and len(ref) == 1 and alt_length == 1 and float(qual) > 30:
                            # Populate the dictionary with the required key: value pairs
                            parsed_vcf_dict[strain_name][pos] = {
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
                            parsed_vcf_dict[strain_name][pos] = {
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
                        elif int(info) != pos and format_dict['MIN_DP'] == '0':
                            # Subtract the starting position (pos) from the final position (info)
                            length = int(info) - pos
                            # Iterate through the range of the deletion, and populate the dictionary for each
                            # position encompassed by this range (add +1 due to needing to include the final
                            # position in the dictionary)
                            for i in range(int(pos), int(info) + 1):
                                parsed_vcf_dict[strain_name][i] = {
                                    'CHROM': ref_genome,
                                    'REF': ref,
                                    'ALT': alt,
                                    'QUAL': qual,
                                    'LENGTH': length,
                                    'FILTER': 'DELETION',
                                    'STATS': format_dict
                                }
        return parsed_vcf_dict, strain_best_ref_dict, strain_best_ref_set_dict

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
    def load_snp_positions(strain_parsed_vcf_dict, strain_consolidated_ref_dict):
        """
        Parse the VCF files, and extract all the query and reference genome-specific SNP locations as well as the
        reference sequence
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param strain_consolidated_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :return: ref_snp_positions: Dictionary of reference chromosome name: absolute position: reference base call
        :return: strain_snp_positions: Dictionary of strain name: all strain-specific SNP positions
        """
        # Initialise dictionaries to store the SNP positions
        consolidated_ref_snp_positions = dict()
        strain_snp_positions = dict()
        ref_snp_positions = dict()
        for strain_name, vcf_dict in strain_parsed_vcf_dict.items():
            best_ref = strain_consolidated_ref_dict[strain_name]
            strain_snp_positions[strain_name] = list()
            # Initialise the reference genome key as required
            if best_ref not in consolidated_ref_snp_positions:
                consolidated_ref_snp_positions[best_ref] = dict()
            # Iterate through all the positions
            for pos, pos_dict in vcf_dict.items():
                if pos_dict['CHROM'] not in ref_snp_positions:
                    ref_snp_positions[pos_dict['CHROM']] = dict()
                # Only consider locations that are called 'PASS' in the dictionary
                if pos_dict['FILTER'] == 'PASS':
                    # Populate the dictionary with the position and the reference sequence at that position
                    consolidated_ref_snp_positions[best_ref][pos] = pos_dict['REF']
                    ref_snp_positions[pos_dict['CHROM']][pos] = pos_dict['REF']
                    strain_snp_positions[strain_name].append(pos)
        return consolidated_ref_snp_positions, strain_snp_positions, ref_snp_positions

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
                        # If the position is in the list of strain-specific SNP positions, add the group to the list
                        if position in snp_positions:
                            strain_groups[strain_name].append(group)
        return strain_groups

    @staticmethod
    def determine_group_snp_positions(strain_snp_positions, strain_groups, strain_species_dict):
        """
        Find all the group-specific SNP positions
        :param strain_snp_positions: type DICT: Dictionary of strain name: all strain-specific SNP positions
        :param strain_groups: type DICT:  Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :return: group_positions_set: Dictionary of species code: group name: set of group-specific SNP positions
        """
        # Initialise a dictionary to store all of the species: group-specific SNP positions
        group_positions_set = dict()
        for strain_name, groups in strain_groups.items():
            # Extract the list of positions from the dictionary
            positions = strain_snp_positions[strain_name]
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            # Initialise the species key in the dictionary if required
            if species not in group_positions_set:
                group_positions_set[species] = dict()
            for group in groups:
                # Initialise the group key in the dictionary if required
                if group not in group_positions_set[species]:
                    group_positions_set[species][group] = set()
                # Add the group-specific positions to the set
                for pos in positions:
                    group_positions_set[species][group].add(pos)
        return group_positions_set

    @staticmethod
    def load_snp_sequence(strain_parsed_vcf_dict, strain_consolidated_ref_dict, group_positions_set, strain_groups,
                          strain_species_dict, consolidated_ref_snp_positions):
        """
        Parse the gVCF-derived dictionaries to determine the strain-specific sequence at every SNP position for every
        group
        :param strain_parsed_vcf_dict: type DICT: Dictionary of strain name: dictionary of parsed VCF data
        :param strain_consolidated_ref_dict: type DICT: Dictionary of strain name: extracted reference genome name
        :param group_positions_set: type DICT: Dictionary of species code: group name: set of group-specific SNP
        positions
        :param strain_groups: type DICT: Dictionary of strain name: list of group(s) for which the strain contains the
        defining SNP
        :param strain_species_dict: type DICT: Dictionary of strain name: species code
        :param consolidated_ref_snp_positions: type DICT: Dictionary of reference name: absolute position: reference base call
        :return: group_strain_snp_sequence: Dictionary of strain name: species code: group name: strain name: position:
            strain-specific sequence
        """
        # Dictionary of degenerate IUPAC codes
        iupac = {
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
        # Initialise a dictionary to store all the SNP locations, and the reference sequence for each reference
        # genome
        group_strain_snp_sequence = dict()
        for strain_name, vcf_dict in strain_parsed_vcf_dict.items():
            # Extract the set of groups to which the strain belongs
            groups = strain_groups[strain_name]
            # Extract the species code from the dictionary
            species = strain_species_dict[strain_name]
            best_ref = strain_consolidated_ref_dict[strain_name]
            # Initialise the dictionary with the species key if required
            if species not in group_strain_snp_sequence:
                group_strain_snp_sequence[species] = dict()
            for group in groups:
                # Add the group and strain name keys to the dictionary
                if group not in group_strain_snp_sequence[species]:
                    group_strain_snp_sequence[species][group] = dict()
                if strain_name not in group_strain_snp_sequence[species][group]:
                    group_strain_snp_sequence[species][group][strain_name] = dict()
                # Extract the set of group-specific SNP positions
                position_set = group_positions_set[species][group]
                # The reference sequence only needs to be determined once. This variable will check to see if it
                # already has been added to the dictionary
                write_ref = False
                if best_ref not in group_strain_snp_sequence[species][group]:
                    group_strain_snp_sequence[species][group][best_ref] = dict()
                    write_ref = True
                for pos in sorted(position_set):
                    # Include the reference position if necessary
                    if write_ref:
                        ref_pos = consolidated_ref_snp_positions[best_ref][pos]
                        group_strain_snp_sequence[species][group][best_ref][pos] = ref_pos
                    # gVCF blocks will compress stretches of normal matches. I haven't added these regions to the
                    # dictionary, so these positions will yield KeyErrors
                    try:
                        # Extract the gVCF dictionary from the position-specific dictionary
                        pos_dict = vcf_dict[pos]
                        # Deletions are recorded as a '-'
                        if pos_dict['FILTER'] == 'DELETION':
                            group_strain_snp_sequence[species][group][strain_name][pos] = '-'
                        else:
                            # Determine if the allele frequency is a mixed population
                            allele_freq = float(pos_dict['STATS']['VAF'].split(',')[0])
                            if allele_freq < 0.8:
                                # Determine the IUPAC code for any multi-allelic sites
                                for code, components in iupac.items():
                                    if sorted(pos_dict['ALT']) == sorted(components):
                                        group_strain_snp_sequence[species][group][strain_name][pos] = code
                            # Otherwise, use the first entry in the 'ALT' field (remember that this field is at least
                            # nt long: alt base,ref base
                            else:
                                group_strain_snp_sequence[species][group][strain_name][pos] = pos_dict['ALT'][0]
                    # If the entry isn't in the dictionary, it is because it matches the reference sequence
                    except KeyError:
                        group_strain_snp_sequence[species][group][strain_name][pos] \
                            = consolidated_ref_snp_positions[best_ref][pos]
        return group_strain_snp_sequence

    @staticmethod
    def create_multifasta(group_strain_snp_sequence, fasta_path):
        """
        Create a multiple sequence alignment in FASTA format for each group from all the SNP positions for the group
        :param group_strain_snp_sequence: type DICT: Dictionary of species: group: strain name: position: sequence
        :param fasta_path: type STR: Absolute path of folder in which alignments are to be created
        :return: group_fasta_dict: Dictionary of group name: FASTA file created for the group
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
                output_dir = os.path.join(fasta_path, species, group)
                make_path(output_dir)
                # Add the group-specific folder to the set of all group folders
                group_folders.add(output_dir)
                for strain_name, sequence_dict in strain_dict.items():
                    # Create a string to store the strain-specific sequence
                    strain_group_seq = str()
                    for pos, sequence in sequence_dict.items():
                        # Append each base to the string
                        strain_group_seq += sequence
                    # Create a SeqRecord from the sequence string in IUPAC ambiguous DNA format. Use the strain name
                    # as the id
                    record = SeqRecord(Seq(strain_group_seq, IUPAC.ambiguous_dna),
                                       id=strain_name,
                                       description='')
                    # Set the name of the FASTA alignment file
                    group_fasta = os.path.join(output_dir, '{group}_alignment.fasta'.format(group=group))
                    # Use SeqIO to append the FASTA sequence to the alignment file
                    with open(group_fasta, 'a+') as fasta:
                        SeqIO.write(record, fasta, 'fasta')
                    # Add the alignment file to the set of all alignment files
                    group_fasta_dict[species][group] = group_fasta
        return group_folders, species_folders, group_fasta_dict

    @staticmethod
    def load_genbank_file(reference_link_path_dict, strain_best_ref_set_dict, dependency_path):
        """
        Use SeqIO to parse the best reference genome GenBank file for annotating SNP locations
        :param reference_link_path_dict: type DICT: Dictionary of strain name: relative path to reference genome
        dependency folder
        :param strain_best_ref_set_dict: type DICT: Dictionary of strain name: set of strain-specific reference genomes
        :param dependency_path: type STR: Absolute path to dependencies
        :return: full_best_ref_gbk_dict: Dictionary of best ref: ref position: SeqIO parsed GenBank file-sourced
        records from closest reference genome for that position
        """
        # Initialise a dictionary to store the SeqIO parsed GenBank files
        best_ref_gbk_dict = dict()
        for strain_name, best_ref_path in reference_link_path_dict.items():
            # Extract the species code from the dictionary
            best_ref_set = strain_best_ref_set_dict[strain_name]
            for best_ref in best_ref_set:
                gbk_file = glob(os.path.join(
                    dependency_path, best_ref_path, '{br}*.gbk'.format(br=os.path.splitext(best_ref)[0])))[0]
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
    def annotate_snps(group_strain_snp_sequence, full_best_ref_gbk_dict, strain_best_ref_set_dict, ref_snp_positions):
        """
        Use GenBank records to annotate each SNP with 'gene', 'locus', and 'product' details
        :param group_strain_snp_sequence: type DICT: Dictionary of species: group: strain name: position: sequence
        :param full_best_ref_gbk_dict: type DICT: Dictionary of best ref: ref position: SeqIO parsed GenBank
        file-sourced records from closest reference genome for that position
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
                for strain_name, pos_dict in strain_dict.items():
                    if strain_name in strain_best_ref_set_dict:
                        # Extract the set of all reference chromosomes used in the reference mapping
                        ref_set = strain_best_ref_set_dict[strain_name]
                        for ref in sorted(ref_set):
                            # As the reference only needs to be added to the dictionary once, ensure that it is not
                            # present before continuing
                            if ref not in species_group_annotated_snps_dict[species][group]:
                                species_group_annotated_snps_dict[species][group][ref] = dict()
                                # gbk_records = best_ref_gbk_dict[ref]
                                gbk_pos_dict = full_best_ref_gbk_dict[ref]
                                for pos in pos_dict:
                                    # Ensure that the position is in the specific chromosome being considered e.g.
                                    # 'NC_017250.1' vs 'NC_017251.1'
                                    if pos in ref_snp_positions[ref]:
                                        species_group_annotated_snps_dict[species][group][ref][pos] = dict()
                                        # Non-coding regions will not be present in the dictionary
                                        try:
                                            # Extract the SeqIO-parsed GenBank feature from the GenBank record
                                            # dictionary
                                            feature = gbk_pos_dict[pos]
                                            # Populate the 'locus', 'gene', and 'product' key: value pairs. Extract
                                            # the values from the feature.qualifiers OrderedDict list
                                            species_group_annotated_snps_dict[species][group][ref][pos]['locus'] = \
                                                feature.qualifiers['locus_tag'][0]
                                            # Not all features have the 'gene' key. Add 'None' if this is the case
                                            try:
                                                gene = feature.qualifiers['gene'][0]
                                            except KeyError:
                                                gene = 'None'
                                            species_group_annotated_snps_dict[species][group][ref][pos]['gene'] \
                                                = gene
                                            species_group_annotated_snps_dict[species][group][ref][pos]['product'] \
                                                = \
                                                feature.qualifiers['product'][0]
                                        # Populate negative key: value pairs if the position is not in the dictionary
                                        except KeyError:
                                            species_group_annotated_snps_dict[species][group][ref][pos]['locus'] \
                                                = 'None'
                                            species_group_annotated_snps_dict[species][group][ref][pos]['gene'] \
                                                = 'None'
                                            species_group_annotated_snps_dict[species][group][ref][pos]['product'] \
                                                = 'None'
        return species_group_annotated_snps_dict
