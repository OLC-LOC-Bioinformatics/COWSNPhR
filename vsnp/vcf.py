#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from vsnp.methods import Methods
import logging
import os
__author__ = 'adamkoziol'


class VCF(object):

    def main(self):
        """

        """
        self.fastq_manipulation()
        self.best_reference_calculation()

    def fastq_manipulation(self):
        """

        """
        logging.info('Locating FASTQ files, creating strain-specific working directories and symlinks to files')
        self.fastq_files = Methods.file_list(path=self.path)
        logging.debug('FASTQ files: \n{fastq_files}'.format(fastq_files='\n'.join(self.fastq_files)))
        self.strain_folder_dict = Methods.strain_list(fastq_files=self.fastq_files)
        self.strain_name_dict = Methods.strain_namer(strain_folders=self.strain_folder_dict)
        logging.debug('Strain names: \n{strain_names}'.format(strain_names='\n'.join(sorted(self.strain_name_dict))))
        self.strain_fastq_dict = Methods.file_link(strain_folder_dict=self.strain_folder_dict,
                                                   strain_name_dict=self.strain_name_dict)
        logging.debug(
            'Strain-specific symlinked FASTQ files: \n{symlinks}'.format(
                symlinks='\n'.join(['{strain_name}: {fastq_files}'.format(strain_name=sn, fastq_files=ff)
                                    for sn, ff in self.strain_fastq_dict.items()])))

    def best_reference_calculation(self):
        """

        """
        logging.info('Running MASH analyses')
        self.fastq_sketch_dict = Methods.call_mash_sketch(strain_fastq_dict=self.strain_fastq_dict,
                                                          strain_name_dict=self.strain_name_dict,
                                                          threads=self.threads,
                                                          logfile=self.logfile)
        logging.info('Parsing MASH outputs to determine closest reference genomes')
        self.mash_output_dict = Methods.call_mash_dist(strain_fastq_dict=self.strain_fastq_dict,
                                                       strain_name_dict=self.strain_name_dict,
                                                       fastq_sketch_dict=self.fastq_sketch_dict,
                                                       ref_sketch_file=os.path.join(
                                                           self.dependencypath, 'mash', 'vsnp_reference.msh'),
                                                       threads=self.threads,
                                                       logfile=self.logfile)


    def __init__(self, path, threads, debug=False):
        """

        :param path: type STR: Path of folder containing FASTQ files
        :param threads: type INT: Number of threads to use in the analyses
        :param debug: type BOOL: Boolean of whether debug level logs are printed to terminal
        """
        SetupLogging(debug=debug)
        # Determine the path in which the sequence files are located. Allow for ~ expansion
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        # Ensure that the path exists
        assert os.path.isdir(self.path), 'Invalid path specified: {path}'.format(path=self.path)
        logging.debug('Supplied sequence path: \n{path}'.format(path=self.path))
        # Initialise class variables
        self.threads = threads
        self.reportpath = os.path.join(self.path, 'reports')
        make_path(self.reportpath)
        assert os.path.isdir(self.reportpath), 'Could not create report path as requested: {rp}'\
            .format(rp=self.reportpath)
        # Extract the path of the folder containing this script
        self.scriptpath = os.path.abspath(os.path.dirname(__file__))
        # Use the script path to set the absolute path of the dependencies folder
        self.dependencypath = os.path.join(os.path.dirname(self.scriptpath), 'dependencies')
        assert os.path.isdir(self.dependencypath), 'Something went wrong with the install. Cannot locate the ' \
                                                   'dependencies folder in: {sp}'.format(sp=self.scriptpath)
        self.logfile = os.path.join(self.path, 'log')
        self.fastq_files = list()
        self.strain_folder_dict = dict()
        self.strain_name_dict = dict()
        self.strain_fastq_dict = dict()
        self.fastq_sketch_dict = dict()
        self.mash_output_dict = dict()

