#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import SetupLogging
from vsnp.vsnp_tree_methods import VSNPTreeMethods
from datetime import datetime
import logging
import os

__author__ = 'adamkoziol'


class VSNPTree(object):

    def main(self):
        """
        Run all the vSNP tree-specific methods
        """
        pass

    def __init__(self, path, threads, debug=False):
        """
        :param path: type STR: Path of folder containing VCF files
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
        self.report_path = os.path.join(self.path, 'reports')
        # Extract the path of the folder containing this script
        self.scriptpath = os.path.abspath(os.path.dirname(__file__))
        # Use the script path to set the absolute path of the dependencies folder
        self.dependencypath = os.path.join(os.path.dirname(self.scriptpath), 'dependencies')
        assert os.path.isdir(self.dependencypath), 'Something went wrong with the install. Cannot locate the ' \
                                                   'dependencies folder in: {sp}'.format(sp=self.scriptpath)
        self.logfile = os.path.join(self.path, 'log')
        self.start_time = datetime.now()
