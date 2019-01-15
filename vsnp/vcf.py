#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import SetupLogging
from vsnp.methods import Methods
import logging
import os
__author__ = 'adamkoziol'


class VCF(object):

    def main(self):
        pass

    def __init__(self, path, debug=None):
        """

        :param path:
        """
        SetupLogging(debug=debug)
        if path.startswith('~'):
            self.path = os.path.abspath(os.path.expanduser(os.path.join(path)))
        else:
            self.path = os.path.abspath(os.path.join(path))
        assert os.path.isdir(self.path), 'Invalid path specified: {path}'.format(path=self.path)
        logging.debug('Path is: {path}'.format(path=self.path))

