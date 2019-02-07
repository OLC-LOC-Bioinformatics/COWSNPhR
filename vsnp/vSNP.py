#!/usr/bin/env python3
from accessoryFunctions.accessoryFunctions import modify_usage_error, SetupLogging
from vsnp.vsnp_vcf_run import VCF
import multiprocessing
from time import time
import click
import sys

__author__ = 'adamkoziol'
__version__ = '0.0.1'

SetupLogging()
start = time()
click_options = [
    click.version_option(version=__version__),
    click.option('-p', '--path',
                 required=True,
                 type=str,
                 help='Specify path of folder containing files to be processed'),
    click.option('-s', '--species',
                 help='OPTIONAL: Used to FORCE species type: af, h37, ab1, ab3, suis1, suis2, suis3, mel1, mel1b, '
                      'mel2, mel3, canis, ceti1, ceti2, ovis, neo, para, typhimurium-atcc13311, typhimurium-14028S, '
                      'typhimurium-LT2, heidelberg-SL476, te_atcc35865, te_09-0932, te_89-0490, te_92-0972, '
                      'te_98-0554, te_mce9, flu, newcaste, belize'),
    click.option('-d', '--debug',
                 is_flag=True,
                 help='Print debugging statements to terminal'),
    click.option('-g', '--get',
                 is_flag=True,
                 help='Debug core functions'),
    click.option('-n', '--no_annotation',
                 is_flag=True,
                 help='Run pipeline without annotation'),
    click.option('-e', '--elite',
                 is_flag=True,
                 help='create a tree with elite sample representation'),
    click.option('-f', '--filter',
                 is_flag=True,
                 help='Find possible positions to filter'),
    click.option('-t', '--threads',
                 type=int,
                 default=multiprocessing.cpu_count() - 1,
                 help='Maximum number of threads to create. Default is number of CPUs - 1'),
    click.option('-q', '--quiet',
                 is_flag=True,
                 help='[**APHIS only**] prevent stats going to cumulative collection'),
    click.option('-m', '--email',
                 type=str,
                 help='[**APHIS only**, Specify your SMTP address for functionality. email options: '
                      'all, s, tod, jess, suelee, chris, email_address'),
    click.option('-u', '--upload',
                 type=str,
                 help='[**APHIS only**, specify own storage for functionality] upload files to the bioinfo drive')
]

# Subcommand-specific options
click_vcf_options = []

click_tree_options = [
    click.option('-a', '--all_vcf',
                 is_flag=True,
                 help='Make tree using all VCFs'),
]


def add_options(options):
    def _add_options(func):
        for option in reversed(options):
            func = option(func)
        return func

    return _add_options


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
def group():
    pass


@group.command()
@add_options(click_options)
@add_options(click_vcf_options)
@add_options(click_tree_options)
def vsnp(**kwargs):
    """
    Full vSNP pipeline
    """
    vsnp_object = VCF(path=kwargs['path'],
                      threads=kwargs['threads'],
                      debug=kwargs['debug'])
    vsnp_object.main()


@group.command()
@add_options(click_options)
@add_options(click_vcf_options)
def vcf(**kwargs):
    """
    VCF creation component of vSNP pipeline
    """
    vcf_object = VCF(path=kwargs['path'],
                     threads=kwargs['threads'],
                     debug=kwargs['debug'])
    vcf_object.main()


@group.command()
@add_options(click_options)
@add_options(click_tree_options)
def tree(**kwargs):
    """
    Phylogenetic tree creation
    """
    tree_object = ''


# Define the list of acceptable sub-programs
program_list = ['vsnp', 'vcf', 'tree']
# Extract the BLAST command to use from the command line arguments
try:
    program = sys.argv[1].lower() if sys.argv[1] in program_list else str()
except IndexError:
    program = str()

# If the program was not specified, call the 'group help'
if not program:
    # Call the help
    group(['--help'])

# Convert the program string to the appropriate subcommand to use when modifying the usage error - ResFinder
# and VirulenceFinder analyses are only available for BLAST programs that use a nt database
subcommand_dict = {
    'vsnp': vsnp,
    'vcf': vcf,
    'tree': tree,
}
try:
    sub_command = subcommand_dict[program]
    # Change the behaviour of click to print the help menu when a subcommand is specified, but is missing arguments
    modify_usage_error(subcommand=sub_command,
                       program_list=program_list)
except KeyError:
    pass

if __name__ == '__main__':
    group()
