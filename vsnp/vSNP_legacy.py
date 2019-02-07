#!/usr/bin/env python

import os
import sys
import multiprocessing
import argparse
import textwrap
import glob
from concurrent import futures
import re
from itertools import repeat as itertools_repeat

import functions

root_dir = str(os.getcwd())


def fix_vcf(each_vcf, arg_options):
    mal = []
    # Fix common VCF errors
    temp_file = each_vcf + ".temp"
    write_out = open(temp_file, 'w') #r+ used for reading and writing to the same file
    initial_file_time_stats = os.stat(each_vcf)
    with open(each_vcf, 'r') as file:
        try:
            for line in file:
                if line.rstrip(): # true if not empty line'^$'
                    line = line.rstrip() #remove right white space
                    line = re.sub('"AC=', 'AC=', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('""', '"', line)
                    line = re.sub('"$', '', line)
                    line = re.sub('GQ:PL\t"', 'GQ:PL\t', line)
                    line = re.sub('[0-9]+\tGT\t.\/.$', '999\tGT:AD:DP:GQ:PL\t1/1:0,80:80:99:2352,239,0', line)
                    line = re.sub('^"', '', line)
                    if line.startswith('##') and line.endswith('"'):
                        line = re.sub('"$', '', line)
                    if line.startswith('##'):
                        line = line.split('\t')
                        line = ''.join(line[0])
                    if not line.startswith('##'):
                        line = re.sub('"', '', line)
                        line = line.split('\t')
                        line = "\t".join(line[0:10])
                        print(line, file=write_out)
                    else:
                        print(line, file=write_out)
        except IndexError:
            print("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            mal.append("##### IndexError: Deleting corrupt VCF file: " + each_vcf)
            os.remove(each_vcf)
        except UnicodeDecodeError:
            print("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            mal.append("##### UnicodeDecodeError: Deleting corrupt VCF file: " + each_vcf)
            os.remove(each_vcf)

    write_out.close()
    os.rename(temp_file, each_vcf)
    # revert timestamp to original allows elites to properly sort on file modification time
    os.utime(each_vcf, times=(initial_file_time_stats.st_mtime, initial_file_time_stats.st_mtime))
    return mal

parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

---------------------------------------------------------
vSNP --> get SNPs, group SNPs, verify SNPs

vSNP is called on a working directory containing FASTQ or VCF files.

See documentation at: https://usda-vs.github.io/snp_analysis/

        Step 1: FASTQs --> VCF

        Step 2: VCFs --> Tables & Trees

-s <OPTIONAL SPECIES TYPES>: af, h37, ab1, ab3, suis1, suis2, suis3, mel1, mel1b, mel2, mel3, canis, ceti1, ceti2, ovis, neo, para, typhimurium-atcc13311, typhimurium-14028S, typhimurium-LT2, heidelberg-SL476, te_atcc35865, te_09-0932, te_89-0490, te_92-0972, te_98-0554, te_mce9, flu, newcaste, belize

'''), epilog='''---------------------------------------------------------''')

#universal
parser.add_argument('-s', '--species', action='store', dest='species', help='OPTIONAL: Used to FORCE species type <see options above>')
parser.add_argument('-d', '--debug', action='store_true', dest='debug_call', help='debug, run without pool.map')
parser.add_argument('-g', '--get', action='store_true', dest='get', help='get, get to the core functions for debugging')
parser.add_argument('-n', '--no_annotation', action='store_true', dest='no_annotation', help='no_annotation, run without annotation')
parser.add_argument('-a', '--all_vcf', action='store_true', dest='all_vcf', help='make tree using all VCFs')
parser.add_argument('-o', '--only_all_vcf', action='store_true', dest='only_all_vcf', help='make tree using all VCFs')
parser.add_argument('-e', '--elite', action='store_true', dest='elite', help='create a tree with on elite sample representation')
parser.add_argument('-f', '--filter', action='store_true', dest='filter_finder', help='Find possible positions to filter')
parser.add_argument('-p', '--processor', action='store', dest='processor', help='max processor usage')
parser.add_argument('-q', '--quiet', action='store_true', dest='quiet', help='[**APHIS only**] prevent stats going to cumlative collection')
parser.add_argument('-m', '--email', action='store', dest='email', help='[**APHIS only**, specify own SMTP address for functionality] email options: all, s, tod, jess, suelee, chris, email_address')
parser.add_argument('-u', '--upload', action='store_true', dest='upload', help='[**APHIS only**, specify own storage for functionality] upload files to the bioinfo drive')
args = parser.parse_args()
if args.only_all_vcf:
    args.all_vcf = True
print("\nSET ARGUMENTS: ")
print(args)
arg_options = {
    "species": args.species,
    "debug_call": args.debug_call,
    "get": args.get,
    "no_annotation": args.no_annotation,
    "all_vcf": args.all_vcf,
    "only_all_vcf": args.only_all_vcf,
    "elite": args.elite,
    "filter_finder": args.filter_finder,
    "processor": args.processor,
    "quiet": args.quiet,
    "upload": args.upload,
}
print("")

email_list = None
if args.email == "all":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, patrick.m.camp@aphis.usda.gov, David.T.Farrell@aphis.usda.gov, Robin.L.Swanson@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov"
elif args.email == "tod":
    email_list = "tod.p.stuber@aphis.usda.gov"
elif args.email == "jess":
    email_list = "Jessica.A.Hicks@aphis.usda.gov"
elif args.email == "suelee":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, Doris.M.Bravo@aphis.usda.gov, kristina.lantz@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
elif args.email == "suelee-":
    email_list = "tod.p.stuber@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov"
elif args.email == "tyler-":
    email_list = "tod.p.stuber@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov"
elif args.email == "chris":
    email_list = "tod.p.stuber@aphis.usda.gov, Jessica.A.Hicks@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov, Suelee.Robbe-Austerman@aphis.usda.gov, eto3@cdc.gov, kristina.lantz@aphis.usda.gov, Tyler.Thacker@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
elif args.email == "chris-":
    email_list = "tod.p.stuber@aphis.usda.gov, Christine.R.Quance@aphis.usda.gov"
elif args.email == "kris":
    email_list = "kristina.lantz@aphis.usda.gov, tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, mary.k.smith@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
elif args.email == "doris":
    email_list = "tod.p.stuber@aphis.usda.gov, jessica.a.hicks@aphis.usda.gov, doris.m.bravo@aphis.usda.gov, suelee.robbe-austerman@aphis.usda.gov, kristina.lantz@aphis.usda.gov, patrick.m.camp@aphis.usda.gov"
else:
    email_list = None

arg_options['email_list'] = email_list
################################################################################################################################################

all_file_types_count = len(glob.glob('*.*'))
fastq_check = len(glob.glob('*fastq.gz'))
vcf_check = len(glob.glob('*vcf'))

# Check that there are either FASTQs or VCFs, not both
if fastq_check > 0:
    fastq_check = True
if vcf_check > 0:
    vcf_check = True
if fastq_check and vcf_check:
    print("\n#####You have a mix of FASTQ and VCF files.  This is not allowed\n\n")
    sys.exit(0)

arg_options['root_dir'] = root_dir

if arg_options['processor']:
    cpu_count = int(args.processor)
else:
    cpu_count = int(multiprocessing.cpu_count())
limited_cpu_count = int(cpu_count / 4)
if limited_cpu_count == 0:
    limited_cpu_count = 1

arg_options['cpu_count'] = cpu_count
arg_options['limited_cpu_count'] = limited_cpu_count

# Check that there an equal number of both R1 and R2 reads
if fastq_check:
    # Pair check
    pair_check = len(glob.glob('*_R2*fastq.gz'))
    if pair_check > 0:
        R1 = glob.glob('*_R1*fastq.gz')
        R2 = glob.glob('*_R2*fastq.gz')
        R1count = len(R1)
        R2count = len(R2)
    else:
        R1 = glob.glob('*fastq.gz')
        R2 = None

    # fastq_count = R1count + R2count
    # if (fastq_count % 2 != 0):
    #     print("\n#####Check paired files.  Unpaired files seen by odd number of counted FASTQs\n\n")
    #     sys.exit(0)
    # if (R1count != R2count):
    #     print("\n#####Check paired files.  R1 files do not equal R2\n\n")
    #     sys.exit(0)
    # if (all_file_types_count != fastq_count):
    #     print("\n#####Only zipped FASTQ files are allowed in directory\n\n")
    #     sys.exit(0)
    # elif (fastq_count > 1):
    if arg_options['all_vcf'] or arg_options['elite'] or arg_options['upload'] or arg_options['filter_finder']:
        print("#####Incorrect use of options when running loop/script 1")
        sys.exit(0)

    print("\n--> RUNNING LOOP/SCRIPT 1\n")
    #Enter script 1 -->
    functions.run_loop(arg_options)
    print("See files, vSNP has finished alignments")
elif vcf_check:
    #fix files
    malformed = []
    vcf_list = glob.glob('*vcf')
    print("Fixing files...\n")
    if arg_options['debug_call'] and not arg_options['get']:
        for each_vcf in vcf_list:
            print(each_vcf)
            mal = fix_vcf(each_vcf, arg_options)
            malformed = list(mal)
    else:
        with futures.ProcessPoolExecutor() as pool:
            mal = pool.map(fix_vcf, vcf_list, itertools_repeat(arg_options))
            malformed = malformed + list(mal)
    malformed = [x for x in malformed if x] # remove blanks
    print("done fixing")
    print(malformed)
    arg_options['malformed'] = malformed

    if not arg_options['species']:
        species = functions.get_species(arg_options)
        if species is None:
            print("\nEXITED\n##### Unable to find a species corresponding to CHROM found in VCF files")
            sys.exit(0)
        arg_options['species'] = species
        print("species %s" % species)
    vcfs_count = len(glob.glob('*vcf'))
    if (all_file_types_count != vcfs_count):
        print("\n##### You have more than just VCF files in your directory.  Only VCF files are allowed if running script 2\n\n")
        sys.exit(0)
    else:
        if arg_options['quiet']:
            print("#####Incorrect use of options when running script 2, when running step 2 -q cannot be used")
            sys.exit(0)
        else:
            if arg_options['species']:
                print("\n--> RUNNING SCRIPT 2\n")
                #Enter script 2 -->
                functions.run_script2(arg_options)
            else:
                print("#####Based on VCF CHROM id (reference used to build VCF) a matching species cannot be found neither was there a -s option given")
                sys.exit(0)
else:
    print("#####Error determining file type.")
    sys.exit(0)
