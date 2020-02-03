#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)


import os
import sys
import optparse
from shutil import move
from tempfile import NamedTemporaryFile

# modified ConPair script adapted to use as a module
# """Program to run GATK Pileup on a single sample"""

def run_gatk_pileup_for_sample(CONPAIR_DIR, GATK, BAM_FILE, OUTFILE, JAVAPATH,
                               MARKER_FILE=None, REFERENCE=None, JAVA_TEMP=None,
                               xmx_java=None):

    if not os.path.exists(BAM_FILE):
        print(('ERROR: Specified bamfile {0} cannot be found.'.format(BAM_FILE)))
        sys.exit(1)

    if not os.path.exists(GATK):
        print(('ERROR: GATK jar {0} cannot be find.'.format(GATK)))
        sys.exit(2)

    if MARKER_FILE == None:
        MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers',
            'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.bed')

    if not os.path.exists(MARKER_FILE):
        print(('ERROR: Marker file {0} cannot be find.'.format(MARKER_FILE)))
        sys.exit(2)

    if REFERENCE == None:
        REFERENCE = os.path.join(CONPAIR_DIR, 'data', 'genomes', 'human_g1k_v37.fa')

    if not os.path.exists(REFERENCE):
        print(('ERROR: Reference genome {0} cannot be find.'.format(REFERENCE)))
        sys.exit(3)

    if JAVA_TEMP != None:
        JAVA_TEMP = "-Djava.io.tmpdir={}".format(JAVA_TEMP)
        if not os.path.isdir(JAVA_TEMP):
            os.makedirs(JAVA_TEMP, 0o770)
    else:
        JAVA_TEMP = ""

    if xmx_java == None:
        xmx_java = '12g'

    command_line = ("{7} {0} -Xmx{1} -jar {2} -T Pileup -R {3} -I {4} -L {5} -o {6} " +
                    "-verbose -rf DuplicateRead --filter_reads_with_N_cigar " +
                    "--filter_mismatching_base_and_quals").format(JAVA_TEMP, xmx_java, GATK, REFERENCE, BAM_FILE,
                                                                  MARKER_FILE, OUTFILE, JAVAPATH)
    print("Running Pileup with the following command:")
    print(command_line)

    os.system(command_line)

'''
    # this is an option form the original script that I haven't used. 
    # instead of dealing with it I am just going to comment it out here.
    # for the current workflow we do not need the functionality of removing
    # chromosome prefixes.
    
    if opts.remove_chr_prefix:
        print("Removing 'chr' prefix...")

        with NamedTemporaryFile(delete=False) as tmp_source:
            with open(opts.outfile) as source_file:
                for line in source_file:
                    if line.startswith("chr"):
                        tmp_source.write(line[3:])
                    else:
                        tmp_source.write(line)

        move(tmp_source.name, source_file.name)
'''
#EOF