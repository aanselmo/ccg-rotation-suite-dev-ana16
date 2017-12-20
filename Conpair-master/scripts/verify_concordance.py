#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.15
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

import sys
import os
import optparse
import math
from collections import defaultdict
from ContaminationModel import *
from ContaminationMarker import *

# modified ConPair script adapted to use as a module
# Program to verify tumor-normal sample concordance prog version 0.15 3/August/2016'
# MMQ = opts.min_mapping_quality
# MBQ = opts.min_base_quality
# AA_BB_only = opts.normal_homozygous_markers_only
def verify_concordance(CONPAIR_DIR, tumor_pileup, normal_pileup,
                       outfile, MARKER_FILE=None, COVERAGE_THRESHOLD=20,
                       MMQ=20, MBQ=20, AA_BB_only=False):

    if MARKER_FILE==None:
        MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers',
                                   'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

    if tumor_pileup == None or normal_pileup == None:
        print('Error: Requires tumor and normal sample pileups')
        sys.exit(1)

    if not os.path.exists(tumor_pileup):
        print('ERROR: Input tumor file {0} cannot be found.'.format(tumor_pileup))
        sys.exit(1)

    if not os.path.exists(normal_pileup):
        print('ERROR: Input normal file {0} cannot be found.'.format(normal_pileup))
        sys.exit(1)


    if not os.path.exists(MARKER_FILE):
        print('ERROR: Marker file {0} cannot be find.'.format(MARKER_FILE))
        sys.exit(2)

    Markers = get_markers(MARKER_FILE)

    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, normal_pileup, min_map_quality=MMQ,
                                                                   min_base_quality=MBQ)
    Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, tumor_pileup, min_map_quality=MMQ,
                                                                  min_base_quality=MBQ)

    concordant = 0
    discordant = 0
    for m in Markers:
        NL = Normal_genotype_likelihoods[m]
        TL = Tumor_genotype_likelihoods[m]
        if NL is None or TL is None:
            continue
        if NL['coverage'] < COVERAGE_THRESHOLD or TL['coverage'] < COVERAGE_THRESHOLD:
            continue
        if AA_BB_only:
            if NL['likelihoods'].index(max(NL['likelihoods'])) == 1:
                continue
        if NL['likelihoods'].index(max(NL['likelihoods'])) == TL['likelihoods'].index(max(TL['likelihoods'])):
            concordant += 1
        else:
            discordant += 1

    if concordant + discordant == 0:
        print(
        'WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(
            COVERAGE_THRESHOLD))
        sys.exit(0)

    if outfile == "-":
        print round(float(concordant) / (concordant + discordant), 3)
        print "Based on " + str(concordant + discordant) + "/" + str(
            len(Markers)) + " markers (coverage per marker threshold: " + str(COVERAGE_THRESHOLD) + " reads)"
        print "Minimum mappinq quality: " + str(MMQ)
        print "Minimum base quality: " + str(MBQ)
    else:
        outfile = open(outfile, 'w')
        outfile.write("Concordance: " + str(round(100.0 * float(concordant) / (concordant + discordant), 2)) + "%\n")
        outfile.write("Based on " + str(concordant + discordant) + "/" + str(
            len(Markers)) + " markers (coverage per marker threshold : " + str(COVERAGE_THRESHOLD) + " reads)\n")
        outfile.write("Minimum mappinq quality: " + str(MMQ) + "\n")
        outfile.write("Minimum base quality: " + str(MBQ) + "\n")
        outfile.close()



