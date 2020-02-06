#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)


import sys
import os
import imp
from collections import defaultdict
import numpy as np
from math import pow

# modified ConPair script adapted to use as a module
# """Program to estimate tumor-normal sample contamination"""
# grid_precision = opts.grid
# MMQ = opts.min_mapping_quality
def estimate_tumor_normal_contamination(CONPAIR_DIR, tumor_pileup, normal_pileup, outfile,
                                        MARKER_FILE=None, grid_precision=0.01, MMQ=20):
    HOMOZYGOUS_P_VALUE_THRESHOLD = 0.999

    ContaminationModel = imp.load_source('/ContaminationModel', CONPAIR_DIR + '/modules/ContaminationModel.py')
    ContaminationMarker = imp.load_source('/ContaminationMarker', CONPAIR_DIR + '/modules/ContaminationMarker.py')
    MathOperations = imp.load_source('/MathOperations', CONPAIR_DIR + '/modules/MathOperations.py')
    Genotypes = imp.load_source('/Genotypes', CONPAIR_DIR + '/modules/Genotypes.py')

    if not tumor_pileup or not normal_pileup:
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(tumor_pileup):
        print(('ERROR: Input tumor file {0} cannot be find.'.format(tumor_pileup)))
        sys.exit(1)

    if not os.path.exists(normal_pileup):
        print(('ERROR: Input normal file {0} cannot be find.'.format(normal_pileup)))
        sys.exit(1)

    if MARKER_FILE==None:
        MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

    if not os.path.exists(MARKER_FILE):
        print(('ERROR: Marker file {0} cannot be find.'.format(MARKER_FILE)))
        sys.exit(2)

    def drange(start, stop, step):
        r = start
        while r < stop:
            yield r
            r += step

    Markers = ContaminationMarker.get_markers(MARKER_FILE)

    Normal_homozygous_genotype = defaultdict(lambda: defaultdict())

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    checkpoints.append(0.5)
    Scores = ContaminationModel.create_conditional_likelihood_of_base_dict(checkpoints)
    checkpoints = [i for i in drange(0.0, 0.5, grid_precision)]
    checkpoints.append(0.5)

    if outfile != "-":
        outfile = open(outfile, 'w')

    ### PARSING THE NORMAL PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(normal_pileup)
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)
        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue

        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        AA_likelihood, AB_likelihood, BB_likelihood = Genotypes.compute_genotype_likelihood(pileup.Quals[marker.ref], pileup.Quals[marker.alt], normalize=True)

        if AA_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.ref, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}
        elif BB_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.alt, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}


        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)


        priors = [lPAA*2, lPAA+lPBB, lPAA+lPAB, lPAB*2, lPAB+lPAA, lPAB+lPBB, lPBB*2, lPBB+lPAA, lPBB+lPAB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

    file.close()
    D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
    ARGMAX = np.argmax(D)
    cont = checkpoints[ARGMAX]

    x1 = max(cont-grid_precision, 0.0)
    x2 = cont
    x3 = min(cont+grid_precision, 1.0)

    if x2 == 0.0:
        x2 += grid_precision/100
    elif x2 == 1.0:
        x2 -= grid_precision/100

    ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

    optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, (x1, x2, x3))

    ### PRINTING THE NORMAL RESULTS

    if outfile == "-":
        print(("Normal sample contamination level: " + str(round(100.0*optimal_val, 3)) + "%"))
    else:
        outfile.write("Normal sample contamination level: " + str(round(100.0*optimal_val, 3)) + "%\n")


    ### PARSING THE TUMOR PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(tumor_pileup)

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)

        try:
            normal_hom_genotype = Normal_homozygous_genotype[pileup.chrom][pileup.pos]['genotype']
        except:
            continue

        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue

        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        Normal_info = Normal_homozygous_genotype[pileup.chrom][pileup.pos]
        AA_likelihood = Normal_info['AA_likelihood']
        AB_likelihood = Normal_info['AB_likelihood']
        BB_likelihood = Normal_info['BB_likelihood']
        nlPAA = MathOperations.log10p(AA_likelihood)
        nlPAB = MathOperations.log10p(AB_likelihood)
        nlPBB = MathOperations.log10p(BB_likelihood)


        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)
        priors = [lPAA+nlPAA, lPBB+nlPAA,lPAB+nlPAA, lPAB+nlPAB, lPAA+nlPAB, lPBB+nlPAB, lPBB+nlPBB, lPAA+nlPBB, lPAB+nlPBB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

    file.close()

    D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
    ARGMAX = np.argmax(D)
    cont = checkpoints[ARGMAX]

    x1 = max(cont-grid_precision, 0.0)
    x2 = cont
    x3 = min(cont+grid_precision, 1.0)

    if x2 == 0.0:
        x2 += grid_precision/100
    elif x2 == 1.0:
        x2 -= grid_precision/100

    ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

    optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, (x1, x2, x3))

    ### PRINTING THE TUMOR RESULTS

    if outfile == "-":
        print(("Tumor sample contamination level: " + str(round(100.0*optimal_val,3)) + "%"))
    else:
        outfile.write("Tumor sample contamination level: " + str(round(100.0*optimal_val,3)) + "%\n")
        outfile.close()
