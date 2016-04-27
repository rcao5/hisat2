#!/usr/bin/env python

#
# Copyright 2016, Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#


import sys, os, subprocess, re
import inspect
import random
import glob
from argparse import ArgumentParser, FileType

"""
"""
def test_HLA_genotyping(reference_type,
                        hla_list,
                        partial,
                        aligners,
                        exclude_allele_list,
                        num_mismatch,
                        sample,
                        threads,
                        verbose):
    # Current script directory
    curr_script = os.path.realpath(inspect.getsourcefile(test_HLA_genotyping))
    ex_path = os.path.dirname(curr_script)

    read_dir = "CAAPA"

    if not os.path.exists(read_dir):
        print >> sys.stderr, "Error: CAAPA is needed"
        sys.exit(1)

    done_fnames = glob.glob("%s/*.done" % read_dir)
    for done_fname in done_fnames:
        done_fname = done_fname.split('/')[-1]
        genome = base_fname = done_fname.split('.')[0]
        fq_fname = "%s/%s.extracted.fq.1.gz" % (read_dir, base_fname)
        fq_fname2 = "%s/%s.extracted.fq.2.gz" % (read_dir, base_fname)
        if not os.path.exists(fq_fname) or not os.path.exists(fq_fname2):
            continue

        print >> sys.stderr, genome
        cmd_aligners = ['.'.join(aligners[i]) for i in range(len(aligners))]
        test_hla_script = os.path.join(ex_path, "hisat2_test_HLA_genotyping.py")
        test_hla_cmd = [test_hla_script,
                        "--reference-type", reference_type,
                        "--hla-list", ','.join(hla_list),
                        "--aligner-list", ','.join(cmd_aligners),
                        "--reads", "%s,%s" % (fq_fname, fq_fname2),
                        # "--best-alleles",
                        # "--exclude-allele-list", ','.join(exclude_allele_list),
                        "--sample", str(sample),
                        "--no-pair-model",
                        "-p", str(threads),
                        "--num-mismatch", str(num_mismatch)]
        if partial:
            test_hla_cmd += ["--partial"]

        if verbose:
            print >> sys.stderr, ' '.join(test_hla_cmd)
            
        proc = subprocess.Popen(test_hla_cmd, stderr=subprocess.PIPE, stdout=open("/dev/null", 'w'))
        for line in proc.stderr:
            line = line.strip()
            if line == "" or line.find("ranked") == -1:
                continue
            ranking, _, allele, abundance = line.split()
            abundance = abundance.split(':')[1]
            abundance = abundance.split('%')[0]
            print >> sys.stdout, genome, allele, abundance
            print >> sys.stderr, "\t", genome, allele, abundance
        proc.communicate()


"""
"""
if __name__ == '__main__':
    parser = ArgumentParser(
        description='test HLA genotyping for CAAPA Genomes')
    parser.add_argument("--reference-type",
                        dest="reference_type",
                        type=str,
                        default="gene",
                        help="Reference type: gene, chromosome, and genome (default: gene)")
    parser.add_argument("--hla-list",
                        dest="hla_list",
                        type=str,
                        default="A,B,C,DQA1,DQB1,DRB1",
                        help="A comma-separated list of HLA genes (default: A,B,C,DQA1,DQB1,DRB1)")
    parser.add_argument('--partial',
                        dest='partial',
                        action='store_true',
                        help='Include partial alleles (e.g. A_nuc.fasta)')
    parser.add_argument("--aligner-list",
                        dest="aligners",
                        type=str,
                        default="hisat2.graph",
                        help="A comma-separated list of aligners (default: hisat2.graph)")
    parser.add_argument("--exclude-allele-list",
                        dest="exclude_allele_list",
                        type=str,
                        default="",
                        help="A comma-separated list of allleles to be excluded")
    parser.add_argument("--num-mismatch",
                        dest="num_mismatch",
                        type=int,
                        default=0,
                        help="Maximum number of mismatches per read alignment to be considered (default: 0)")
    parser.add_argument("--sample",
                        dest="sample",
                        type=float,
                        default=1.01,
                        help="Reads sample rate (default: 1.01)")
    parser.add_argument("-p", "--threads",
                        dest="threads",
                        type=int,
                        default=1,
                        help="Number of threads")
    parser.add_argument('-v', '--verbose',
                        dest='verbose',
                        action='store_true',
                        help='also print some statistics to stderr')

    args = parser.parse_args()

    if not args.reference_type in ["gene", "chromosome", "genome"]:
        print >> sys.stderr, "Error: --reference-type (%s) must be one of gene, chromosome, and genome." % (args.reference_type)
        sys.exit(1)
    args.hla_list = args.hla_list.split(',')
    if args.aligners == "":
        print >> sys.stderr, "Error: --aligners must be non-empty."
        sys.exit(1)    
    args.aligners = args.aligners.split(',')
    for i in range(len(args.aligners)):
        args.aligners[i] = args.aligners[i].split('.')
    args.exclude_allele_list = args.exclude_allele_list.split(',')

    test_HLA_genotyping(args.reference_type,
                        args.hla_list,
                        args.partial,
                        args.aligners,
                        args.exclude_allele_list,
                        args.num_mismatch,
                        args.sample,
                        args.threads,
                        args.verbose)
