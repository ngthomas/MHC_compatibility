
#!/usr/bin/env python

"""
a python script that replaces filtered/consensus haplotype back into reference sequence
"""

from __future__ import print_function
import os
import sys
import argparse
import csv
import string
import itertools
from itertools import izip
import collections

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = itertools.tee(iterable)
    next(b, None)
    return izip(a, b)

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-inFasta', help="input reference fasta file", type=argparse.FileType('r'))
    parser.add_argument('-inHaplo', help="input post-processed microhaplotype csv file", type=argparse.FileType('r'))
    parser.add_argument('-inVcf', help="input vcf file", type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="Output full length fasta file",
                        default=sys.stdout, type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    if args.inFasta is None:
        args.inFasta = open('/Users/tng/Projects/MHC_SH/input/Ref_MHC_steelhead_mod.fasta', 'r')
        args.inHaplo = open('/Users/tng/Projects/MHC_SH/input/reported_haplotype-5.csv', 'r')
        args.inVcf = open('/Users/tng/Projects/MHC_SH/input/SH_MHC_gtseq2_cleaned_RD100.recode.vcf', 'r')
        args.outfile = open('/Users/tng/Projects/MHC_SH/output/refined_haplotype.fasta', 'w')

    # assume that the second line of the fasta file contains the reference sequence you need
    refseq = ""
    for fastq in args.inFasta:
        if (not fastq.startswith(">",0,1)):
           refseq = fastq

    #print(refseq[17:26])

    variantPos = [0]

   # grab the entry position
    for vcfLine in args.inVcf:
       if (not vcfLine.startswith("#",0,1)):
          vcfField = vcfLine.split("\t")
          variantPos.append(int(vcfField[1]))

    variantPos.append(len(refseq))

    print(variantPos)


    with args.inHaplo as csvfile:
        haploEntry = csv.reader(csvfile, delimiter=',', doublequote=True)
        for col in haploEntry:
            if (col[1] == "group"):
                continue

            indivId = col[3]
            haplo1 = list(col[4])
            haplo1.append("")
            args.outfile.write(str(">"+indivId+"_1\n"))
            strbits = [refseq[a:b-1] for a,b in pairwise(variantPos)]
            args.outfile.write(string.join([j for i in zip(strbits, haplo1) for j in i], sep=""))
            args.outfile.write("\n")
            #print(strbits)
            haplo2 = list(col[5])
            haplo2.append("")
            args.outfile.write(str(">"+indivId+"_2\n"))
            strbits = [refseq[a:b-1] for a,b in pairwise(variantPos)]
            args.outfile.write(string.join([j for i in zip(strbits, haplo2) for j in i], sep=""))
            args.outfile.write("\n")


            #print(row[1])

    #print(args)

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
