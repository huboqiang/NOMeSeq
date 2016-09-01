#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os

import cPickle as pickle
import numpy as np
from optparse   import OptionParser
import gzip

def prepare_optparser():
    usage ="""usage: %s [options] sample_input_information_xls

Get chromosomes lists from fai file and then sort it in ascii order.

Then merge it into one file, unzipped.

Example:
    python %s --ref /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.fai --cutSites GCA.GCC.GCT /datd/huboqiang/test_NOM/test/singleC/merge_150610-mES-gWBS-12345/singleC

    """ % (sys.argv[0],sys.argv[0])

    description = "Merge singleC results."

    optparser = OptionParser(
        version="%s v0.1 2015.8" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )

    optparser.add_option(
        "-r", "--ref", default="mm9",
        help="\nReference genome. [default: %default]"
    )
    optparser.add_option(
        "-c", "--cutSites", default="GCA.GCC.GCT",
        help="\nMethylation sites. [default: %default]"
    )
    return optparser

def get_chrom(ref):
    f_ref = open(ref, "r")
    l_chrom = []
    for line in f_ref:
        line = line.strip()
        f    = line.split()
        l_chrom.append(f[0])
    
    return sorted(l_chrom)


def merge_gz(l_infile, outfile, outfile_bdg):
    f_outfile_bdg = open(outfile_bdg, "w")
    f_outfile = open(outfile, "w")
    for infile in l_infile:
        f_in = gzip.open(infile, "rb")
        for line in f_in:
            f_outfile.write(line)
            line = line.strip()
            print >>f_outfile_bdg, "%s" % (line)
        f_in.close()
    
    f_outfile.close()
    f_outfile_bdg.close()

def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        file_prefix = args[0]
        ref = options.ref
        cut_sites= options.cutSites

    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)

    l_chrom = get_chrom(ref)
    outfile = "%s/all.%s.bed" % (file_prefix, cut_sites)
    outfile_bdg = "%s/all.%s.bedGraph" % (file_prefix, cut_sites)
    l_infile = ["%s/%s.%s.bed.gz" % (file_prefix, chrom, cut_sites) for chrom in l_chrom]
    merge_gz(l_infile, outfile, outfile_bdg)
    
if __name__ == '__main__':
    main()