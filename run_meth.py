#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os

import cPickle as pickle
import numpy as np
from optparse   import OptionParser

import MethGC.frame.module01_mapping_from_raw   as m01
import MethGC.frame.module02_chromatinNetwork2   as m02
import MethGC.utils.module_create_database  as m_db

def prepare_optparser():
    usage ="""usage: %s [options] sample_input_information_xls

Reference fasta file and refGene file should be all put in dictionary
self.Database, which were defined in settings/projpath.py

If not put in self.Database, this program will downloading from UCSC.
Detail information could be reached in utils/module_create_database.py

Suport genome includes:
    http://hgdownload.soe.ucsc.edu/goldenPath


Using -h or --help for more information

Example:
    python %s --ref hg19 --cutSites GCA.GCC.GCT,ACG.TCG samp_test.xls

    """ % (sys.argv[0],sys.argv[0])

    description = " NOM-seq analysis pipeline "

    optparser = OptionParser(
        version="%s v0.1 2015.7" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )

    optparser.add_option(
        "-r", "--ref", default="hg19",
        help="\nReference genome. [default: %default]"
    )

    optparser.add_option(
        "-c", "--cutSites", default="GCA.GCC.GCT,ACG.TCG",
        help="\nMethylation sites. [default: %default]"
    )
    return optparser

def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        sam_Info = args[0]
        ref = options.ref
        l_cut_sites= options.cutSites.split(",")
        
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)
    
#    part0 = m_db.DataBaseInit(ref, sam_Info, is_debug = 0)
#    part0.check_files(l_cut_sites)
    
    part1 = m01.MapFromRaw(ref, sam_Info, l_cut_sites, is_debug = 0)
    Depth_Pos = 3
    part1.s01_Trim()
    part1.s02_Bismark()
    part1.s03_Bam2SingleC()
    part1.s04_statMeth(Depth_Pos)
    part1.s05_NDR()
    part1.s06_merge_singleC()
    part1.s07_mergeSample(Depth_Pos)
    part1.s08_plotDist()
    part1.s09_NDR_IGV()
    part1.s10_NDR_flanking()
    part1.s11_NDR_motif()

    Depth_Pos = 1
    part1.s01_Trim()
    part1.s02_Bismark()
    part1.s03_Bam2SingleC()
    part1.s04_statMeth(Depth_Pos)
    part1.s05_NDR()
    part1.s06_merge_singleC()
    part1.s07_mergeSample(Depth_Pos)

if __name__ == '__main__':
    main()
