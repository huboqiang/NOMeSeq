# -*- coding: utf-8 -*-
from __future__ import division
import gzip
import os,sys
import numpy as np
import scipy.stats
import tabix
from collections import Counter

def show_help():
    print >>sys.stderr,      "\n\tpython ",sys.argv[0]," /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hheart_Week5_female_R1/singleC/all.GCA.GCC.GCT.NDR.bed /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hheart_Week5_female_R1/singleC/GCA.GCC.GCT.log /date/huboqiang/NORM_seq_PGC/human/02.SingleC.used/rep_in_NDR/hheart_Week5_femaile_R1.nuc.rep.bed /date/huboqiang/NORM_seq_PGC/human/02.SingleC.used/rep_in_NDR/hheart_Week5_femaile_R1.nuc.rep.log /date/huboqiang/NORM_seq_PGC/human/02.SingleC.used/rep/hheart_Week5_female_R1_1/all.GCA.GCC.GCT.bed.gz /date/huboqiang/NORM_seq_PGC/human/02.SingleC.used/rep/hheart_Week5_female_R1_2/all.GCA.GCC.GCT.bed.gz"


def getPval(np_cnt, np_RATIO):
    np_obs = np_cnt
    np_exp = np_RATIO * np.sum(np_obs)
    chisquare = np.sum((np_obs - np_exp)**2/np_exp)
    pval = 0
    if np_obs[0] < np_exp[0]:  
        pval = -1*np.log10(1-scipy.stats.chi2.cdf(chisquare, 1))
    return pval

def main():
    try:
        infile_NDR_bed = sys.argv[1]
        infile_cnt_log = sys.argv[2]
        outfile_NDR_bed = sys.argv[3]
        outfile_cnt_log = sys.argv[4]

        l_GCH_bed = sys.argv[5:]
    except IndexError:
        show_help()
        sys.exit(1)
    
    l_tb_GCH = [ tabix.open(GCH_bed) for GCH_bed in l_GCH_bed ]
    
    f_outfile_NDR_bed = open(outfile_NDR_bed, "w")
    f_outfile_cnt_log = open(outfile_cnt_log, "w")
    
    count_u = 0
    count_m = 0
    with open(infile_cnt_log, "r") as f_infile_cnt_log:
        for line in f_infile_cnt_log:
            line = line.strip()
            f = line.split()
            if f[0] == "all":
                count_u += int(f[1])
                count_m += int(f[2])

    np_RATIO = np.array([count_u/(count_u+count_m), count_m/(count_u+count_m)])
    
    l_indicator = []
    
    with open(infile_NDR_bed, "r") as f_infile_NDR_bed:
        for line in f_infile_NDR_bed:
            line = line.strip()
            f = line.split()
            chrom = f[0]
            begin = int(f[1])
            endin = int(f[2])
            l_cnt = []
            l_pval = []
            for tb_GCH in l_tb_GCH:
                record = tb_GCH.query(chrom, begin, endin)
                np_cnt = np.zeros(2)
                for rec in record:
                    umt_cnt = int(rec[3])
                    met_cnt = int(rec[4])
                    np_cnt[0] += umt_cnt
                    np_cnt[1] += met_cnt
                
                l_cnt.append(np_cnt)
                pval = getPval(np_cnt, np_RATIO)
                l_pval.append(pval)
            
            idx = map(lambda x: "1" if x>=5 else "0", l_pval)
            l_indicator.append( "".join(idx) )
            
            str_umt  = map(lambda x: "%d"    % x[0], l_cnt)
            str_met  = map(lambda x: "%d"    % x[1], l_cnt)
            str_pval = map(lambda x: "%1.2f" % x,    l_pval)
            
            print >>f_outfile_NDR_bed, "%s\t%s\t%s\t%s\t%s\t%s\t%s" % (f[0], f[1], f[2], ",".join(str_umt), ",".join(str_met),",".join(str_pval), "".join(idx))
    
    print >>f_outfile_cnt_log, "\n".join(["%s\t%d" % (key, val) for key, val in sorted(dict(Counter(l_indicator)).iteritems())])
    
if __name__ == '__main__':
    main()
