# -*- coding: utf-8 -*-
from __future__ import division
import gzip
import os,sys
import numpy as np
import scipy.stats
import tabix

def show_help():
   print >>sys.stderr,      "\n\tpython ",sys.argv[0]," /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/singleC/chr10.ACG.TCG.bed.gz COUNT_U COUNT_M 15"

class NDR(object):
    def __init__(self, chrom, count_u, count_m, tb_file, cutoff = 1e-5, depth = 3, bin_len=100, step_len=20, dist_len=140 ):
        self.chrom    = chrom
        self.bin_len  = bin_len
        self.step_len = step_len
        self.dist_len = dist_len
        self.cutoff   = cutoff
        self.depth = depth

        self.np_RATIO = np.array( [count_u/(count_u+count_m), count_m/(count_u+count_m)] )
        self.tb_file = tb_file
        self.NDR_info = {}
        
        self.l_umt_met = {}
        self.l_rat_umt_met = {}
        
    def parse_line(self, pos_beg, cnt_umt, cnt_met):
        self.__get_bins(pos_beg)
        self.__NDR_append(cnt_umt, cnt_met, pos_beg)
        self.__NDR_bin_detect()
    
    def __NDR_append(self, cnt_umt, cnt_met, pos_beg):
        """Load count in a bin"""
        for i in range(self.bin_idx_min, self.bin_idx_max+1):
            if i not in self.l_umt_met:
                self.l_umt_met[i] = np.zeros(2, dtype="int")
                self.l_rat_umt_met[i] = []

            self.l_umt_met[i][0] += cnt_umt
            self.l_umt_met[i][1] += cnt_met
            ratio = cnt_met/(cnt_umt + cnt_met)
            self.l_rat_umt_met[i].append(ratio)
    
    def __NDR_bin_detect(self):
        """检视所有 cover 的 bins."""
        for idx in sorted(self.l_umt_met):
            self.__bin_info(idx)
            
            """对这种可能显著并且不会继续落 reads 的 bin"""
            if idx < self.bin_idx_min:
                np_obs    = self.l_umt_met[idx]
                np_exp    = self.np_RATIO * np.sum(np_obs)
                chisquare = np.sum((np_obs - np_exp)**2/np_exp)
                pval = -1*np.log10(1-scipy.stats.chi2.cdf(chisquare, 1))
                if np_obs[1] < np_exp[1]:  
                    pval = -1 * pval
                
                bin_up   = int((self.bin_begin + self.bin_endin)/2 - self.step_len*10)
                bin_down = int((self.bin_begin + self.bin_endin)/2 + self.step_len*10)
                np_val_ext = np.zeros(2)
                l_rat_umt_met_ext = []
                try:
                    record = self.tb_file.query(self.chrom, bin_up, bin_down)
                    for rec in record:
                        total = int(rec[3]) + int(rec[4])
                        if total >= self.depth:
                            np_val_ext[0] += int(rec[3])
                            np_val_ext[1] += int(rec[4])
                            ratio = int(rec[4])/(int(rec[3]) + int(rec[4]))
                            l_rat_umt_met_ext.append(ratio)
                except:
                    pass
                
                bin_sur1_up   = self.bin_begin - 60
                bin_sur1_down = self.bin_begin
                bin_sur2_up   = self.bin_endin
                bin_sur2_down = self.bin_endin + 60
                
                
                l_rat_umt_met_sur1 = []
                l_rat_umt_met_sur2 = []
                try:
                    record = self.tb_file.query(self.chrom, bin_sur1_up, bin_sur1_down)
                    for rec in record:
                        total = int(rec[3]) + int(rec[4])
                        if total >= self.depth:
                            ratio = int(rec[4])/(int(rec[3]) + int(rec[4]))
                            l_rat_umt_met_sur1.append(ratio)
                except:
                    pass
                
                try:
                    record = self.tb_file.query(self.chrom, bin_sur2_up, bin_sur2_down)
                    for rec in record:
                        total = int(rec[3]) + int(rec[4])
                        if total >= self.depth:
                            ratio = int(rec[4])/(int(rec[3]) + int(rec[4]))
                            l_rat_umt_met_sur2.append(ratio)
                except:
                    pass
                
                
                np_obs_ext    = np_val_ext
                np_exp_ext    = self.np_RATIO * np.sum(np_obs_ext)
                val_reg = np.array(self.l_rat_umt_met[idx]).mean()
                val_ext = np.array(l_rat_umt_met_ext).mean()
                val_sur1 = np.array(l_rat_umt_met_sur1).mean()
                val_sur2 = np.array(l_rat_umt_met_sur2).mean()
                
                chisquare_ext = np.sum((np_obs_ext - np_exp_ext)**2/np_exp_ext)
                pval_ext = -1*np.log10(1-scipy.stats.chi2.cdf(chisquare_ext, 1))
                if np_obs_ext[1] < np_exp_ext[1]:  
                    pval_ext = -1 * pval_ext
                
                pval_ttest = -1*np.log10(scipy.stats.ttest_ind(self.l_rat_umt_met[idx], l_rat_umt_met_ext, equal_var=False)[1])
                if val_reg < val_ext:
                    pval_ttest = -1 * pval_ttest

                pval_ttest2l = -1*np.log10(scipy.stats.ttest_ind(self.l_rat_umt_met[idx], l_rat_umt_met_sur1, equal_var=False)[1])
                if val_reg < val_sur1:
                    pval_ttest2l = -1 * pval_ttest2l

                pval_ttest2r = -1*np.log10(scipy.stats.ttest_ind(self.l_rat_umt_met[idx], l_rat_umt_met_sur2, equal_var=False)[1])
                if val_reg < val_sur2:
                    pval_ttest2r = -1 * pval_ttest2r

                
                print "%s\t%d\t%d\t%d\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%d,%d\t%d,%d\t%1.2f,%1.2f\t%1.2f,%1.2f\t%1.2f,%1.2f,%1.2f,%1.2f" % (self.chrom, self.bin_begin, self.bin_endin,  (self.bin_begin+self.bin_endin)/2, pval, pval_ext, pval_ttest, pval_ttest2l, pval_ttest2r, self.l_umt_met[idx][0], self.l_umt_met[idx][1], np_val_ext[0], np_val_ext[1], np_exp[0], np_exp[1], np_exp_ext[0], np_exp_ext[1], val_reg, val_ext, val_sur1, val_sur2)
#                print "%s\t%d\t%d\t%1.2f\t%s\t%s" % (self.chrom, self.bin_begin, self.bin_endin, pval_ttest, val_reg, val_ext)
                del self.l_umt_met[idx]
                del self.l_rat_umt_met[idx]
    

    def __get_bins(self, pos):
        """pos to bin"""
        self.bin_idx_min = int( (pos - self.bin_len)/self.step_len )
        self.bin_idx_max = int( pos/self.step_len )
#        print pos, self.bin_idx_min, (self.bin_idx_min) * self.step_len + 1, (self.bin_idx_min+1) * self.step_len + self.bin_len
#        print pos, self.bin_idx_max, (self.bin_idx_max) * self.step_len + 1, (self.bin_idx_max+1) * self.step_len + self.bin_len
        if self.bin_idx_min < 0:
            self.bin_idx_min = 0
        
    def __bin_info(self, in_idx):
        """bin to pos"""
        self.bin_begin = (in_idx) * self.step_len + 1
        self.bin_endin = (in_idx+1) * self.step_len + self.bin_len
#        print "BinInfo", in_idx, self.bin_begin, self.bin_endin

def main():
    try:
        infile = sys.argv[1]
        COUNT_U= int(sys.argv[2])
        COUNT_M= int(sys.argv[3])
        CUTOFF = int(sys.argv[4])
    except IndexError:
        show_help()
        sys.exit(1)
    
#    l_line=  f_cntfile.readlines()
    
    CHROM =infile.split("/")[-1].split(".")[0]
    
    DEPTH = 3
    tb_file = tabix.open(infile)
    
    m_NDR =  NDR(CHROM, COUNT_U, COUNT_M, tb_file, CUTOFF, DEPTH, bin_len=100, step_len=20, dist_len=140 )
    
    f_infile = gzip.open(infile, "rb")
    for line in f_infile:
        line = line.strip("\n")
        f    = line.split()

        pos_beg = int(f[1])
        cnt_umt = int(f[3])
        cnt_met = int(f[4])
        cnt_tot = cnt_umt + cnt_met
        
        if cnt_tot < DEPTH:
            continue
        
        m_NDR.parse_line(pos_beg, cnt_umt, cnt_met)
    
    f_infile.close()
    
if __name__ == '__main__':
    main()
