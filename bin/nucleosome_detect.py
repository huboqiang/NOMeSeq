# -*- coding: utf-8 -*-
from __future__ import division
import gzip
import os,sys
import numpy as np
import scipy.stats
import tabix



def show_help():
   print >>sys.stderr, "\n\tpython ",sys.argv[0]," /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/singleC/chr10.ACG.TCG.bed.gz COUNT_U COUNT_M 15"

class NDR(object):
    def __init__(self, chrom, count_u, count_m, cutoff, tb_file, bin_len=100, step_len=20, dist_len=140 ):
        self.chrom    = chrom
        self.bin_len  = bin_len
        self.step_len = step_len
        self.dist_len = dist_len
        self.cutoff   = cutoff

        self.np_RATIO = np.array( [count_u/(count_u+count_m), count_m/(count_u+count_m)] )
        self.tb_file = tb_file
        self.NDR_info = {}
        
        self.l_umt_met = {}
        
    def parse_line(self, pos_beg, cnt_umt, cnt_met):
        self.__get_bins(pos_beg)
        self.__NDR_append(cnt_umt, cnt_met, pos_beg)
        self.__NDR_bin_detect()
    
    def last_line(self):
        if len(self.NDR_info) != 0:
            if self.__NDR_len() > self.dist_len:
                self.__NDR_report()

    def __NDR_append(self, cnt_umt, cnt_met, pos_beg):
        """Load count in a bin"""
        for i in range(self.bin_idx_min, self.bin_idx_max+1):
            if i not in self.l_umt_met:
                self.l_umt_met[i] = np.zeros(3, dtype="int")

            self.l_umt_met[i][0] += cnt_umt
            self.l_umt_met[i][1] += cnt_met
            self.l_umt_met[i][2] += 1
            
    
    def __NDR_bin_detect(self):
        """检视所有 cover 的 bins."""
        for idx in sorted(self.l_umt_met):
            self.__bin_info(idx)
            
            """对这种可能显著并且不会继续落 reads 的 bin"""
            if idx < self.bin_idx_min:
                """确定NDR和现有bin 不overlap，看看是否足够长，打印出来，否则可以扔掉"""
                if len(self.NDR_info) != 0 and self.bin_begin > self.NDR_info['end']:
                    if self.__NDR_len() >= self.dist_len-1:
                        self.__NDR_report()
                            
                    self.__NDR_remove()
                        
                    del self.l_umt_met[idx]
                    continue
                                    
                pval, pval_local = self.__chisq_bin(idx)
                """如果bin显著，无NDR则该bin作为NDR的基点，有NDR则和现有bin必有overlap（否则上一步被截住），对NDR elongate即可"""
                if max(pval, pval_local)< self.cutoff:
#                    self.__NDR_remove()
                    del self.l_umt_met[idx]
                    continue
                else:
                    if 'beg' not in self.NDR_info:
                        self.NDR_info = {
                            'beg': self.bin_begin, 'end': self.bin_endin, 'idx':[idx], 'pval':[pval], 'pval_local':[pval_local], 
                            'umt':[self.l_umt_met[idx][0]], 'met':[self.l_umt_met[idx][1]], 'cnt':[self.l_umt_met[idx][2]]
                        }
                    else:
                        if idx not in self.NDR_info['idx']:
                            self.NDR_info['end'] = self.bin_endin
                            self.NDR_info['idx'].append(idx)
                            self.NDR_info['pval'].append(pval)
                            self.NDR_info['pval_local'].append(pval_local)
                            self.NDR_info['umt'].append(self.l_umt_met[idx][0])
                            self.NDR_info['met'].append(self.l_umt_met[idx][1])
                            self.NDR_info['cnt'].append(self.l_umt_met[idx][2])
                        
                    
#                    print self.NDR_info

    
    def __NDR_remove(self):
        if 'idx' in self.NDR_info:
            for idx_rm in self.NDR_info['idx']:
                if idx_rm in self.l_umt_met:
                    del self.l_umt_met[idx_rm]
        
        self.NDR_info = {}
#        del self.l_umt_met[idx]
    
    def __NDR_len(self):
        return self.NDR_info['end'] - self.NDR_info['beg']
    
    def __NDR_report(self):
#        str_p = ",".join( ["%1.2e" % (pval) for pval in self.NDR_info['pval']] )

#        l_val = np.unique(np.array(sorted(["%d" % (idx) for idx in self.NDR_info['idx']]),dtype="int"))
#        np_idx = np.array([ self.NDR_info['idx'].index(i) for i in l_val ],dtype="int")

        np_idx = np.array(self.NDR_info['idx'], dtype="int")

        np_pval= np.array(self.NDR_info['pval'])
        np_pval_local= np.array(self.NDR_info['pval_local'])
        np_umt = np.array(self.NDR_info['umt'])
        np_met = np.array(self.NDR_info['met'])
        np_cnt = np.array(self.NDR_info['cnt'])
        str_p = ",".join( ["%1.2f" % (pval) for pval in np_pval] )
        str_pl = ",".join( ["%1.2f" % (pval) for pval in np_pval_local] )
        str_u = ",".join( ["%d" % (umt)  for umt  in np_umt ] )
        str_m = ",".join( ["%d" % (met)  for met  in np_met ] )
        str_r = ",".join( ["%d" % (np_met[i]/(np_met[i]+np_umt[i])*100)  for i in xrange(len(np_umt)) ] )
        str_c = ",".join( ["%d" % (cnt)  for cnt  in np_cnt ] )

#        if np_pval.max() > 3 and np_pval_local.max() > 5 and np_pval_local.mean() > 3:
        print "%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (self.chrom, self.NDR_info['beg'], self.NDR_info['end'], str_p, str_pl, str_u, str_m, str_c, str_r)

    def __chisq_bin(self, idx):
        np_obs    = self.l_umt_met[idx][0:2]
        np_exp    = self.np_RATIO * np.sum(np_obs)
        chisquare = np.sum((np_obs - np_exp)**2/np_exp)
        pval = 1
        pval_local = 1
        '''
            只考虑 umet的reads 足够多的情况，这种情况下 reads 全未被甲基化，即被组蛋白所覆盖，nucleosome.
        '''
        if np_obs[0] > np_exp[0]:
            bin_up   = int((self.bin_begin + self.bin_endin)/2 - self.step_len*10)
            bin_down = int((self.bin_begin + self.bin_endin)/2 + self.step_len*10)
            np_val = np.zeros(2)
            try:
                record = self.tb_file.query(self.chrom, bin_up, bin_down)
                for rec in record:
                    np_val[0] += int(rec[3])
                    np_val[1] += int(rec[4])
            except:
                pass
            
            np_exp_local = np_val/np_val.sum()*np_obs.sum()
            chisquare_local = np.sum((np_obs - np_exp_local)**2/np_exp_local)
            pval = -1*np.log10(1-scipy.stats.chi2.cdf(chisquare, 1))
            pval_local = -1*np.log10(1-scipy.stats.chi2.cdf(chisquare_local, 1))
#        print self.bin_begin, self.bin_endin, np_obs, np_exp, pval, idx
#        print np_obs[0] / np_exp[0], self.bin_begin, self.bin_endin
        return pval, pval_local
        
    

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
    tb_file = tabix.open(infile)
    DEPTH = 3
#    CUTOFF = 10**(-1*CUTOFF)
    m_NDR =  NDR(CHROM, COUNT_U, COUNT_M, CUTOFF, tb_file, bin_len=40, step_len=20, dist_len=60 )
    
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
    
    m_NDR.last_line()        
    
if __name__ == '__main__':
    main()
