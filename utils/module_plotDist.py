#-*- coding:utf-8 -*-
from __future__ import division
import re,sys,os,gzip
import subprocess
import time
import string
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns
sns.set(style="ticks")
import brewer2mpl
from optparse   import OptionParser

##chr	beg	end	strand	repBody_ext1000_inf	repBody_ext1000_bin	genomicBody_ext1000_inf	genomicBody_ext1000_bin	geneBody_ext2000_inf	geneBody_ext2000_bin	geneTSS_ext2000_inf	geneTSS_ext2000_bin	spliceTSS_ext1000_inf	spliceTSS_ext1000_bin	repTSS_ext400_inf	repTSS_ext400_bin	repTES_ext400_inf	repTES_ext400_bin


class RegionInfo(object):
    def __init__(self):
        super(RegionInfo, self).__init__()
        self.type = ['repBody', 'genomicBody', 'geneBody', 'geneTSS', 'spliceTSS', 'repTSS', 'repTES']
        self.ext  = [ 1000,      1000,          5000,       2000,      1000,        400,      400    ]
        self.tag  = ["12",      "12",          "Body",     "TSS",     "12",        "12",     "12"    ]
        self.bin  = [ 301,       301,           301,        201,       201,         201,      201    ]

class SampleInfo(object):
    def __init__(self):
        super(SampleInfo, self).__init__()
        
        # stage level
        # stage level
        self.l_samp_stage = range(len(self.l_samp))
        self.l_stage  = self.l_samp
        self.l_color = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors

class BinInfo(RegionInfo, SampleInfo):
    def __init__(self, l_sam):
        self.l_samp = l_sam
#        print self.l_samp
        super(BinInfo, self).__init__()
        self.__basic_info()
        
    def __basic_info(self):
        self.tag_rat = {}
        self.tag_hit = {}
        self.tag12_rat = {}
        self.tag12_hit = {}

        
        # 4 + 2*len(self.type) + 3 + len(self.l_samp_stage)
    
    def read_files(self, infile, dbfile):
        self.f_infile = open(infile, "r")
        self.f_dbfile = open(dbfile, "r")
        
        while 1:
            line = self.f_infile.readline()
            line_db = self.f_dbfile.readline()
            if len(line_db) == 0 or len(line) == 0:
                break
            
            line_db = line_db.strip()
            f_db    = line_db.split()
            
            line   = line.strip()
            f_rat  = line.split()
            
            np_rat = np.zeros(len(self.l_stage), dtype="float")
            np_hit = np.zeros(len(self.l_stage), dtype="int")

            """
                装载这一位点各修饰类型的ratio。
            """
            
            for j in xrange(0, len(self.l_samp)):
                val  = f_rat[1+j]
                if val == "NA":
                    continue
                else:
                    val = float(val)

                np_rat[ self.l_samp_stage[j] ] += val
                np_hit[ self.l_samp_stage[j] ] += 1
            
            for i,ltype in enumerate(self.type):
#                print 
                tag = self.tag[i]
                idx = f_db[4+2*i+1]
                if idx == "NA":
                    continue
                
                idx = int(idx)
                if tag == "12":
                    tag1 = self.type[i]
                    tag2 = f_db[4+2*i].split("__")[0]
                    tag3 = f_db[4+2*i].split("__")[1]
                    
                    if tag1 not in self.tag12_rat:
                        self.tag12_rat[tag1]    = {}
                        self.tag12_hit[tag1]    = {}
                    if tag2 not in self.tag12_rat[tag1]:
                        self.tag12_rat[tag1][tag2]    = {}
                        self.tag12_hit[tag1][tag2]    = {}

                    if tag3 not in self.tag12_rat[tag1][tag2]:
                        self.tag12_rat[tag1][tag2][tag3]    = np.zeros([len(self.l_stage), self.bin[i] ])
                        self.tag12_hit[tag1][tag2][tag3]    = np.zeros([len(self.l_stage), self.bin[i] ])
                        
                    self.tag12_rat[tag1][tag2][tag3][:,idx]    += np_rat.T
                    self.tag12_hit[tag1][tag2][tag3][:,idx]    += np_hit.T
                    
                else:
                    if tag not in self.tag_rat:
                        self.tag_hit[tag]    = np.zeros([len(self.l_stage), self.bin[i]])
                        self.tag_rat[tag]    = np.zeros([len(self.l_stage), self.bin[i]])
                    
                    self.tag_rat[tag][:,idx] += np_rat.T
                    self.tag_hit[tag][:,idx] += np_hit.T

        self.f_infile.close()
        self.f_dbfile.close()
    
    def out_xls(self, out_prefix):
        if not os.path.isdir(out_prefix):
            os.mkdir(out_prefix)
        
        for tag in self.tag_rat:
            i = self.tag.index(tag)
            f_infile_log    = open("%s/Ratio_plot.%s.%s.xls"  % (out_prefix, self.type[i], tag), "w")
            f_infile_hitlog = open("%s/Hitcnt_plot.%s.%s.xls" % (out_prefix, self.type[i], tag), "w")
            for j, stage in enumerate(self.l_stage):
                np_rat = self.tag_rat[tag][j,:]
                np_hit = self.tag_hit[tag][j,:]
                val = np_rat/np_hit
                print >>f_infile_log,    "%s\t%s" % (stage, "\t".join(np.array(   val, dtype="string")) )
                print >>f_infile_hitlog, "%s\t%s" % (stage, "\t".join(np.array(np_hit, dtype="string")) )
                
            f_infile_log.close()
            f_infile_hitlog.close()
            self.plot_xls("%s/Ratio_plot.%s.%s.xls"  % (out_prefix, self.type[i], tag) )
        
        for tag1 in self.tag12_rat:
            for tag2 in self.tag12_rat[tag1]:
                for tag3 in self.tag12_rat[tag1][tag2]:
                    f_infile_log    = open("%s/Ratio_plot.%s.%s_%s.xls"  % (out_prefix, tag1, tag2, tag3), "w")
                    f_infile_hitlog = open("%s/Hitcnt_plot.%s.%s_%s.xls" % (out_prefix, tag1, tag2, tag3), "w")
                    for i, stage in enumerate(self.l_stage):
                        np_rat = self.tag12_rat[tag1][tag2][tag3][i,:]
                        np_hit = self.tag12_hit[tag1][tag2][tag3][i,:]
                        val = np_rat/np_hit
                        print >>f_infile_log,    "%s\t%s" % (stage, "\t".join(np.array(   val, dtype="string")))
                        print >>f_infile_hitlog, "%s\t%s" % (stage, "\t".join(np.array(np_hit, dtype="string")))
                    
                    f_infile_log.close()
                    f_infile_hitlog.close()
                    self.plot_xls("%s/Ratio_plot.%s.%s_%s.xls"  % (out_prefix, tag1, tag2, tag3))
            
    def plot_xls(self, infile):
        ltype = infile.split("/")[-1].split(".")[1]
        idx   = self.type.index(ltype)
        
        out_pdf  = "%s.pdf" % (".".join(infile.split(".")[:-1]))
        print out_pdf
        
        f_infile = open(infile, "r")
        fig = plt.figure(figsize=(9,5))
        ax  = plt.subplot( 1, 1, 1 )
        i = 0
        length = 0
        for line in f_infile:
            line = line.strip()
            f    = line.split()
            length = len(f)
            np_val = np.array(f[1:],dtype="float")

            ax.plot( range(0, len(f)-1), np_val, linewidth=1, label=f[0], color=self.l_color[i] )
            i += 1
    
        ax.legend(loc="upper right",fontsize=9)
        
        
    
        x_pos = [     0,      50,    100,   150,   199 ]
        x_lab = ["-%d" % self.ext[idx],"-%d" % int(self.ext[idx]/2),    "0", "%d" % int(self.ext[idx]/2), "%d" % self.ext[idx]]
    
        if self.bin[idx] == 301:
            x_pos = [     0,     50,     100,   150,   200,   250,   299 ]
            x_lab = ["-%d" % self.ext[idx],"-%d" % int(self.ext[idx]/2),   "5'", "50%",  "3'", "%d" % int(self.ext[idx]/2), "%d" % self.ext[idx]]
                
        ax.get_xaxis().set_ticks(x_pos)
        ax.get_xaxis().set_ticklabels(x_lab)
#        ax.set_ylim(0,0.15)
        ax.set_ylabel('Normalized m6A occupancy')
        ax.set_xlabel("Distance from TSS(bp)")
#        sns.despine(trim=True)
    
        fig.savefig(out_pdf, format="pdf")

            
        
        


def prepare_optparser():
    usage ="""usage: %s [options] 
    
    Using -h or --help for more information
    
Example:
     python %s /date/huboqiang/NORM_seq_PGC/mouse/StatInfo/MouseSample_list_20150825_meth.tab/ACG.TCG/RatioMatrix.xls /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG.anno.bed /date/huboqiang/NORM_seq_PGC/mouse/StatInfo/MouseSample_list_20150825_meth.tab/ACG.TCG/plot /date/huboqiang/NORM_seq_PGC/mouse/StatInfo/MouseSample_list_20150825_meth.tab.ACG.TCG.allsamp.xls
    
    """ % (sys.argv[0],sys.argv[0])
    
    description = "Calculating Distribution plot for methylation level. "
    
    optparser = OptionParser(version="%s v0.2 20150124" % (sys.argv[0]),description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",      action="help",       help="\nShow this help message and exit.")
    return optparser

def main():
    prepare_optparser()
    (options, args) = prepare_optparser().parse_args()
    try:
        infile      = args[0]
        dbfile      = args[1]
        out_prefix  = args[2]
        samfile     = args[3]
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)
    
    print samfile
    f_samfile = open(samfile, "r")
    l_sam = []
    f_samfile.readline()
    for line in f_samfile:
        line = line.strip()
        f = line.split()
        l_sam.append(f[0])
    
    f_samfile.close()
    
    print l_sam
    m_process = BinInfo(l_sam)
    m_process.read_files(infile, dbfile)
    m_process.out_xls(out_prefix)
    
if __name__ == '__main__':
    main()
    
