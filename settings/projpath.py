#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import numpy   as np
import cPickle as pickle
import pandas  as pd

class LoadSamp(object):
    def __init__(self):
        super( LoadSamp,self ).__init__()
        
    def load_RNA_samInfo(self,RNA_infile):
        self.samInfo_RNA_infile  = RNA_infile
        self.samInfo_pd_RNA  = pd.read_csv(self.samInfo_RNA_infile , sep="\t")
        
      
class DirSystem(object):
    def __init__(self):
        super( DirSystem,self ).__init__()
        self.home_dir               = os.path.abspath('./')
        path_file                   = os.path.realpath(__file__)
        self.dir_raw_data           = "%s/00.0.raw_data"               % (self.home_dir)
        self.dir_trim_data          = "%s/00.1.trim_data"              % (self.home_dir)
        self.dir_bismark            = "%s/01.bam"                      % (self.home_dir)
        self.dir_singleC            = "%s/02.SingleC"                  % (self.home_dir)
        self.dir_NDR                = "%s/03.NDR"                      % (self.home_dir)
        self.dir_NDR_flanking       = "%s/04.Flank_plot"               % (self.home_dir)
        self.dir_NDR_Mrg            = "%s/05.NDR_merge"                % (self.home_dir)
        
        self.dir_StatInfo           = "%s/StatInfo"                    % (self.home_dir)
        self.path                   = "%s"     % ("/".join(path_file.split('/')[:-2]))
        self.bin                    = "%s/bin" % ("/".join(path_file.split('/')[:-2]))

        if not os.path.exists( self.dir_StatInfo ):
            os.mkdir( self.dir_StatInfo )
        
        #####################################
        ######## Revise this path!!! ########
        #####################################
        self.Database               = "/data/Analysis/huboqiang/Database_Meth"

        
    def define_scripts(self, s_idx):
        dir_script = "%s/scripts"      % (self.home_dir)
        self.scripts = "%s/scripts/%s" % (self.home_dir, s_idx)

        if not os.path.exists(dir_script):
            os.mkdir(dir_script)

        if not os.path.exists(self.scripts):
            os.mkdir(self.scripts)

        

class UsedSoftware(object):
    def __init__(self):
        super( UsedSoftware,self ).__init__()
        
        ##############################################
        ######## Revise the following path!!! ########
        ##############################################
        self.sftw_py        = "/data/Analysis/huboqiang/software/anaconda/bin/python"
        self.sftw_pl        = "/data/Analysis/huboqiang/lib/local_perl/bin/perl"
        
        self.sftw_samtools  = "/data/Analysis/huboqiang/software/samtools-0.1.18/samtools"
        self.sftw_bedtools  = "/data/Analysis/huboqiang/software/bedtools-2.17.0/bin/bedtools"
        self.sftw_bgzip     = "/data/Analysis/huboqiang/bin/bgzip"
        self.sftw_tabix     = "/data/Analysis/huboqiang/bin/tabix"
        self.sftw_ucsc_dir  = "/data/Analysis/huboqiang/software/UCSC"
        self.sftw_igvtools  = "/data/Analysis/huboqiang/bin/igvtools"
        self.sftw_hommer    = "/data/Analysis/huboqiang/project/human_embryo_sequencing/07.CHIP-seq/test_hommer/bin/findMotifsGenome.pl"
                
        self.sftw_trim      = "/data/Analysis/lilin/bin/trim_galore"
        self.sftw_bismark   = "/data/bin/bismark_v0.7.6"
        self.sftw_changeID  = "/datc/guohongshan/PBAT/scBS/ChangeReadID.pl"
        self.sftw_bowtie_dir= "/data/Analysis/huboqiang/software/bowtie-1.0.0"
        

class ProjInfo(LoadSamp,DirSystem,UsedSoftware):
    def __init__(self):
        super( ProjInfo,self ).__init__()

