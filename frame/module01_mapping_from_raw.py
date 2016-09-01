#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import MethGC.utils.module_running_jobs as m_jobs
import MethGC.settings.projpath         as m_proj
import MethGC.settings.scripts          as m_scpt

def make_dir(l_args):
    """docstring for make_dir"""
    if len(l_args) == 1:
        """ only dictionary """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
    
    elif len(l_args) == 2:
        """ dictionary/sample """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s" % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s" % (l_args[0],l_args[1]) )

    elif len(l_args) == 3:
        """ dictionary/sample/subname """
        if not os.path.isdir( l_args[0] ):
            os.mkdir( l_args[0] )
        if not os.path.isdir( "%s/%s"    % (l_args[0],l_args[1]) ):
            os.mkdir(         "%s/%s"    % (l_args[0],l_args[1]) )
        if not os.path.isdir( "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) ):
            os.mkdir(         "%s/%s/%s" % (l_args[0],l_args[1],l_args[2]) )


class MapFromRaw(m_scpt.Scripts):
    def __init__(self, ref, sam_RNAinfo, l_sites=["GCA.GCC.GCT", "ACG.TCG"], is_debug=1):
        super(MapFromRaw, self).__init__()
        
        self.s_idx = ".".join(sam_RNAinfo.split("/")[-1].split(".")[:-1])
        self.load_RNA_samInfo(sam_RNAinfo)
        self.define_scripts(self.s_idx)
        self.l_query = l_sites
        self.ref = ref
        self.is_debug = is_debug
        self.define_files(self.ref)
        self.M_CvtEnd = {"PE":2,"SE":1}
        
    def s01_Trim(self):
        sh_file      = "%s/s01.Trim.sh"      % (self.scripts)
        sh_work_file = "%s/s01.Trim_work.sh" % (self.scripts)

        l_sh_info = self.s_01_trim()
        l_sh_work = []
        
        for samp in self.samInfo_pd_RNA['sample']:
            print samp
            idx        =(self.samInfo_pd_RNA['sample'] == samp)
            brief_name = self.samInfo_pd_RNA[ idx ]['brief_name'].values[0]

            make_dir([self.dir_trim_data, brief_name])
            for i in os.walk( '%s/%s' % (self.dir_raw_data, samp)):
                in_fq1 = "%s/%s" % (i[0],i[2][0] )
                in_fq2 = "%s/%s" % (i[0],i[2][1] )
                if in_fq1.split('_')[-2] == "R2":
                    in_fq2 = "%s/%s" % (i[0],i[2][0] )
                    in_fq1 = "%s/%s" % (i[0],i[2][1] )
                
                l_sh_work.append("sh %s %s %s %s" % \
                    (sh_file, in_fq1, in_fq2, brief_name)
                )
                
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)

    def s02_Bismark(self):
        sh_file      =  "%s/s02.Bismark.sh"      % (self.scripts)
        sh_work_file =  "%s/s02.Bismark_work.sh" % (self.scripts)
        
        l_sh_info = self.s_02_bismark()
        l_sh_work = []
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            print brief_name
            idx    =(self.samInfo_pd_RNA['brief_name'] == brief_name)
            sample = self.samInfo_pd_RNA[idx]['sample'].values[0]
            make_dir([self.dir_bismark, brief_name])
            for i in os.walk( '%s/%s' % (self.dir_raw_data, sample)):
                in_fq1 = "%s_val_1.fq.gz" % (i[2][0].split(".")[0] )
                in_fq2 = "%s_val_2.fq.gz" % (i[2][1].split(".")[0] )
                if in_fq1.split('_')[-2] == "R2":
                    in_fq2 = "%s_val_2.fq.gz" % (i[2][0].split(".")[0] )
                    in_fq1 = "%s_val_1.fq.gz" % (i[2][1].split(".")[0] )
                
                l_sh_work.append("sh %s %s %s %s %s" % \
                    (sh_file, brief_name, in_fq1, in_fq2, self.ref)
                )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=3, is_debug = self.is_debug)

    def s03_Bam2SingleC(self):
        sh_file      =  "%s/s03.Bam2SingleC.sh"      % (self.scripts)
        sh_work_file =  "%s/s03.Bam2SingleC_work.sh" % (self.scripts)
        
        l_sh_info = self.s_03_Bam2SingleC()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        self.l_chrom   = []
        f_fai_input = open(fai_input, "r")
        for line in f_fai_input:
            line = line.strip("\n")
            f    = line.split()
            self.l_chrom.append(f[0])
        f_fai_input.close()
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            make_dir([self.dir_singleC, brief_name, "bam"])
            make_dir([self.dir_singleC, brief_name, "singleC"])
            for chrom in self.l_chrom:
                l_sh_work.append("sh %s %s %s %s %s" % (sh_file, 
                    brief_name, chrom, self.ref, " ".join(self.l_query))
                )
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)

    def __load_ChrInfo(self):
        self.M_QuerySamChr = {}
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            for query in self.l_query:
                infile = "%s/%s/%s/%s.log" % (
                    self.dir_singleC, brief_name, "singleC", query
                )
                
                if query not in self.M_QuerySamChr:
                    self.M_QuerySamChr[query] = {}
                
                if brief_name not in self.M_QuerySamChr[query]:
                    self.M_QuerySamChr[query][brief_name] = {}
                
                f_infile = open(infile, "r")
                for line in f_infile:
                    line = line.strip()
                    f = line.split()
                    self.M_QuerySamChr[query][brief_name][f[0]] = {
                        'umt':int(f[1]), 'met':int(f[2]), 'cnt':int(f[3])
                    }
#                    print brief_name, query, self.M_QuerySamChr[query][brief_name]
                f_infile.close()
    
    def s04_statMeth(self, depth=5):
        sh_file      =  "%s/s04.StatLog.sh"      % (self.scripts)
        sh_work_file =  "%s/s04.StatLog_work.sh" % (self.scripts)
        
        l_sh_info = self.s_04_StatLog()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            for query in self.l_query:
                l_sh_work.append("sh %s %s %s %d %s" % (sh_file, brief_name,
                     self.ref, depth, query))

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
        
        self.__load_ChrInfo()
        print self.M_QuerySamChr
        
        for query in self.l_query:
            outfile = "%s/%s.%s.allsamp.xls" % (self.dir_StatInfo, self.s_idx, query)
            f_outfile=open(outfile, "w")
            print >>f_outfile, "sample\tall_umt\tall_met\tall_cnt\tlambda_umt\tlambda_met\tlambda_cnt"
            for brief_name in self.samInfo_pd_RNA['brief_name']:
                print >>f_outfile, "%s\t%d\t%d\t%d\t%d\t%d\t%d" % (brief_name, 
                    self.M_QuerySamChr[query][brief_name]['all']['umt'],
                    self.M_QuerySamChr[query][brief_name]['all']['met'],
                    self.M_QuerySamChr[query][brief_name]['all']['cnt'],
                    self.M_QuerySamChr[query][brief_name]['lambda']['umt'],
                    self.M_QuerySamChr[query][brief_name]['lambda']['met'],
                    self.M_QuerySamChr[query][brief_name]['lambda']['cnt'])
            
            f_outfile.close()
    
    def s05_NDR(self, pvalue = 15):
        sh_file      =  "%s/s05.NDR.sh"      % (self.scripts)
        sh_work_file =  "%s/s05.NDR_work.sh" % (self.scripts)
    
        l_sh_info = self.s_05_NDR()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        self.l_chrom   = []
        f_fai_input = open(fai_input, "r")
        for line in f_fai_input:
            line = line.strip("\n")
            f    = line.split()
            self.l_chrom.append(f[0])
        f_fai_input.close()
        
        self.__load_ChrInfo()
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            for query in self.l_query:
                for chrom in self.l_chrom:
                    l_sh_work.append("sh %s %s %s %s %d %d %d" % (sh_file, 
                        brief_name, query, chrom, 
                        self.M_QuerySamChr[query][brief_name]['all']['umt'],
                        self.M_QuerySamChr[query][brief_name]['all']['met'],
                        pvalue )
                    )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)

    def s06_merge_singleC(self):
        sh_file      =  "%s/s06.mergeSC.sh"      % (self.scripts)
        sh_work_file =  "%s/s06.mergeSC_work.sh" % (self.scripts)
    
        l_sh_info = self.s_06_mergeSC()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            for query in self.l_query:
                l_sh_work.append("sh %s %s %s %s" % (sh_file, brief_name, 
                    query, self.ref)
                )

        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
    
    def s07_mergeSample(self, depth=5):
        for query in self.l_query:
            make_dir([ self.dir_StatInfo, self.s_idx, query ])
            mergeRatio = "%s/%s/%s/RatioMatrix.xls" % (self.dir_StatInfo, self.s_idx, query)
            f_mergeRatio = open(mergeRatio, "w")
            l_file = [ "%s/%s/singleC/all.%s.bed" % 
                (self.dir_singleC, sam, query) 
                for sam in self.samInfo_pd_RNA['brief_name']
            ]
            
            print >>f_mergeRatio, "chrpos\t%s" % ("\t".join(self.samInfo_pd_RNA['brief_name']))
            shell_info = " paste %s " % (" ".join(l_file))
            p=subprocess.Popen(shell_info,stdout=subprocess.PIPE,shell=True)
            for line in p.stdout:
                line = line.strip('\n')
                f   = line.split()
                info = "%s:%s-%s" % (f[0], f[1], f[2])
                l_out = []
                for i in xrange(len(self.samInfo_pd_RNA['brief_name'])):
                    cnt_umt = int(f[3+i*5])
                    cnt_met = int(f[4+i*5])
                    val = "NA"
                    if cnt_umt+cnt_met >= depth:
                        val = str(cnt_met/(cnt_met+cnt_umt))
                    l_out.append(val)
                
                print >>f_mergeRatio, "%s\t%s" % (info, "\t".join(l_out))
            
            f_mergeRatio.close()
        
    def s08_plotDist(self):
        sh_file      =  "%s/s08.plotDist.sh"      % (self.scripts)
        sh_work_file =  "%s/s08.plotDist_work.sh" % (self.scripts)
        
        l_sh_info = self.s_08_plotDist()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        
        for query in self.l_query:
            mergeRatio = "%s/%s/%s/RatioMatrix.xls" % (
                self.dir_StatInfo, self.s_idx, query
            )
            db_file = "%s/%s/%s_lambda.fa.%s.anno.bed" % (
                self.Database, self.ref, self.ref, query
            )
            outprefix = "%s/%s/%s/plot" % (
                self.dir_StatInfo, self.s_idx, query
            )
            sam_file = "%s/%s.%s.allsamp.xls" % (
                self.dir_StatInfo, self.s_idx, query
            )
            l_sh_work.append("sh %s %s %s %s %s" % 
                (sh_file, mergeRatio, db_file, outprefix, sam_file)
            )
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
        
        
    def s09_NDR_IGV(self):
        sh_file      =  "%s/s09.NDR_IGV.sh"      % (self.scripts)
        sh_work_file =  "%s/s09.NDR_IGV_work.sh" % (self.scripts)
        
        l_sh_info = self.s_09_NDR_IGV()
        l_sh_work = []
        
        fai_input = "%s/%s/%s_lambda.fa.fai" % (
            self.Database, self.ref, self.ref
        )
        make_dir([self.dir_NDR])        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            for query in self.l_query:
                l_sh_work.append("sh %s %s %s %s" % (sh_file, brief_name, 
                    query, self.ref)
                )
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
        
    
    def s10_NDR_flanking(self):
        sh_file      =  "%s/s10.NDR_flanking.sh"      % (self.scripts)
        sh_work_file =  "%s/s10.NDR_flanking_work.sh" % (self.scripts)
        
        l_sh_info = self.s_10_NDR_flanking()
        l_sh_work = []
        
        make_dir([self.dir_NDR_flanking])
        make_dir([self.dir_NDR, "distal_region"])
        make_dir([self.dir_NDR, "TSS_region"])
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            l_sh_work.append("sh %s %s %s" % (sh_file, brief_name, self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)
    
    def s11_NDR_motif(self):
        sh_file      =  "%s/s11.NDR_motif.sh"      % (self.scripts)
        sh_work_file =  "%s/s11.NDR_motif_work.sh" % (self.scripts)
        
        l_sh_info = self.s_11_NDR_motif()
        l_sh_work = []
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            l_sh_work.append("sh %s %s %s" % (sh_file, brief_name, self.ref))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)


    def s09_2_Het_long(self):
        sh_file      =  "%s/s09.2.Het_long.sh"      % (self.scripts)
        sh_work_file =  "%s/s09.2.Het_long_work.sh" % (self.scripts)
        
        l_sh_info = self.s_09_2_HetLong_motif()
        l_sh_work = []
        
        for brief_name in self.samInfo_pd_RNA['brief_name']:
            l_sh_work.append("sh %s %s" % (sh_file, brief_name))
        
        my_job = m_jobs.run_jobs(sh_file, sh_work_file, l_sh_info, l_sh_work)
        my_job.running_multi(cpu=8, is_debug = self.is_debug)


    
