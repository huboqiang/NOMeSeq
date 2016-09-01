#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import gzip
import subprocess
import time
import string

import numpy as np
import pandas as pd
from scipy import stats
import functools
from pymongo import MongoClient
from bson.objectid import ObjectId
from optparse import OptionParser
import logging
import tabix

logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    stream=sys.stderr
)

class DbGene(object):
    """docstring for DbGene"""
    def __init__(self):
        super(DbGene, self).__init__()

    def loadExp_from_file(self,
            infile_exp = "/data/Analysis/huboqiang/tmp/lilin_20150310/05.cuffnorm/genes.fpkm_table.rev.clean"
        ):
        self.pd_exp = pd.read_csv(infile_exp, sep="\t", index_col=[0])
        self.__load_samInfo()
        self.__load_exp()
    
    def __load_samInfo(self):
        self.struct = {}
        l_samples = self.pd_exp.columns
        M_samInfo = {}
        for i,sam in enumerate(l_samples):
            sex = sam.split("_")[0]
            tissue = sam.split("_")[1]
            stage = sam.split("_")[2]
            if sex == "ICM":
                sex = "U"
                tissue = sam.split("_")[0]
                stage = sam.split("_")[1]

            if tissue not in self.struct:
                self.struct[tissue] = {}
            if stage not in self.struct[tissue]:
                self.struct[tissue][stage] = {}
            if sex not in self.struct[tissue][stage]:
                self.struct[tissue][stage][sex] = []
            
            self.struct[tissue][stage][sex].append(i)
    
    def __load_exp(self):
        for gene in self.pd_exp.index:
            gene_rev = gene.replace(".", "_")
            INFO = {'gene': gene_rev}
            for tissue in self.struct:
                for stage in self.struct[tissue]:
                    for sex in self.struct[tissue][stage]:
                        np_idx = np.array(self.struct[tissue][stage][sex], dtype="int")
                        print np_idx
                        tag = "%s__%s__%s" % (tissue, stage, sex)
                        INFO[tag] = list(self.pd_exp.loc[gene].values[np_idx])
            
            self.posts_PGC_exp.insert_one(INFO)
        
        self.posts_PGC_exp.create_index([("gene", 1)])


    def loadExpIsoform_from_file(self,
            infile_exp = "/date/huboqiang/NORM_seq_PGC/RNA_more/isoformInfoRNA/human/mergeIsoformFPKM.xls"
        ):
        self.pd_exp = pd.read_csv(infile_exp, sep="\t", index_col=[0, 1, 2])
        logging.info("%s Reading Done." % (infile_exp))
        self.__load_samInfo()
        self.__load_exp_isoform()
    
    def __load_exp_isoform(self):
        for gene in self.pd_exp.index:
            gene_rev = gene[0]
            TSS = gene[1]
            iso = gene[2]
            INFO = {'gene': gene_rev, "TSS" : TSS, "Isoform" : iso}
            for tissue in self.struct:
                for stage in self.struct[tissue]:
                    for sex in self.struct[tissue][stage]:
                        np_idx = np.array(self.struct[tissue][stage][sex], dtype="int")
                        tag = "%s__%s__%s" % (tissue, stage, sex)
#                        print gene_rev
#                        print self.pd_exp.loc[gene_rev].values.shape, self.pd_exp.loc[gene_rev].values
#                        print np_idx, np_idx.shape
                        INFO[tag] = list(self.pd_exp.loc[gene_rev].values[0][np_idx])
            
            self.posts_PGC_exp.insert_one(INFO)
        
        self.posts_PGC_exp.create_index([("gene", 1)])
            

class DbTF(object):
    """docstring for DbTF"""
    def __init__(self):
        super(DbTF, self).__init__()
    
    def load_from_file(self, 
            infile_Freq="/data/Analysis/huboqiang/software/encode-motifs-v1.3/motifs.txt", 
            infile_Bed="/data/Analysis/huboqiang/software/encode-motifs-v1.3/matches.txt"
        ):
        logging.info("Initiating Frequency file from %s" % (infile_Freq))
        self.__load_freq(infile_Freq)
        logging.info("Frequency file initiation done!")
        
        logging.info("Initiating interval file from %s" % (infile_Bed))
        self.__load_interval(infile_Bed)
        logging.info("Initiating interval file done!")
    
    def __upload_TF(self):
        if 'TF_seq' in self.INFO:
            self.INFO['TF_seq'] = "".join(self.INFO['TF_seq'])
            self.posts_TF.insert_one(self.INFO)
    
    def __load_freq(self, infile_Freq):
        f_infile_Freq = open(infile_Freq, "r")
        self.INFO = {}
        for line in f_infile_Freq:
            line = line.strip()
            f = line.split()
            if line[0] == ">":
                self.__upload_TF()
                self.INFO = { 'TF_name':f[0][1:], 'TF_info': f[1], 'TF_seq':[], 'TF_freq':[], 'TF_cnt':0, 'prob_NDR':[] }
            else:
                self.INFO['TF_seq'].append(f[0])
                self.INFO['TF_freq'].append(f[1:])
        
        # the last line
        self.__upload_TF()
        f_infile_Freq.close()
        
        logging.info("Creating index for Transcriptional factor database.")
        self.posts_TF.create_index([("TF_name", 1)])
        logging.info("Creating index for Transcriptional factor database done.")
    
    def __load_interval(self, infile_Bed):
        M_TF_cnt = {}
        logging.info("Reading interval file %s ..." % (infile_Bed))
        f_infile_Bed = gzip.open(infile_Bed, "rb")
        for line in f_infile_Bed:
            line = line.strip()
            f = line.split()
            if f[0] not in M_TF_cnt:
                M_TF_cnt[f[0]] = 0
                
            M_TF_cnt[f[0]] += 1
        
        f_infile_Bed.close()
        
        logging.info("Reading interval file done.")
        logging.info("Uploading the frequency information to database.")
        for TF_name, TF_cnt in M_TF_cnt.items():
            self.posts_TF.update({'TF_name': TF_name}, {"$set" :{'TF_cnt': TF_cnt}})

class NDR_info(object):
    """Loading NDR information."""
    def __init__(self):
        super(NDR_info, self).__init__()
    
    def init_db(self, 
            gene_file = "/datd/huboqiang/test_hESC/database/refGene.up2000_down2000.promoter.Bsorted.longestTid.bed",
            motifBed_file = "/data/Analysis/huboqiang/software/encode-motifs-v1.3/matches.txt.gz"
        ):
        """ reference file used. """
        self.file_geneTSS_tb = tabix.open(gene_file)
        self.file_motifBed_tb = tabix.open(motifBed_file)
    
    def read_bed_file(self, infile_prefix):
        """ read distal region and TSS region respectly."""
        infile_distal = "%s.distal.merge.bed" % (infile_prefix)
        infile_TSS = "%s.TSS.merge.bed" % (infile_prefix)

        with open(infile_distal, "r") as f_infile:
            for line in f_infile:
                line = line.strip()
                self.__parse_line(line, "d")

        with open(infile_TSS, "r") as f_infile:
            for line in f_infile:
                line = line.strip()
                self.__parse_line(line, "t")
        
        self.posts_NDR_reg.create_index([("region", 1)])
    
    def __parse_line(self, line, ltype = "d"):
        f = line.split()
        chrom = f[0]
        if chrom != "chrM" and chrom != "lambda":
            beg = int(f[1])
            end = int(f[2])
            region = "%s:%s-%s" % (chrom, beg, end)
            print region
            try:
                len(self.posts_NDR_reg.find_one({"region":region}))
                
            except:
                l_sam_idx = np.array(f[3].split(","), dtype="string")
                sam_idx = (l_sam_idx == "1")
                np_sam = np.array(self.samInfo_pd_RNA['brief_name'], dtype="string")
                self.info = {
                    "region" : region, 
                    "type" : ltype,
                    "TF" : [],
                    "TF_id" : {},
                    "gene" : "",
                    "gene_100k" : {},
                    "sample" : list(np_sam[sam_idx])
                }
                ### list example:
                ### blogpost1.tags = ['tag1', 'tag2', 'tag3', 'tag4', 'tag5']
                ### db.blogpost.find({ 'tags' : 'tag1'})
                
                ### hash example:
                ### {foo : 1, bar : 4, baz : {a : 1, b : 2 ,c : "fafofu"}}
                ### db.my_collection.find({"baz.a" : 1, "baz.b" : 2})
                
                self.__NDR_gene_overlap(chrom, beg, end, ltype)
                self.__NDR_TF_overlap(chrom, beg, end)
#                print self.info
                self.posts_NDR_reg.insert_one(self.info)
    
    def __NDR_gene_overlap(self, chrom, beg, end, ltype):
        """ 
            For distal region, searching the 100kbp region up/down stream. For TSS region, searching 2kb around.
            As distal region shares no overlap with TSS region in previous setting, the two region would not overlap.
            Get the gene name, and their key in self.posts_PGC_exp
        """
        def score_onlyDist(input_beg, input_end, rec_beg, rec_end):
            return abs( int( (input_beg+input_end)/2-(rec_beg+rec_end)/2 ) )
        
        if ltype == "d":
            input_beg = beg - 100000
            input_end = end + 100000
        elif ltype == "t":
            input_beg = beg - 1000
            input_end = end + 1000
        
        if input_beg <= 0:
            input_beg = 1
        
        try:
            record = self.file_geneTSS_tb.query(chrom, input_beg, input_end)
            M_score_best = {}
            best_score = 100000
            for rec in record:
                rec_beg = int(rec[1])
                rec_end = int(rec[2])
                gene = rec[6]
                gene = gene.replace(".", "_")
#                gene_id_find = self.posts_PGC_exp.find_one({"gene":gene})
#                print gene_id_find, gene
                score = score_onlyDist(beg, end, rec_beg, rec_end)
                self.info["gene_100k"][gene] = score
                if score < best_score:
                    M_score_best['gene'] = gene
                    best_score = score
#                    print gene, score, best_score, beg, end, rec_beg, rec_end
            
            if "gene" in M_score_best:
                self.info['gene'] = M_score_best['gene']
#                print self.info['gene']

        except:
            pass
#            logging("No query in Genome TSS database for %s:%d-%d" % (chrom, input_beg, input_end))
            
    def __NDR_TF_overlap(self, chrom, beg, end):
        try:
            record = self.file_motifBed_tb.query(chrom, beg, end)
            for rec in record:
                TF_name = rec[4]
                TF_id_find = self.posts_TF.find_one({"TF_name":TF_name})
                if not (TF_id_find is None):
                    TF_id = TF_id_find['_id']
                    self.info['TF'].append(TF_name)
                    self.info["TF_id"][TF_name] = TF_id
        
        except:
            pass
#            logging("No query in TF database for %s:%d-%d" % (chrom, beg, end))
    
    
        
    