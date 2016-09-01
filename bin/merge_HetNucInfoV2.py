#-*- coding:utf-8 -*-
from __future__ import division
import sys

import numpy as np
import tabix
from optparse   import OptionParser


def prepare_optparser():
    usage ="""usage: %s [options] input.bed

Merge a bed file with the information(value) in a given column.
In this class, the interval of the bed file could be merged based on
 the information of:
 1. Seed region. (level 1)
 2. Extend region. (level 2)

Seed region required a higher value, while the Extend region were based
on the Seed region. Extend region should be located between Seed regions
in order to connect seed regions.
        
Example:
    python %s -c 4 -cutoff_min_len 3000 --cutoff_max_avgVal -0.7 --max_ext_len 500 --max_avg_val 0.2 --min_seed_val 1 --min_ext_val 1.5  --maxVal 15.92  /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hheart_Week5_female_R1/singleC/all.GCA.GCC.GCT.GlobalPvalue.bed

    """ % (sys.argv[0],sys.argv[0])

    description = "Merge singleC results."

    optparser = OptionParser(
        version="%s v0.1 2016.1.15" % (sys.argv[0]),
        description=description,
        usage=usage,
        add_help_option=False
    )

    optparser.add_option(
        "-c", "--val_col", default=4,
        help="\nColumn for values for merging cutoff. [default: %default]"
    )
    optparser.add_option(
        "--cutoff_min_len", default=3000,
        help="\nMin length for reporting. [default: %default]"
    )
    optparser.add_option(
        "--cutoff_max_avgVal", default=-0.7,
        help="\nMax average value for reporting. [default: %default]"
    )
    optparser.add_option(
        "--max_ext_len", default=500,
        help="\nMax length for Extending region between Seed regions. [default: %default]"
    )
    optparser.add_option(
        "--max_avg_val", default=0.2,
        help="\nMax average value for consider as a Seed+Extend region. [default: %default]"
    )
    optparser.add_option(
        "--min_seed_val", default=1,
        help="\nMin value for a given interval for considering as a Seed region. [default: %default]"
    )
    optparser.add_option(
        "--min_ext_val", default=1.5,
        help="\nMin value for a given interval for considering as a Extend region. [default: %default]"
    )
    optparser.add_option(
        "--maxVal", default=15.92,
        help="\nValue for replacing the inf. [default: %default]"
    )
    return optparser

func_format = lambda x: "%1.2f" % x


class LevelCondition(object):
    '''
        The condition class were the father class of the merge class.
        Condition for merge were defined here and can be modified.
    '''
    def __init__(self, **kwarg):
        self.cutoff_min_len = kwarg['cutoff_min_len']
        self.cutoff_max_avgVal = kwarg['cutoff_max_avgVal']
        self.max_ext_len = kwarg['max_ext_len']
        self.max_avg_val = kwarg['max_avg_val']
        self.min_seed_val = kwarg['min_seed_val']
        self.min_ext_val = kwarg['min_ext_val']
        self.maxVal = kwarg['maxVal']

    def cond_Seed_region(self, value):
        '''
            For Seed region, it means that 
            3. The value of this line is so low that below the min_seed_val
        '''
        pass_idx = 0
        if (value < self.min_seed_val) and ~np.isnan(value):
            pass_idx = 1
        return pass_idx

    def cond_Ext_region(self, end_pass, end, l_value, Chi_global_p):
        '''
            For Extending region, it means that 
            1. The end position for this region from the begginning position
               should be within max_ext_len,
            2. The average value should below the max_avg_val.
            3. The value of this line should below the min_ext_val, but higher
               than min_seed_val of course.
        '''
        delta_pass = end - end_pass
        pass_idx = 0
        np_chi_global = np.array(l_value + [Chi_global_p], dtype="float")
        if (end_pass != 0) and (delta_pass < self.max_ext_len)              and\
            (np_chi_global.mean() < self.max_avg_val) and (Chi_global_p < self.min_ext_val):
            pass_idx = 1
        return pass_idx

    def cond_Elong_region(self, beg, end_pass):
        '''
            For elongation, the begin position should be within the max_ext len
            from the end_pass position.
        '''
        pass_idx = 0
        if beg < end_pass + self.max_ext_len:
            pass_idx = 1
        return pass_idx
    
    def cond_Print_region(self, length, mean_val):
        pass_idx = 0
        if (length > self.cutoff_min_len) and (mean_val < self.cutoff_max_avgVal):
            pass_idx = 1
        return pass_idx

    def read_value(self, in_str):
        Chi_global_p = float(in_str)
        if in_str == "inf":
            Chi_global_p = self.maxVal
        elif in_str == "-inf":
            Chi_global_p = -self.maxVal
        return Chi_global_p

class TwoLeverBedMerge(LevelCondition):
    '''
        Merge a bed file with the information(value) in a given column.
        In this class, the interval of the bed file could be merged based on
         the information of:
         1. Seed region. (level 1)
         2. Extend region. (level 2)
        
        Seed region required a higher value, while the Extend region were based
        on the Seed region. Extend region should be located between Seed regions
        in order to connect seed regions.
    '''
    def __init__(self, infile, val_col, **kwarg):
        
        super(TwoLeverBedMerge, self).__init__(**kwarg)
        self.f_input = open(infile, "r")
        self.val_col = val_col - 1
        self.M_interval = {}
    
    def process_file(self):
        self.__Delete()
#        print self.M_interval
        for line in self.f_input:
            line = line.strip()
            f = line.split()
            beg = int(f[1])
            end = int(f[2])
            
            Chi_global_p = self.read_value(f[self.val_col])
            
#            Chi_global_p = float(f[self.val_col])
#            if f[self.val_col] == "inf":
#                Chi_global_p = maxVal
#            elif f[self.val_col] == "-inf":
#                Chi_global_p = -maxVal

            pass_idx = 0
            end_pass = self.M_interval['end_pass']
            l_value = self.M_interval["Chi_global_p"]
            
            if self.cond_Seed_region(Chi_global_p):
                pass_idx = 1
            elif self.cond_Ext_region(end_pass, end, l_value, Chi_global_p):
                pass_idx = 2
            
            if (pass_idx == 1) or (pass_idx == 2):
                if len(self.M_interval['chr']) == 0:
                    if pass_idx == 1:
                        self.__Init(f)
                else:
                    if self.cond_Elong_region(beg, end_pass):
                        self.__Elong(f, pass_idx)
                    else:
                        self.__Output()
                        self.__Delete()
                        self.__Init(f)
            else:
                if len(self.M_interval['chr']) > 0:
                    self.__Output()
                    self.__Delete()

        if len(self.M_interval['chr']) > 0:
            self.__Output()
            self.__Delete()

        self.f_input.close()

    def __Delete(self):
        self.M_interval = {
            "chr": "",
            "beg": 0,
            "end": 0,
            "end_pass": 0,
            "Chi_global_p" : [],
        }
    
    def __Init(self, f):
        beg = int(f[1])
        end = int(f[2])
        self.M_interval = {
            "chr": f[0],
            "beg": beg,
            "end": end,
            "end_pass": end,
            "Chi_global_p" : [ float(f[self.val_col]) ],
            "Chi_global_p_pass" : [ float(f[self.val_col]) ],
        }
    
    def __Elong(self, f, pass_idx):
        beg = int(f[1])
        end = int(f[2])
        self.M_interval['end'] = end
        self.M_interval["Chi_global_p"].append(float(f[self.val_col]))
        if pass_idx == 1:
            self.M_interval['end_pass'] = self.M_interval['end']
            self.M_interval["Chi_global_p_pass"] = self.M_interval["Chi_global_p"]
    
    def __Output(self):
        l_Chi_global_p  = ",".join(map(func_format, self.M_interval["Chi_global_p"]))
        out_str = "\t".join([l_Chi_global_p])
        length = self.M_interval['end_pass'] - self.M_interval['beg']
        np_rec = np.array(self.M_interval["Chi_global_p"], dtype="float")
        mean_val = np_rec.mean()
        if self.cond_Print_region(length, mean_val):
            print "%s\t%d\t%d\t%d\t%f\t%s" % (self.M_interval['chr'], self.M_interval['beg'], self.M_interval['end'], length, np_rec.mean(), out_str)



def main():
    prepare_optparser()
    (options,args) = prepare_optparser().parse_args()
    try:
        file_input = args[0]
        val_col = int(options.val_col)
        cutoff_min_len = int(options.cutoff_min_len)
        cutoff_max_avgVal = float(options.cutoff_max_avgVal)
        max_ext_len = int(options.max_ext_len)
        max_avg_val = float(options.max_avg_val)
        min_seed_val = float(options.min_seed_val)
        min_ext_val = float(options.min_ext_val)
        maxVal = float(options.maxVal)
    except IndexError:
        prepare_optparser().print_help()
        sys.exit(1)
    
    M_info = {
        'cutoff_min_len' : cutoff_min_len,
        'cutoff_max_avgVal' : cutoff_max_avgVal,
        'max_ext_len' : max_ext_len,
        'max_avg_val' : max_avg_val,
        'min_seed_val' : min_seed_val,
        'min_ext_val' : min_ext_val,
        'maxVal' : maxVal,
    }
    
    m_mrg = TwoLeverBedMerge(file_input, val_col, **M_info)
    m_mrg.process_file()

if __name__ == '__main__':
    main()
