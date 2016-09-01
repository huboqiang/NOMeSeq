#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "split.h"
#include "Fasta.h"

#ifndef POSINFO_H_
#define POSINFO_H_

class PosInfo{
    public:
        uint64_t info_pos;
        int info_umt;
        int info_met;
        std::string info_chrom;
        
        ofstream f_outfile;
        
        PosInfo(std::string chrom);
        ~PosInfo();
        void parse_pos(uint64_t pos, uint64_t pos2, int umt, int  met);
        void report_pos();
};
#endif