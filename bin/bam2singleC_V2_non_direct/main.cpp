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
#include "PosInfo.h"

using namespace std;

char alphabet[128] =
{
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

string reverse_complement(string &in_str)
{
    char c_bases[5] ={
        'T', 'G', 'C', 'A', 'N'
    };
    string out_str;
    for (int i=in_str.size()-1; i>=0; i--)
    {   
        out_str.push_back(c_bases[alphabet[in_str[i]]]);
    }
    return out_str;
}

void usage()
{
    /*
        awk '{OFS="\t";print $1,$2,$2+1}' /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG/chr10.bed | samtools mpileup  -O  -f /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa -l /dev/stdin /datd/huboqiang/test_NOM/02.SingleC/mESC_gWBS1_1/bam/mESC_gWBS1_1.sort.rmdup.chr10.bam  | awk '{OFS="\t";print $1,$2,$2,$3,$4,$5,$6,$7}' |  bedtools intersect -sorted -loj -a /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG/chr10.rev.bed -b /dev/stdin | le
    */
    
    
   cout << "awk '{OFS=\"\\t\";print $1,$2,$2+1}' /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG/chr10.bed | samtools mpileup  -O  -f /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa -l /dev/stdin /datd/huboqiang/test_NOM/02.SingleC/mESC_gWBS1_1/bam/mESC_gWBS1_1.sort.rmdup.chr10.bam  | awk '{OFS=\"\\t\";print $1,$2,$2,$3,$4,$5,$6,$7}' |  bedtools intersect -sorted -loj -a /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG/chr10.rev.bed -b /dev/stdin | ./pileup2singleC /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa /dev/stdin >/datd/huboqiang/test_NOM/02.SingleC/mESC_gWBS1_1/singleC/chr10.ACG.TCG.bed " << endl;
   cout << "  -h  get help information"   << endl;
   exit (0);
}


int main(int argc, char *argv[])
{
    int c;
    int d = 5;
    while ( (c=getopt(argc,argv,"d:h")) != -1 ){
     switch(c)
     {
       case 'h' : usage();break;
       case 'd' : d=atoi(optarg);break;
       default : usage();
     }
    }
    if (argc < 4) usage();
    string file_name = argv[optind++];   
    string file_pileup = argv[optind++];
    string out_file_singleC= argv[optind++];
    
//    string out_file_totCnt = out_file_singleC + ".CT_total.xls";
    
    
    FastaReference *fr = new FastaReference;
    bool memmap = true;
    fr->open(file_name, memmap);
     
     
    ifstream infile;
    infile.open(file_pileup.c_str() );
    if (!infile){
        cerr << "fail to open input file" << file_pileup << endl;
        exit(0);
    }
    
//    ofstream f_outfile_tot;
//    f_outfile_tot.open(out_file_totCnt.c_str());
//    if ( ! f_outfile_tot ){
//        cerr << "Cannot find interval file" << out_file_totCnt << endl;
//        exit(0);
//    }
    
    
    string lineStr;

    uint64_t cnt_umt_tot = 0;
    uint64_t cnt_met_tot = 0;
    uint64_t pos = 0, pos2 = 0, idx_max = 0;
    string  chr, ref, base;
    int     dep, met, umt;
    
    PosInfo pos_bed(out_file_singleC);
    
    while ( getline(infile,lineStr,'\n') ){
        if (lineStr[0] == ' ' || lineStr[0] == '\n' ){
            continue;
        }
        vector<string> lineVec;
        boost::split(lineVec,lineStr, boost::is_any_of(" \t\n"), boost::token_compress_on);
        
        chr = lineVec[0];
        pos = boost::lexical_cast<uint64_t>(lineVec[1]);
        pos2= boost::lexical_cast<uint64_t>(lineVec[4]);
        
        if (lineVec[3] == "."){
            pos_bed.parse_pos(pos, 0, 0, 0);
        }
        else{
            ref = lineVec[6];     
            dep = boost::lexical_cast<uint64_t>(lineVec[7]);
        
            boost::to_upper(ref);
        
            size_t seqLength = fr->sequenceLength(chr);
        
            base = lineVec[8];
        
            uint64_t beg = pos - 1;
            if ( ref == "C" ){
                met = count( base.begin(),base.end(),'.' );     // 正链 C
                umt = count( base.begin(),base.end(),'T' ) + count( base.begin(),base.end(),'t' );     // 正负链 T
            }
            else{
                met = count( base.begin(),base.end(),',' );     // 负链 G
                umt = count( base.begin(),base.end(),'A' ) + count( base.begin(),base.end(),'a' );     // 正负链 A
            }
            pos_bed.parse_pos(pos, pos2, umt, met);
            if (met + umt >= d){
                cnt_umt_tot += umt;
                cnt_met_tot += met;
            }
            idx_max = pos;
        }
    }

//    f_outfile_tot << cnt_umt_tot << "\t" << cnt_met_tot << endl;
    infile.close();
}
