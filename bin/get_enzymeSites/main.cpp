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

void load_fai(string &file_name, map<string,uint64_t> &ChrLen){
   string file_fai = file_name + ".fai";
   ifstream infile;
   infile.open(file_fai.c_str());
   if ( ! infile ){
      cerr << "Cannot find index file" << file_fai << endl;
      exit(0);
   }
	
   string lineStr;
	while (getline(infile,lineStr,'\n')){
		if (lineStr[0] == ' ' || lineStr[0] == '\n'){
			continue;
		}
		vector<string> lineVec;
		boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
		ChrLen[lineVec[0]] = boost::lexical_cast<uint64_t>(lineVec[1]);
	}
	infile.close();
}

void usage()
{
   cout << "./get_enzymeSite -s CCGG.AGCT /data/Analysis/huboqiang/database/mm9_tophat/genome.fa" << endl;
   cout << "  -s  Sequence for cutting site of restrict enzyme. default = CCGG.AGCT"   << endl;
   cout << "  -h  get help information"   << endl;
   exit (0);
}


string l_site = "CCGG.AGCT";
int main(int argc, char *argv[])
{
   int c;
   while ( (c=getopt(argc,argv,"s:h")) != -1 ){
      switch(c)
      {
         case 's' : l_site = optarg;break;
         case 'h' : usage();break;
         default : usage();
      }
   }
   if (argc < 2) usage();
   string file_name  	= argv[optind++];	
   
   transform(l_site.begin(), l_site.end(),l_site.begin(), ::toupper);
   
   /* Files for enzyme Sites information */
	string file_out_bed = file_name + "." + l_site + ".bed";
	ofstream outfile_bed;
   outfile_bed.open(file_out_bed.c_str());
   if ( ! outfile_bed ){
      cerr << "fail to open input file" << outfile_bed << endl;
      exit(0);
   }
	
	string file_out_dist = file_name + "." + l_site + ".dist";
	ofstream outfile_dist;
   outfile_dist.open(file_out_dist.c_str());
   if ( ! outfile_dist ){
      cerr << "fail to open input file" << file_out_dist << endl;
      exit(0);
   }
   
   /* Put enzyme sites into vector */
	vector<string> V_site;
	boost::split( V_site, l_site,boost::is_any_of("."), boost::token_compress_on);
   int max_len = 0;
   for (int i=0;i<V_site.size();++i){
      if ( V_site[i].length()>max_len ){
         max_len = V_site[i].length();
      }
   }
   
   /* Put Chromsome length into map */
   map<string,uint64_t> 			  ChrLen;
	map<string,uint64_t >::iterator map_size;
   load_fai( file_name, ChrLen);
	
   /* Load Chromosome sequence */
   FastaReference *fr = new FastaReference;
   bool memmap = true;
   fr->open(file_name, memmap);
	
   /* 
      Iterate each chromosome to find enzyme sites.
      For a given site in chromosome chr,
        scan along the chromosome for position j,
        and get get cis-trans sequence for length of enzyme-sites.
   */

	for (map_size=ChrLen.begin();map_size!=ChrLen.end();map_size++){  
		string   chr = map_size->first;
	   uint64_t len = map_size->second;
		int pre_site = 0;
		for (int i=0; i<len-max_len; i++){
         string tmp = fr->getSubSequence( chr, i, max_len);
         for (int j=0; j<V_site.size();j++){
            string tmp1 = tmp.substr( 0,V_site[j].size() );
            transform(tmp1.begin(), tmp1.end(),tmp1.begin(), ::toupper);
            string tmp2 = reverse_complement( tmp1 );

   			if (tmp1 == V_site[j]){
   				outfile_bed << chr << "\t" << i+1 << "\t" << i+tmp1.size() << "\t" << "+" << "\t" << tmp1 << "\t" ;
   				if (pre_site != 0){
                  outfile_bed  << i - pre_site << endl;
   					outfile_dist << i - pre_site << endl;
   				}
               else{
                  outfile_bed  << 0 << endl;
               }
   				pre_site = i;
   			}
            
			if (tmp2 == V_site[j]){
				outfile_bed << chr << "\t" << i+1 << "\t" << i+tmp2.size() << "\t" << "-" << "\t" << tmp2 << "\t" << i - pre_site << endl;
				if (pre_site != 0){
					outfile_dist << i - pre_site << endl;
				}
				pre_site = i;
			}
            
         }
      }
   }
   
   outfile_bed.close();
   outfile_dist.close();
   
}
