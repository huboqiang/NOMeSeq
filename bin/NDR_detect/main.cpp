#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <zlib.h>
#include <pthread.h>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include "gzstream.h"
#include "NDR.h"

using namespace std;

int depth = 5;
double cutoff = 0.00001;
int bin = 100;
int step_len = 20;
int dist_len = 140;

string file_count = "/datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/singleC/chr10.GCA.GCC.GCT.bed.CT_total.xls";
void usage()
{   cout << "./NDR_detect -n /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/singleC/chr10.GCA.GCC.GCT.bed.CT_total.xls -d 5 /datd/huboqiang/test_NOM/02.SingleC/mESC_gWBS1_1/singleC/chr10.GCA.GCC.GCT.bed.gz " << endl;
    cout << "   -n <file>   file for CT_total count information, default = " << file_count << endl;
    cout << "   -d <int>    Minimun depth for a site to take into consideration, default = " << depth << endl;
    cout << "   -c <double> Threshold for passing Chisquare-test for a bin to take into consideration, default = " << cutoff << endl;
    cout << "   -b <int> length for a bin, default = " << bin << endl;
    cout << "   -s <int> length for sliding window step, default = " << step_len << endl;
    cout << "   -l <int> minimum length for a NDR region, default = " << dist_len << endl;
    cout << "   -h         get help information" << endl;
    exit(0);
}

int main(int argc, char *argv[]){
    while( (c=getopt(argc, argv, "n:d:c:b:s:l:h")) != -1 ){
        switch(c){
            case 'n': file_count=optarg; break;
            case 'd': depth=atod(optarg); break;
            case 'c': cutoff=(double)(optarg); break
            case 'b': bin=atoi(optarg); break;
            case 's': step_len=atoi(optarg); break;
            case 'l': dist_len=atoi(optarg); break
            case 'h': usage(); break;
            default: usage();
        }
    }
    
    if (argc < 2) usage();
    
    string bed_file = argv[optind++];
    
    
    
    
}