#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <zlib.h>
#include <pthread.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "gzstream.h"

using namespace std;

int depth = 5;
int threadNum = 8;
int tmp_i = 0;

string file_prefix;
string cut_site;
string file_out;
ofstream outfile;


uint64_t bufferNum = 100000;
uint64_t readNum;
int *job_processing;

string *inBlock;
string *outBlock;
int **Array_umt_meth;

vector<string> V_chrom;


void usage() {
    cout << "/home/huboqiang/module_tangpipe/MethGC/bin/depth_stat/depth_stat -p 8 -d 5 /data/Analysis/huboqiang/Database_Meth/hg19/hg19_lambda.fa.fai /date/huboqiang/NORM_seq_PGC/human/02.SingleC/hPGC_Week10_rep1/singleC GCA.GCC.GCT" << endl;
    cout << "   -p <int>   threads number, default=" << threadNum << endl;
    cout << "   -d <int>   depth, default=" << depth << endl;
    cout << "   -h         get help information" << endl;
    exit(0);
}

void loadRef(string infai, vector<string> &V_chrom){
    ifstream infile;
    infile.open(infai.c_str());
    if ( ! infile ){
        cerr << "Cannot find interval file" << infai << endl;
        exit(0);
    }
    string lineStr;
    
    while(getline(infile, lineStr, '\n')){
        vector<string> lineVec;
        boost::split(lineVec, lineStr,boost::is_any_of("\t "), boost::token_compress_on);
        V_chrom.push_back(lineVec[0]);
    }
    infile.close();
}

void *thread_parse_line(void* threadId_p){
    int threadId = *((int*)threadId_p);
    while (1){
        usleep(1);
        if (job_processing[threadId] == 1){
            for (uint64_t i=0; i<readNum; i++){
                if (i%threadNum == threadId){
                    vector<string> lineVec;
                    outBlock[i] = inBlock[i];
                    boost::split(lineVec, inBlock[i],boost::is_any_of("\t"), boost::token_compress_on);

                    int umt = boost::lexical_cast<int>(lineVec[3]);
                    int met = boost::lexical_cast<int>(lineVec[4]);
                    if (umt+met >= depth){
                        Array_umt_meth[tmp_i][3*threadId+0] += umt;
                        Array_umt_meth[tmp_i][3*threadId+1] += met;
                        Array_umt_meth[tmp_i][3*threadId+2] += 1;
                    }
                }
            }
            job_processing[threadId] = 0;
        }
        else if (job_processing[threadId] == 2){
            return NULL;
        }
    }
}


void read_chrom_depth_phread(){
    string chrom = V_chrom[tmp_i];
    igzstream infile;
    string chrom_file = file_prefix + "/" + chrom + "." + cut_site + ".bed.gz";
    infile.open(chrom_file.c_str());
    if ( ! infile ){
        cerr << "Cannot find interval file " << chrom_file << endl;
        exit(0);
    }
    
    inBlock = new string[bufferNum];
    outBlock= new string[bufferNum];
    pthread_t *pthread = new pthread_t[threadNum];
    int *pthreadId = new int[threadNum];
    job_processing = new int[threadNum];
    for (int i=0; i<threadNum; i++){   
        job_processing[i] = 0;
        pthreadId[i] = i;
        pthread_create((pthread+i), NULL, thread_parse_line, (void*)(pthreadId+i));
    }
//    cerr << threadNum << " threads creation done!\n" << endl;
    
    while (1){
        //读取bufferNum行
        readNum = 0;
        while (readNum < bufferNum && getline(infile, inBlock[readNum], '\n')){
            readNum ++;
        }
//        cerr << "Reading done." << endl;
        //线程运行，直到Loop_id等于readNum, 然后进入子线程等待状态
        //放行主线程，进入下一个循环，读数据
        for (int i=0; i<threadNum; i++){
            job_processing[i] = 1;
        }

        //等子线程，直到job_processing都是0，说明任务完成，可以继续读下一批数据
        while (1){
            usleep(1);
            int i=0;
            for (; i<threadNum; i++){
                if (job_processing[i] == 1){
                    break;
                }
            }
            if (i == threadNum){
                break;
            }
        }

        //如果已经到了文件尾，则终止全部子线程，并且退出读文件的循环
        if (readNum < bufferNum){
            for (int i=0; i<threadNum; i++){   
                job_processing[i] = 2;
            }
            break;//读到文件尾退出
        }
    }
    
    //等待全部子线程结束

    for (int i=0; i<threadNum; i++){
        pthread_join(pthread[i], NULL);
    }
    delete[] pthread;
    delete[] pthreadId;
    delete[] job_processing;
    delete[] inBlock;
}

int main(int argc, char *argv[]){
    int c;

    while ( (c=getopt(argc,argv,"p:d:h")) != -1 ){
        switch(c){
            case 'p': threadNum = atoi(optarg); break;
            case 'd': depth = atoi(optarg); break;
            case 'h': usage(); break;
            default:  usage();
        }
    }
    if (argc < 3) usage();

    string infai = argv[optind++];

    file_prefix = argv[optind++];
    cut_site = argv[optind++];
    file_out = file_prefix + "/" + cut_site + ".log";
    outfile.open(file_out.c_str());
    if ( ! outfile ){
        cerr << "Cannot find interval file " << file_out << endl;
        exit(0);
    }

    loadRef(infai, V_chrom);
    
    Array_umt_meth = new int *[V_chrom.size()];

    uint64_t allumt_tot = 0;
    uint64_t allmet_tot = 0;
    uint64_t allcnt_tot = 0;

    for(int i=0; i<V_chrom.size(); ++i){
        tmp_i = i;
        Array_umt_meth[i] = new int[3*threadNum];
        for(int j=0; j<3*threadNum; ++j){
            Array_umt_meth[i][j] = 0;
        }

        read_chrom_depth_phread();
        uint64_t umt_tot = 0;
        uint64_t met_tot = 0;
        uint64_t cnt_tot = 0;
        for(int j=0; j<threadNum; ++j){
            umt_tot += Array_umt_meth[i][3*j];
            met_tot += Array_umt_meth[i][3*j+1];
            cnt_tot += Array_umt_meth[i][3*j+2];
        }
        outfile << V_chrom[i] << "\t" << umt_tot << "\t" << met_tot << "\t" << cnt_tot << endl;
        if(V_chrom[i].length()>3){
            string tmp = V_chrom[i].substr(0, 3);
            if(tmp == "chr"){
                if(V_chrom[i][3] != 'M'){
                    allumt_tot += umt_tot;
                    allmet_tot += met_tot;   
                    allcnt_tot += cnt_tot;   
                }
            }
        }
    }
    outfile << "all\t" << allumt_tot << "\t" << allmet_tot << "\t" << allcnt_tot << endl;
    outfile.close();
    return 1;
}