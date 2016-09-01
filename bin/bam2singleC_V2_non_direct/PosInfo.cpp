#include "PosInfo.h"

PosInfo::PosInfo(string out_file_singleC){    
    info_pos = 0;
    info_umt = 0;
    info_met = 0;
    
    f_outfile.open(out_file_singleC.c_str());
    if ( ! f_outfile ){
        cerr << "Cannot find interval file" << out_file_singleC << endl;
        exit(0);
    }


    vector<string> tmp_vec;
    vector<string> tmp_vec2;
    boost::split(tmp_vec, out_file_singleC, boost::is_any_of("/"));
    boost::split(tmp_vec2, tmp_vec[ tmp_vec.size()-1 ], boost::is_any_of("."));
    info_chrom = tmp_vec2[0];
}

PosInfo::~PosInfo(){
    report_pos();
    f_outfile.close();
}

void PosInfo::parse_pos(uint64_t pos, uint64_t pos2, int umt=0, int  met=0){
    if (info_pos == pos){
        if(pos == pos2){
            info_umt = umt;
            info_met = met;
        }
    }
    else{
        report_pos();
        info_umt = 0;
        info_met = 0;
        if(pos == pos2){
            info_umt = umt;
            info_met = met;
        }   
    }
    info_pos = pos;
}

void PosInfo::report_pos(){
    if (info_pos != 0){
        f_outfile << info_chrom << "\t" << info_pos << "\t" << info_pos << "\t" << info_umt << "\t" << info_met << endl;
    }
}

