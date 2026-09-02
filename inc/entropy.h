#include "common.h"

#ifndef _ENTROPY_H_
#define _ENTROPY_H_

#define NSGRP 8
using namespace std;

inline uint32_t getRowCmpSize(vector<uint8_t>& rRow){
    uint32_t rsza[4];
    for(int k=0; k < 4; ++k){
        rsza[k] = rRow[k];
    }
    uint32_t rsz = rsza[0] | (rsza[1]<<8)|(rsza[2]<<16)|(rsza[3]<<24);
    rsz -= 4;
    return rsz;
}
void enc_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info, 
                 const ProcessInfo& proc_info, const QuantInfo& quant_info );
void dec_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info,
                 const ProcessInfo& proc_info, const QuantInfo& quant_info );

#endif