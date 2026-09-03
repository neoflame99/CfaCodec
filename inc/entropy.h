#include "common.h"

#ifndef _ENTROPY_H_
#define _ENTROPY_H_

#define NSGRP 8
#define SZSZ  2  // byte size of a row bitstream
#define KBSZ  4  // bit size of K to store
using namespace std;

inline uint32_t getRowCmpSize(vector<uint8_t>& rRow){
    uint32_t rsza[SZSZ];
    for(int k=0; k < SZSZ; ++k){
        rsza[k] = rRow[k];
    }    
    uint32_t rsz = 0;
    for(int k=0; k < SZSZ; ++k){
         rsz |= rsza[k]<<(8*k);
    }
    //rsz -= 2;
    return rsz;
}
void enc_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info, 
                 const ProcessInfo& proc_info, const QuantInfo& quant_info );
void dec_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info,
                 const ProcessInfo& proc_info, const QuantInfo& quant_info );
void enc_entropy_gr(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info );
void dec_entropy_gr(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info );
#endif
