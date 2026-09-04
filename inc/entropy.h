#include "common.h"

#ifndef _ENTROPY_H_
#define _ENTROPY_H_

#define NSGRP 8
#define SZSZ  2  // byte size of a row bitstream
#define KBSZ  4  // bit size of K to store
using namespace std;

struct GR{
    uint32_t q;
    uint32_t r;
    GR():q(0), r(0){}
};

class Entropy{
public:
    MsstInfo    msst_info;
    ProcessInfo proc_info;
    Entropy(const MsstInfo& _msst_info, const ProcessInfo& _proc_info):
        msst_info(_msst_info), proc_info(_proc_info)
    {
    }

    ~Entropy(){}

    inline uint32_t getRowCmpSize(vector<uint8_t>& rRow){
        uint32_t rsza[SZSZ];
        for(int k=0; k < SZSZ; ++k){
            rsza[k] = rRow[k];
        }
        uint32_t rsz = 0;
        for(int k=0; k < SZSZ; ++k){
             rsz |= rsza[k]<<(8*k);
        }
        return rsz;
     }
    void enc_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;
    void dec_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;

    void enc_entropy_bc(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;
    void dec_entropy_bc(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;
    void enc_entropy_gr(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;
    void dec_entropy_gr(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv) ;
};

#endif
