
#ifndef _NLT_H_
#define _NLT_H_

#include "common.h" 
using namespace std;

inline void nlt(vector<cfapix>& vcfa, int32_t mxv, float gm){
    double np;
    int32_t w = vcfa.size();
    for(int32_t k=0; k < w; ++k){
        np = (double)vcfa[k]; 
        np /= mxv; 
        np = pow( np, gm);
        np *= mxv; 
        vcfa[k] = (cfapix) np;
    }
}

inline void inlt(vector<cfapix>& vcfa, int32_t mxv, float gm){
    double np;
    double igm = 1.0/gm;
    int32_t w = vcfa.size();
    for(int32_t k=0; k < w; ++k){
        np = (double)vcfa[k]; 
        np /= mxv; 
        np = pow( np, igm);
        np *= mxv; 
        vcfa[k] = (cfapix) np;
    }
}
#endif