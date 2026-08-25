
#ifndef _NLT_H_
#define _NLT_H_

#include "common.h" 
using namespace std;

class NLT{
public:
    int32_t mxv;
    float gm;
    NLT(int32_t _mxv, float _gm):mxv(_mxv), gm(_gm){}
    ~NLT(){}
    void fnlt(vector<cfapix>& vcfa);
    void inlt(vector<cfapix>& vcfa);
};
inline void NLT::fnlt(vector<cfapix>& vcfa){
    if(gm == 1.0f){
        return;
    }
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

inline void NLT::inlt(vector<cfapix>& vcfa){
    if(gm == 1.0f){
        return;
    }
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