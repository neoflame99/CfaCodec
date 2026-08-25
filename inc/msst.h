#include "common.h"
#include <vector>

#ifndef __MSST_H__
#define __MSST_H__
using namespace std;

class MSST{
public:
    BayerInfo bayer_info; 
    ProcessInfo proc_info;
    MSST(){}
    MSST(const BayerInfo& _bayer_info, const ProcessInfo& _proc_info):bayer_info(_bayer_info), proc_info(_proc_info){}
    ~MSST(){}
    void enc_msst(vector<msst>& dst, const vector<cfapix>& src);
    void dec_msst(vector<cfapix>& dst, const vector<msst>& src);
};
void fycc(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info);
void iycc(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
    const BayerInfo& bayer_info);
#endif
