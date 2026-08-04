#include "common.h"
#include <vector>

#ifndef __MSST_H__
#define __MSST_H__
using namespace std;

void fmsst(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info);
void imsst(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
    const BayerInfo& bayer_info);
void fycc(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info);
void iycc(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
    const BayerInfo& bayer_info);
#endif
