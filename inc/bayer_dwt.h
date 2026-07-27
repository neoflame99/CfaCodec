#include <string>
#include "common.h"
#include "nlt.h"
#include "msst.h"
#include "dwt.h"
#include "quantize.h"

#ifndef _BAYER_DWT_H_
#define _BAYER_DWT_H_

using namespace std;

void bayer_fdwt(vector<msstSm>& dst, vector<cfapix>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info);

void bayer_idwt(vector<cfapix>& dst, vector<msstSm>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info);

#endif //_BAYER_DWT_H_