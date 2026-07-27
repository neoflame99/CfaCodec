#include "common.h"

#ifndef _QUNATIZE_H_
#define _QUNATIZE_H_

using namespace std;

void  quant_msstv(vector<msstSm>& msstv, const MsstInfo& msst_info, const int32_t dwt_lv, const QuantInfo& quant_info);
void iquant_msstv(vector<msstSm>& msstv, const MsstInfo& msst_info, const int32_t dwt_lv, const QuantInfo& quant_info);

#endif