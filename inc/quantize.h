#include "common.h"

#ifndef _QUNATIZE_H_
#define _QUNATIZE_H_

using namespace std;

class QUANT{
public:
    QuantInfo quant_info;
    MsstInfo msst_info;
    int32_t  dwt_lv;
    QUANT(){}
    QUANT(const QuantInfo& _quant_info, const MsstInfo& _msst_info, const int32_t _dwt_lv):
        quant_info(_quant_info), msst_info(_msst_info), dwt_lv(_dwt_lv){}
    ~QUANT(){}
    void fquant_msstv(vector<msstSm>& msstv);
    void iquant_msstv(vector<msstSm>& msstv);
};

#endif