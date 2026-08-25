
#include "common.h"

#ifndef __DWT_H__
#define __DWT_H__
using namespace std;

class DWT{
public:
    MsstInfo msst_info;
    ProcessInfo proc_info;
    DWT(){}
    DWT(const MsstInfo& _msst_info, const ProcessInfo& _proc_info): msst_info(_msst_info), proc_info(_proc_info){}
    ~DWT(){}
    void enc_dwt53(vector<msst>& dst );
    void dec_dwt53(vector<msst>& dst );
};

#endif