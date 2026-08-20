#include "bayer_dwt.h"
#include "monitor.h"
static inline void sat(vector<cfapix>& dst, const uint32_t mxv){
    for(size_t k=0; k < dst.size(); ++k){
        if(dst[k] < 0){
            dst[k] = 0;
        }else if(dst[k] > mxv){
            dst[k] = mxv;
        }
    }
}
void bayer_fdwt(vector<msstSm>& dst, vector<cfapix>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info)
{

    
    MsstInfo msst_info;
    msst_info.w = bayer_info.w >> 1;
    msst_info.h = bayer_info.h >> 1;
    vector<msst> msstv( msst_info.w * msst_info.h, msst());

    //-- Perform NLT on the source Bayer image if gamma is not 1.0
    int32_t mxv = (1U << bayer_info.bpp) - 1;
    if(proc_info.g != 1.0f){
        nlt(src, mxv, proc_info.g);
    }

    //fprintf(stdout, "bayer_val[888962] = %d\n", src[888962]);
    //-- Perform forward MSST on the Bayer image
    int32_t dp = 0, sp = 0;
    //if( proc_info.sel_ycc ){
    //    for(int32_t h = 0; h < bayer_info.h; h+=2){
    //        fycc(msstv, dp, src, sp, bayer_info);
    //        sp += bayer_info.w*2;
    //    }
    //}else{
        for(int32_t h = 0; h < bayer_info.h; h+=2){
            fmsst(msstv, dp, src, sp, bayer_info, proc_info);
            sp += bayer_info.w*2;
        }
    //}
    cout << "bayer_fdwt, fmsst: dp = " << dp << ", sp = " << sp << endl;

    //-- Perform forward DWT on the MSST data
    if(msstv.size() != msst_info.w * msst_info.h){
        msstv.resize(msst_info.w * msst_info.h, msst());
    }
    int32_t S, N;
    for(int32_t lv =0 ; lv < proc_info.dwt_lv; lv++){
        N = msst_info.w >> lv;
        for(int32_t h = 0; h < msst_info.h; h++){
            S = h * msst_info.w;
            fdwt53(msstv, S, N);
        }
    }
    
    //dumpmsst(msstv,"msstv.txt", msst_info);
    //-- Converting msst to msstSm
    ConvMsst2MsstSm(dst, msstv);
    
    //-- Quantization
    quant_msstv(dst, msst_info, proc_info.dwt_lv, quant_info);

}

void bayer_idwt(vector<cfapix>& dst, vector<msstSm>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info)
{
    MsstInfo msst_info;
    msst_info.w = bayer_info.w >> 1;
    msst_info.h = bayer_info.h >> 1;
    vector<msst> msstv(msst_info.w * msst_info.h, msst());

    if(src.size() != msst_info.w * msst_info.h){
        cerr << "Error: Source msst size does not match expected dimensions." << endl;
        return;
    }

    //-- Inverse Quantization
    iquant_msstv(src, msst_info, proc_info.dwt_lv, quant_info);

    //-- convert msstSm to msst
    ConvMsstSm2Msst( msstv, src);

    //-- Perform inverse DWT on the msst data
    int32_t S, N;
    for(int32_t lv = proc_info.dwt_lv-1; lv >= 0; --lv){
        N = msst_info.w >> lv;
        for(int32_t h = 0; h < msst_info.h; h++){
            S = h * msst_info.w;
            idwt53(msstv, S, N);
        }
    }
    //-- Perform inverse MSST on the DWT data
    int32_t dp = 0, sp = 0;
    //if( proc_info.sel_ycc ){
    //    for(int32_t h = 0; h < bayer_info.h; h+=2){
    //        iycc(dst, dp, msstv, sp, bayer_info);
    //        dp += bayer_info.w*2;
    //    }
    //}else{
        for(int32_t h = 0; h < bayer_info.h; h+=2){
            imsst(dst, dp, msstv, sp, bayer_info, proc_info);
            dp += bayer_info.w*2;
        }
    //}
    uint32_t mxv = (1U << bayer_info.bpp) - 1;
    sat(dst, mxv);
    //fprintf(stdout, "bayer_val[888962] = %d\n", dst[888962]);
    cout << "bayer_idwt, imsst: dp = " << dp << ", sp = " << sp << endl;
    //-- Perform inverse NLT on the reconstructed Bayer image if gamma is not 1.0
    //int32_t mxv = (1U << bayer_info.bpp) - 1;
    if(proc_info.g != 1.0f){
        inlt(dst, mxv, proc_info.g);
    }

}
