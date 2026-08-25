#include "bayer_dwt.h"
#include "monitor.h"
//static inline void sat(vector<cfapix>& dst, const uint32_t mxv){
//    for(size_t k=0; k < dst.size(); ++k){
//        if(dst[k] < 0){
//            dst[k] = 0;
//        }else if(dst[k] > mxv){
//            dst[k] = mxv;
//        }
//    }
//}
void bayer_fdwt(vector<msstSm>& dst, vector<cfapix>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info)
{

    MsstInfo msst_info;
    msst_info.w = bayer_info.w >> 1;
    msst_info.h = bayer_info.h >> 1;
    vector<msst> msstv( msst_info.w * msst_info.h, msst());

    int32_t mxv = (1U << bayer_info.bpp) - 1;
    NLT   nlt(mxv, proc_info.g);
    MSST  cmsst(bayer_info, proc_info);
    DWT   cdwt(msst_info, proc_info);
    QUANT cquant(quant_info, msst_info, proc_info.dwt_lv);

    //-- Perform NLT on the source Bayer image if gamma is not 1.0
    nlt.fnlt(src);

    //-- Perform forward MSST on the Bayer image
    cmsst.enc_msst(msstv, src);

    //-- Perform forward DWT on the MSST data
    if(msstv.size() != msst_info.w * msst_info.h){
        msstv.resize(msst_info.w * msst_info.h, msst());
    }
    cdwt.enc_dwt53(msstv);
    
    //dumpmsst(msstv,"msstv.txt", msst_info);
    //-- Converting msst to msstSm
    ConvMsst2MsstSm(dst, msstv);
    
    //-- Quantization
    cquant.fquant_msstv(dst);

}

void bayer_idwt(vector<cfapix>& dst, vector<msstSm>& src, const BayerInfo& bayer_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info)
{
    MsstInfo msst_info;
    msst_info.w = bayer_info.w >> 1;
    msst_info.h = bayer_info.h >> 1;
    vector<msst> msstv(msst_info.w * msst_info.h, msst());

    int32_t mxv = (1U << bayer_info.bpp) - 1;
    NLT   nlt(mxv, proc_info.g);
    MSST  cmsst(bayer_info, proc_info);
    DWT   cdwt(msst_info, proc_info);
    QUANT cquant(quant_info, msst_info, proc_info.dwt_lv);

    if(src.size() != msst_info.w * msst_info.h){
        cerr << "Error: Source msst size does not match expected dimensions." << endl;
        return;
    }

    //-- Inverse Quantization
    cquant.iquant_msstv(src);

    //-- convert msstSm to msst
    ConvMsstSm2Msst( msstv, src);

    //-- Perform inverse DWT on the msst data
    cdwt.dec_dwt53(msstv);

    //-- Perform inverse MSST on the DWT data
    cmsst.dec_msst(dst, msstv);

    //-- Perform inverse NLT on the reconstructed Bayer image if gamma is not 1.0
    nlt.inlt(dst);

}
