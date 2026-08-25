#include "quantize.h"

static inline void quant(vector<msstSm>& msstv, const int32_t S, const int32_t N, const int32_t B, const QuantInfo& quant_info){
    // B : band, HH 0, HL 1, LH 2, HH 3 ...
    // S : band start index
    // N : band length 

    int32_t Tb = quant_info.Tb[B];
    for(int32_t k =S; k < N+S; ++k){
        //if( k == 221761){
        //    printf("k = %d, Tb = %d\n", k, Tb);
        //    printf("mY = %d, mYdDg = %d, mCrCo = %d, mCbCg = %d\n", msstv[k].mY, msstv[k].mYdDg, msstv[k].mCrCo, msstv[k].mCbCg);
        //}
        msstv[k].mY  >>= Tb;
        msstv[k].mYdDg >>= Tb; 
        msstv[k].mCrCo >>= Tb;
        msstv[k].mCbCg >>= Tb;
        //if( k == 221761){
        //    printf("mY = %d, mYdDg = %d, mCrCo = %d, mCbCg = %d\n", msstv[k].mY, msstv[k].mYdDg, msstv[k].mCrCo, msstv[k].mCbCg);
        //    int a=10;
        //}
    }
    // truncation bits and quantization step to be determined 
}
static inline void iquant(vector<msstSm>& msstv, const int32_t S, const int32_t N, const int32_t B, const QuantInfo& quant_info){
    // B : band, HH 0, HL 1, LH 2, HH 3 ...
    // S : band start index
    // N : band length 
    int32_t Tb = quant_info.Tb[B];
    for(int32_t k =S; k < N+S; ++k){
        //if( k == 221761){
        //    printf("k = %d, Tb = %d\n", k, Tb);
        //    printf("mY = %d, mYdDg = %d, mCrCo = %d, mCbCg = %d\n", msstv[k].mY, msstv[k].mYdDg, msstv[k].mCrCo, msstv[k].mCbCg);
        //}
        msstv[k].mY  <<= Tb;
        msstv[k].mYdDg <<= Tb; 
        msstv[k].mCrCo <<= Tb;
        msstv[k].mCbCg <<= Tb;
        //if( k == 221761){
        //    printf("mY = %d, mYdDg = %d, mCrCo = %d, mCbCg = %d\n", msstv[k].mY, msstv[k].mYdDg, msstv[k].mCrCo, msstv[k].mCbCg);
        //    int a=10;
        //}
    }
    // truncation bits and quantization step to be determined 
}

void QUANT::fquant_msstv(vector<msstSm>& msstv){
    int32_t N, S, B, O;
    O = 0;
    for(int32_t m=0; m < msst_info.h; ++m, O+= msst_info.w){
        B = 0;
        for(int32_t n=1; n <= dwt_lv; ++n){
            N = msst_info.w >> n;
            S = O+N; 
            quant(msstv, S, N, B, quant_info );
            B++;
        }
        N = msst_info.w >> dwt_lv;
        S = O;
        quant(msstv, S, N, B, quant_info);
        // dwt_lv = 1
        // N = w >> 1
        // B=0: S = N = w/2
        // dwt_lv = 2
        // N = w >> 2
        // B=1: S = N = w/4
        // dwt_lv = 2
        // N = w >> 2
        // B=2: S = 0
    }
}
void QUANT::iquant_msstv(vector<msstSm>& msstv){
    int32_t N, S, B, O;
    O = 0;
    for(int32_t m=0; m < msst_info.h; ++m, O+= msst_info.w){
        B = 0;
        for(int32_t n=1; n <= dwt_lv; ++n){
            N = msst_info.w >> n;
            S = O+ N; 
            iquant(msstv, S, N, B, quant_info );
            B++;
        }
        N = msst_info.w >> dwt_lv;
        S = O;
        iquant(msstv, S, N, B, quant_info);
    }
}