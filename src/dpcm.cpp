#include "dpcm.h"
#include <cassert>

void dpcm_mixc(cfa_fmt& ocfa, const cfa_fmt& icfa){
    cfapix* src = icfa.pCfa;
    cfapix* dst = ocfa.pCfa;
    assert(src!= nullptr && dst != nullptr );
    uint32_t w = icfa.w;
    ocfa.w   = w;
    ocfa.bpp = icfa.bpp;

    dst[0] = src[0];
    dst[1] = src[1];
    for(uint32_t k=2; k < w; ++k){
        dst[k] = src[k] - src[k-2];
    }
}

void dpcm_monoc(cfa_fmt& ocfa, const cfa_fmt& icfa){
    cfapix* src = icfa.pCfa;
    cfapix* dst = ocfa.pCfa;
    assert(src!= nullptr && dst != nullptr );
    uint32_t w = icfa.w;
    ocfa.w   = w;
    ocfa.bpp = icfa.bpp;

    dst[0] = src[0];
    for(uint32_t k=1; k < w; ++k){
        dst[k] = src[k] - src[k-1];
    }
}

void idpcm_mixc(cfa_fmt& ocfa, const cfa_fmt& icfa){
    cfapix* src = icfa.pCfa;
    cfapix* dst = ocfa.pCfa;
    assert(src!= nullptr && dst != nullptr );
    uint32_t w = icfa.w;
    ocfa.w   = w;
    ocfa.bpp = icfa.bpp;

    dst[0] = src[0];
    dst[1] = src[1];
    for(uint32_t k=2; k < w; ++k){
        dst[k] = dst[k-2] + src[k];
    }
}

void idpcm_monoc(cfa_fmt& ocfa, const cfa_fmt& icfa){
    cfapix* src = icfa.pCfa;
    cfapix* dst = ocfa.pCfa;
    assert(src!= nullptr && dst != nullptr );
    uint32_t w = icfa.w;
    ocfa.w   = w;
    ocfa.bpp = icfa.bpp;

    dst[0] = src[0];
    for(uint32_t k=1; k < w; ++k){
        dst[k] = dst[k-1] + src[k];
    }
}