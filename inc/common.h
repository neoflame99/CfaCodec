#include <cstdio>
#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <cassert>
#include <cmath>

using namespace std;

#ifndef _CFA_COMMON_H_
#define _CFA_COMMON_H_

//-- DWT level
#define DWTL 1   // upto 4
//-- Maximum Number of bands
//   1   |-------|------|
//   2   |---|---|------|
//   3   |-|-|---|------|
#define MNB (4+1)
//-- Maximum Bit Count
#define MBC (16+DWTL)
//-- Bit Count Width
#define BCW 5
// SGP: signs of group
#define SGP 4
// SGPCD: code length for signs of group, 1 or 2
#define SGPCD 1 

struct msstSm;

typedef int32_t cfapix;
struct cfa_fmt{
    cfapix *pCfa;
    uint32_t w;
    uint32_t bpp;
};
struct BayerInfo{
    int32_t w;
    int32_t h;
    int32_t cfa_pat;
    int32_t bpp;
};
struct BayerImgInfo{
    std::string filename;
    int32_t w;
    int32_t h;
    int32_t bpp;
    bool csi2_style = false; // false: compact bitstream, true: CSI-2 RAW packing
    bool excess_bytes_as_offset = true; // true: skip leading excess, false: drop trailing excess
};
struct SaveInfo{
    string bitstream_fname;
    string codec_report_fname;
    string dec_cfa_fname;
    bool dec_cfa_csi2_style;
    bool save_cfaimg;
    bool save_cfaproc;
    bool save_msstv_enc;
    bool save_msstv_dec;
};
struct ProcessInfo{
    float   g;     // for NLT
    uint32_t dwt_lv;
    uint32_t ngrp; // number of pixels in a group, 4 or 8
    uint32_t bcw;  // Bit Count Width, 5
    uint32_t mbc;  // Maximum Bit Count (16+DWTL)
    uint32_t sgpcd;// code length for signs of a group, 1 or 2
    bool     sel_ycc; // false: msst, true: ycc
    bool     sel_entp;// false: entropy_bc, true: entropy_gr
};

struct MsstInfo{
    int32_t w;
    int32_t h;
};
struct QuantInfo{
    int32_t Qp;
    int32_t Rp;
    int32_t Pb[MNB];
    int32_t Gb[MNB];
    int32_t Tb[MNB];
    QuantInfo(): Qp(0), Rp(0){
        for(int32_t k=0; k < MNB; ++k){
            Pb[k] = 0;
            Gb[k] = 0;
        }
    }
    void getTb(){
        for(int32_t k=0; k < MNB; ++k){
            int32_t r ; 
            if( Pb[k] < Rp){
                r = 1;
            }else{
                r = 0;
            }
            Tb[k] = Qp - Gb[k]-r;
            Tb[k] = (Tb[k] < 0) ? 0 : (Tb[k] > 15) ? 15: Tb[k];
        }
    }
};
struct msst{
    int32_t Y ; // Y
    int32_t YdDg; // Yd  
    int32_t CrCo; // Cr = R - G1 or R - ((G1 + G2)>>1)
    int32_t CbCg; // Cb = B - G2 or G - ((G1 + G2)>>1)
    msst():Y(0),YdDg(0),CrCo(0),CbCg(0){}

    msst(const msst& a): Y(a.Y), YdDg(a.YdDg), CrCo(a.CrCo), CbCg(a.CbCg){}
    msst(const msstSm& a){
        fromMsstSm(a);
    }

    inline void operator=(const msstSm& a){
        fromMsstSm(a);
    }
    inline void fromMsstSm(const msstSm& a);
    friend inline msst operator+(const int32_t lhs, const msst& rhs){
        msst ret;
        ret.Y  = lhs + rhs.Y;
        ret.YdDg = lhs + rhs.YdDg;
        ret.CrCo = lhs + rhs.CrCo;
        ret.CbCg = lhs + rhs.CbCg;
        return ret;
    }
    friend inline msst operator-(const int32_t lhs, const msst& rhs){
        msst ret;
        ret.Y  = lhs - rhs.Y;
        ret.YdDg = lhs - rhs.YdDg;
        ret.CrCo = lhs - rhs.CrCo;
        ret.CbCg = lhs - rhs.CbCg;
        return ret;
    }
    friend inline msst operator*(const int32_t lhs, const msst& rhs){
        msst ret;
        ret.Y  = lhs * rhs.Y;
        ret.YdDg = lhs * rhs.YdDg;
        ret.CrCo = lhs * rhs.CrCo;
        ret.CbCg = lhs * rhs.CbCg;
        return ret;
    }
    inline msst operator+(const msst& rhs) const{
        msst ret;
        ret.Y  = Y  + rhs.Y;
        ret.YdDg = YdDg + rhs.YdDg;
        ret.CrCo = CrCo + rhs.CrCo;
        ret.CbCg = CbCg + rhs.CbCg;
        return ret;
    }
    inline msst operator-(const msst& rhs) const{
        msst ret;
        ret.Y  = Y  - rhs.Y;
        ret.YdDg = YdDg - rhs.YdDg;
        ret.CrCo = CrCo - rhs.CrCo;
        ret.CbCg = CbCg - rhs.CbCg;
        return ret;
    }
    inline msst operator+(const int32_t& rhs) const{
        msst ret;
        ret.Y  = Y  + rhs;
        ret.YdDg = YdDg + rhs;
        ret.CrCo = CrCo + rhs;
        ret.CbCg = CbCg + rhs;
        return ret;
    }
    inline msst operator-(const int32_t& rhs) const{
        msst ret;
        ret.Y  = Y  - rhs;
        ret.YdDg = YdDg - rhs;
        ret.CrCo = CrCo - rhs;
        ret.CbCg = CbCg - rhs;
        return ret;
    }
    inline msst operator*(const int32_t& rhs) const{
        msst ret;
        ret.Y  = Y  * rhs;
        ret.YdDg = YdDg * rhs;
        ret.CrCo = CrCo * rhs;
        ret.CbCg = CbCg * rhs;
        return ret;
    }
    inline msst operator/(const int32_t& rhs) const{
        msst ret;
        ret.Y  = Y  / rhs;
        ret.YdDg = YdDg / rhs;
        ret.CrCo = CrCo / rhs;
        ret.CbCg = CbCg / rhs;
        return ret;
    }
    inline msst& operator+=(const msst& rhs){
        Y  += rhs.Y;
        YdDg += rhs.YdDg;
        CrCo += rhs.CrCo;
        CbCg += rhs.CbCg;
        return *this;
    }
    inline msst& operator-=(const msst& rhs){
        Y  -= rhs.Y;
        YdDg -= rhs.YdDg;
        CrCo -= rhs.CrCo;
        CbCg -= rhs.CbCg;
        return *this;
    }
    inline msst& operator*=(const int32_t& rhs){
        Y    *= rhs;
        YdDg *= rhs;
        CrCo *= rhs;
        CbCg *= rhs;
        return *this;
    }
    inline msst& operator/=(const int32_t& rhs){
        assert(rhs != 0);
        Y    /= rhs;
        YdDg /= rhs;
        CrCo /= rhs;
        CbCg /= rhs;
        return *this;
    }
    inline bool operator==(const msst& rhs) const{
        return (Y==rhs.Y && YdDg==rhs.YdDg && CrCo==rhs.CrCo && CbCg==rhs.CbCg);
    }
    inline bool operator!=(const msst& rhs) const{
        return !(*this==rhs);
    }
    inline msst& operator=(const msst& rhs){
        Y  = rhs.Y;
        YdDg = rhs.YdDg;
        CrCo = rhs.CrCo;
        CbCg = rhs.CbCg;
        return *this;
    }
    inline msst& operator>>(const int32_t& rhs){
        Y  >>= rhs;
        YdDg >>= rhs;
        CrCo >>= rhs;
        CbCg >>= rhs;
        return *this;
    }
    inline msst& operator<<(const int32_t& rhs){
        Y  <<= rhs;
        YdDg <<= rhs;
        CrCo <<= rhs;
        CbCg <<= rhs;
        return *this;
    }
};

struct msstSm{
    uint32_t mY;
    uint32_t mYdDg;
    uint32_t mCrCo;
    uint32_t mCbCg;
    bool    sY;
    bool    sYdDg;
    bool    sCrCo;
    bool    sCbCg;
    msstSm():mY(0), mYdDg(0), mCbCg(0), mCrCo(0), sY(0), sYdDg(0), sCbCg(0), sCrCo(0)
    {}
    msstSm(const msst& a){
        fromMsst(a);
    }
    inline void fromMsst(const msst& a){
        mY = abs(a.Y );
        mYdDg= abs(a.YdDg);
        mCrCo= abs(a.CrCo);
        mCbCg= abs(a.CbCg);
        sY = a.Y  >> (sizeof(a.Y )*8-1);
        sYdDg= a.YdDg >> (sizeof(a.YdDg)*8-1);
        sCrCo= a.CrCo >> (sizeof(a.CrCo)*8-1);
        sCbCg= a.CbCg >> (sizeof(a.CbCg)*8-1);
    }
    inline void operator=(const msst& a){
        fromMsst(a);
    }
};

inline void msst::fromMsstSm(const msstSm& a){
    Y  = a.sY ? -a.mY : a.mY;
    YdDg = a.sYdDg? -a.mYdDg: a.mYdDg;
    CrCo = a.sCrCo? -a.mCrCo: a.mCrCo;
    CbCg = a.sCbCg? -a.mCbCg: a.mCbCg;
}

inline void ConvMsst2MsstSm(vector<msstSm>& dst, const vector<msst>& src){
    size_t dlen = dst.size();
    size_t slen = src.size();
    if( dlen != slen){
        dst.resize(slen);
    }
    for(size_t k=0; k < slen; ++k){
        dst[k] = src[k];
    }
}
inline void ConvMsstSm2Msst(vector<msst>& dst, const vector<msstSm>& src){
    size_t dlen = dst.size();
    size_t slen = src.size();
    if( dlen != slen){
        dst.resize(slen);
    }
    for(size_t k=0; k < slen; ++k){
        dst[k] = src[k];
    }
}


#endif
