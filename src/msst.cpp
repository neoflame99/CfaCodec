#include "msst.h"

void fmsst(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    if(bayer_info.cfa_pat&0x2){
        bl = sp;
        rl = (sp+bayer_info.w);
    }else{
        rl = sp;
        bl = (sp+bayer_info.w);
    }

    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    if(bayer_info.cfa_pat&0x1){
        for(int k=0; k < bayer_info.w; k+=2){
            g1 = src[rl+k]; r = src[rl+k+1 ];
            b  = src[bl+k]; g2= src[bl+k+1 ];
            co = r - b;
            dg = g2 - g1;
            u  = b + (co >> 1);
            v  = g1+ (dg >> 1);
            cg = v - u;
            y  = u + (cg >> 1);
            dst[dp].Y  = y;
            dst[dp].YdDg = dg;
            dst[dp].CrCo = co;
            dst[dp].CbCg = cg;
            dp++;
        }
    }else{
        for(int k=0; k < bayer_info.w; k+=2){
            r = src[rl+k  ]; g1 = src[rl+k+1];
            g2= src[bl+k  ];  b = src[bl+k+1];
            co = r - b;
            dg = g2 - g1;
            u  = b + (co >> 1);
            v  = g1+ (dg >> 1);
            cg = v - u;
            y  = u + (cg >> 1);
            dst[dp].Y  = y;
            dst[dp].YdDg = dg;
            dst[dp].CrCo = co;
            dst[dp].CbCg = cg;
            dp++;
        }
    }
}

void imsst(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
         const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    if(bayer_info.cfa_pat&0x2){
        bl = dp;
        rl = (dp+bayer_info.w);
    }else{
        rl = dp;
        bl = (dp+bayer_info.w);
    }

    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    if(bayer_info.cfa_pat&0x1){
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y ;
            dg = src[sp].YdDg;
            co = src[sp].CrCo;
            cg = src[sp].CbCg;
            sp++;

            u  = y - (cg >> 1); 
            v  = u + cg;
            g1 = v - (dg >> 1);
            b  = u - (co >> 1);
            g2 = g1 + dg;
            r  = b + co;

            dst[rl+k] = g1 ; dst[rl+k+1 ] = r ; 
            dst[bl+k] = b  ; dst[bl+k+1 ] = g2;  
        }
    }else{
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y ;
            dg = src[sp].YdDg;
            co = src[sp].CrCo;
            cg = src[sp].CbCg;
            sp++;

            u  = y - (cg >> 1); 
            v  = u + cg;
            g1 = v - (dg >> 1);
            b  = u - (co >> 1);
            g2 = g1 + dg;
            r  = b + co;

            dst[rl+k  ]= r ;  dst[rl+k+1]= g1;
            dst[bl+k  ]= g2;  dst[bl+k+1]=  b;
        }
    }
}
void fycbcr(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    if(bayer_info.cfa_pat&0x2){
        bl = sp;
        rl = (sp+bayer_info.w);
    }else{
        rl = sp;
        bl = (sp+bayer_info.w);
    }

    bool sel= false;
    int32_t r, b, g1, g2, gg1, gg2, gg;
    int32_t cr, cb, y1, y2, y, yd, cc;
    if(bayer_info.cfa_pat&0x1){
        for(int k=0; k < bayer_info.w; k+=2){
            g1 = src[rl+k]; r = src[rl+k+1 ];
            b  = src[bl+k]; g2= src[bl+k+1 ];
            gg = (g1+g2) >> 1;
            gg1 = (sel)? g1: gg;
            gg2 = (sel)? g2: gg;
            cr = r - gg1;
            cb = b - gg2;
            cc = (cr+cb)>>1;
            y1 = g1 + cc;
            y2 = g2 + cc;
            yd = y1 - y2;
            y  = y1;
            dst[dp].Y  = y ;
            dst[dp].YdDg = yd;
            dst[dp].CrCo = cr;
            dst[dp].CbCg = cb;
            dp++;
        }
    }else{
        for(int k=0; k < bayer_info.w; k+=2){
            r = src[rl+k  ]; g1 = src[rl+k+1];
            g2= src[bl+k  ];  b = src[bl+k+1];
            gg = (g1+g2) >> 1;
            gg1 = (sel)? g1: gg;
            gg2 = (sel)? g2: gg;
            cr = r - gg1;
            cb = b - gg2;
            cc = (cr+cb)>>1;
            y1 = g1 + cc;
            y2 = g2 + cc;
            yd = y1 - y2;
            y  = y1;
            dst[dp].Y  = y ;
            dst[dp].YdDg = yd;
            dst[dp].CrCo = cr;
            dst[dp].CbCg = cb;
            dp++;
        }
    }
}
void iycbcr(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
         const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    if(bayer_info.cfa_pat&0x2){
        bl = dp;
        rl = (dp+bayer_info.w);
    }else{
        rl = dp;
        bl = (dp+bayer_info.w);
    }

    bool sel = false;
    int32_t r, b, g1, g2;
    int32_t y, yd, y1, y2, cr, cb, cc, gg, gg1, gg2;
    if(bayer_info.cfa_pat&0x1){
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y ;
            yd = src[sp].YdDg;
            cr = src[sp].CrCo;
            cb = src[sp].CbCg;
            sp++;

            y1 = y;
            y2 = y-yd;
            cc = (cr+cb)>>1;
            g1 = y1 - cc; 
            g2 = y2 - cc;
            gg = (g1+g2) >> 1;
            gg1 = sel ? g1 : gg;
            gg2 = sel ? g2 : gg;
            r  = cr + gg1;
            b  = cb + gg2; 

            dst[rl+k] = g1 ; dst[rl+k+1 ] = r ; 
            dst[bl+k] = b  ; dst[bl+k+1 ] = g2;  
        }
    }else{
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y ;
            yd = src[sp].YdDg;
            cr = src[sp].CrCo;
            cb = src[sp].CbCg;
            sp++;

            y1 = y;
            y2 = y-yd;
            cc = (cr+cb)>>1;
            g1 = y1 - cc; 
            g2 = y2 - cc;
            gg = (g1+g2) >> 1;
            gg1 = sel ? g1 : gg;
            gg2 = sel ? g2 : gg;
            r  = cr + gg1;
            b  = cb + gg2; 

            dst[rl+k  ]= r ;  dst[rl+k+1]= g1;
            dst[bl+k  ]= g2;  dst[bl+k+1]=  b;
        }
    }
}