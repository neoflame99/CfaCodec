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
void fycc(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    if(bayer_info.cfa_pat&0x2){  // blue line first
        bl = sp;
        rl = (sp+bayer_info.w);
    }else{                       // red line first
        rl = sp;
        bl = (sp+bayer_info.w);
    }

    int32_t r, b, g1, g2, gg1, gg2;
    int32_t g1l,g1r, g2l, g2r;
    int32_t cr, cb, y, dg, cc ;
    if(bayer_info.cfa_pat&0x1){  // green ahead of R
        for(int k=0; k < bayer_info.w; k+=2){
            g1 = src[rl+k]; r = src[rl+k+1 ];
            b  = src[bl+k]; g2= src[bl+k+1 ];

            g1l = g1;
            g1r = k== bayer_info.w-2? g1 : src[rl+k+2];
            g2l = k==0 ? g2: src[bl+k-1];
            g2r = g2;

            gg1 = (g1l+g1r)>>1;
            gg2 = (g2l+g2r)>>1;
            cr = r - gg1;
            cb = b - gg2;
            cc = (cr+cb)>>1;
            y  = g1 + cc;
            dg = g2 - g1;
            dst[dp].Y    = y ;
            dst[dp].YdDg = dg;
            dst[dp].CrCo = cr;
            dst[dp].CbCg = cb;
            dp++;
        }
    }else{                       // green ahead of B
        for(int k=0; k < bayer_info.w; k+=2){
            r = src[rl+k  ]; g1 = src[rl+k+1];
            g2= src[bl+k  ];  b = src[bl+k+1];

            g1l = k==0 ? g1: src[rl+k-1];
            g1r = g1;
            g2l = g2;
            g2r = k== bayer_info.w-2? g2 : src[bl+k+2];

            gg1 = (g1l+g1r)>>1;
            gg2 = (g2l+g2r)>>1;
            cr = r - gg1;
            cb = b - gg2;
            cc = (cr+cb)>>1;
            y  = g1 + cc;
            dg = g2 - g1;
            dst[dp].Y    = y ;
            dst[dp].YdDg = dg;
            dst[dp].CrCo = cr;
            dst[dp].CbCg = cb;
            dp++;
        }
    }
}
void iycc(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
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

    int32_t r, b, g1, g2, g1l, g1r, g2l, g2r;
    int32_t y, dg, cr, cb, cc, gg, gg1, gg2;
    int32_t sp2;
    sp2 = sp; 
    if(bayer_info.cfa_pat&0x1){ // G1, R / B, G2
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y ;
            dg = src[sp].YdDg;
            cr = src[sp].CrCo;
            cb = src[sp].CbCg;
            sp++;

            cc = (cr+cb)>>1;
            g1 = y - cc; 
            g2 = dg + g1;

            dst[rl+k  ] = g1 ; 
            dst[bl+k+1] = g2;  
        }
        for(int k=0; k < bayer_info.w; k+=2){
            cr = src[sp2].CrCo;
            cb = src[sp2].CbCg;
            sp2++;

            g1 = dst[rl+k  ]; 
            g2 = dst[bl+k+1];

            g1l= g1;
            g1r= k== bayer_info.w-2? g1: dst[rl+k+2];
            g2l= k== 0? g2: dst[bl+k-1];
            g2r= g2;

            gg1 = (g1l+g1r)>>1;
            gg2 = (g2l+g2r)>>1;

            r  = cr + gg1;
            b  = cb + gg2; 

            dst[rl+k+1] = r ; 
            dst[bl+k  ] = b ; 
        }

    }else{   // R,G1 / G2,B 
        for(int k=0; k < bayer_info.w; k+=2){
            y  = src[sp].Y   ;
            dg = src[sp].YdDg;
            cr = src[sp].CrCo;
            cb = src[sp].CbCg;
            sp++;

            cc = (cr+cb)>>1;
            g1 = y - cc; 
            g2 = dg + g1;

            dst[rl+k+1]= g1;
            dst[bl+k  ]= g2;
        }
        for(int k=0; k < bayer_info.w; k+=2){
            cr = src[sp2].CrCo;
            cb = src[sp2].CbCg;
            sp2++;

            g1 = dst[rl+k+1]; 
            g2 = dst[bl+k  ];

            g1l= k== 0? g1: dst[rl+k-1];
            g1r= g1;
            g2l= g2;
            g2r= k== bayer_info.w-2? g2: dst[bl+k+2];

            gg1 = (g1l+g1r)>>1;
            gg2 = (g2l+g2r)>>1;

            r  = cr + gg1;
            b  = cb + gg2; 

            dst[rl+k  ]= r; 
            dst[bl+k+1]= b;
        }
    }
}
