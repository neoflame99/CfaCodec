#include "msst.h"

static inline void fmsst(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info, const ProcessInfo& proc_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    const BayerPhase phase = bayer_phase(bayer_info.cfa_pat);
    if(phase.red_on_odd_row){
        bl = sp;
        rl = (sp+bayer_info.w);
    }else{
        rl = sp;
        bl = (sp+bayer_info.w);
    }

    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    if(phase.red_on_odd_col){
        for(int k=0; k < bayer_info.w; k+=2){
            g1 = src[rl+k]; r = src[rl+k+1 ];
            b  = src[bl+k]; g2= src[bl+k+1 ];
            co = r - b;
            dg = g2 - g1;
            u  = b + (co >> 1);
            v  = g1+ (dg >> 1);
            cg = v - u;
            y  = u + (cg >> 1);
            if(proc_info.sel_ycc){
                y = g1;
                dg= g2;
                co= r ;
                cg= b ;
            }
            dst[dp].Y    = y;
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
            if(proc_info.sel_ycc){
                y = g1;
                dg= g2;
                co= r ;
                cg= b ;
            }
            dst[dp].Y    = y;
            dst[dp].YdDg = dg;
            dst[dp].CrCo = co;
            dst[dp].CbCg = cg;
            dp++;
        }
    }
}

static inline void imsst(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
         const BayerInfo& bayer_info, const ProcessInfo& proc_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    const BayerPhase phase = bayer_phase(bayer_info.cfa_pat);
    if(phase.red_on_odd_row){
        bl = dp;
        rl = (dp+bayer_info.w);
    }else{
        rl = dp;
        bl = (dp+bayer_info.w);
    }

    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    if(phase.red_on_odd_col){
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
            if(proc_info.sel_ycc){
                g1 = y;
                g2 = dg;
                r  = co;
                b  = cg;
            }

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
            if(proc_info.sel_ycc){
                g1 = y;
                g2 = dg;
                r  = co;
                b  = cg;
            }

            dst[rl+k  ]= r ;  dst[rl+k+1]= g1;
            dst[bl+k  ]= g2;  dst[bl+k+1]=  b;
        }
    }
}
static inline void fmsst_ir(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info, const ProcessInfo& proc_info, int32_t ro){
    //cfa_pat:
    //4: rgbg/gi, 5: gbrg/ig,  6: gi/rgbg,  7: ig/gbgr
    //8: bgrg/gi, 9: grgb/ig, 10: gi/bgrg, 11: ig/grgb
    //
    //   0  1  2  3     0  1  2  3     0  1  2  3     0  1  2  3 
    // 0 R  G  B  G   0 G  B  G  R   0 G  I  G  I   0 I  G  I  G
    // 1 G  I  G  I   1 I  G  I  G   1 R  G  B  G   1 G  B  G  R
    // 2 B  G  R  G   2 G  R  G  B   2 G  I  G  I   3 I  G  I  G
    // 3 G  I  G  I   3 I  G  I  G   3 B  G  R  G   2 G  R  G  B
    //    00b:00b        00b:11b        01b:00b        01b:11b  
    //
    //   0  1  2  3     0  1  2  3     0  1  2  3     0  1  2  3 
    // 0 B  G  R  G   0 G  R  G  B   0 G  I  G  I   0 I  G  I  G  
    // 1 G  I  G  I   1 I  G  I  G   1 B  G  R  G   1 G  R  G  B  
    // 2 R  G  B  G   2 G  B  G  R   2 G  I  G  I   3 I  G  I  G  
    // 3 G  I  G  I   3 I  G  I  G   3 R  G  B  G   2 G  B  G  R  
    //    10b:00b        10b:11b        11b:00b        11b:11b  
    //
    //  x, y are the coordinates of the top-left pixel of a 4x4 block.
    //  ypat_code, xpat_code are the 2-bit codes for the y and x coordinates.
    //  ycode = ypat_code ^ y, xcode = xpat_code ^ x
    //  color position code =  ycode ^ xcode.
    
    int32_t pat = 0, rpat, rcode, cpat;
    pat =static_cast<int32_t>(bayer_info.cfa_pat); 
    if( pat < 4 || pat > 11){
        cout << "fmsst_ir: invalid cfa_pat = " << static_cast<int32_t>(bayer_info.cfa_pat) << endl;
        return;
    }
    pat -= 4;
    switch(pat){
        case 0: rpat = 0x0; cpat = 0x0; break;
        case 1: rpat = 0x0; cpat = 0x3; break; 
        case 2: rpat = 0x1; cpat = 0x0; break; 
        case 3: rpat = 0x1; cpat = 0x3; break; 
        case 4: rpat = 0x2; cpat = 0x0; break;
        case 5: rpat = 0x2; cpat = 0x3; break;
        case 6: rpat = 0x3; cpat = 0x0; break;
        case 7: rpat = 0x3; cpat = 0x3; break;
        default: break;
    }
    rcode = (ro & 0x3) ^ rpat;
    int32_t ri, bi, g1i, g2i;
    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    for(int k=0; k < bayer_info.w; k+=4){
        // k's lsb 2 bits are always 00b because k starts 0 and increments by 4.
        ri  = k | (rcode ^ (0x0 ^ cpat));
        g1i = k | (rcode ^ (0x1 ^ cpat));
        bi  = k | (rcode ^ (0x2 ^ cpat));
        g2i = k | (rcode ^ (0x3 ^ cpat));

        r  = src[ri ];
        g1 = src[g1i]; 
        b  = src[bi ];
        g2 = src[g2i];
        co = r - b;
        dg = g2 - g1;
        u  = b + (co >> 1);
        v  = g1+ (dg >> 1);
        cg = v - u;
        y  = u + (cg >> 1);
        dst[dp].Y    = y;
        dst[dp].YdDg = dg;
        dst[dp].CrCo = co;
        dst[dp].CbCg = cg;
        dp++;
    }
}

static inline void imsst_ir(vector<cfapix>& dst, int32_t& dp, const vector<msst>& src, int32_t& sp,
         const BayerInfo& bayer_info, const ProcessInfo& proc_info, int32_t ro){
    //cfa_pat:
    //4: rgbg/gi, 5: gbrg/ig,  6: gi/rgbg,  7: ig/gbgr
    //8: bgrg/gi, 9: grgb/ig, 10: gi/bgrg, 11: ig/grgb
    //
    //   0  1  2  3     0  1  2  3     0  1  2  3     0  1  2  3 
    // 0 R  G  B  G   0 G  B  G  R   0 G  I  G  I   0 I  G  I  G
    // 1 G  I  G  I   1 I  G  I  G   1 R  G  B  G   1 G  B  G  R
    // 2 B  G  R  G   2 G  R  G  B   2 G  I  G  I   3 I  G  I  G
    // 3 G  I  G  I   3 I  G  I  G   3 B  G  R  G   2 G  R  G  B
    //    00b:00b        00b:11b        01b:00b        01b:11b  
    //
    //   0  1  2  3     0  1  2  3     0  1  2  3     0  1  2  3 
    // 0 B  G  R  G   0 G  R  G  B   0 G  I  G  I   0 I  G  I  G  
    // 1 G  I  G  I   1 I  G  I  G   1 B  G  R  G   1 G  R  G  B  
    // 2 R  G  B  G   2 G  B  G  R   2 G  I  G  I   3 I  G  I  G  
    // 3 G  I  G  I   3 I  G  I  G   3 R  G  B  G   2 G  B  G  R  
    //    10b:00b        10b:11b        11b:00b        11b:11b  
    //
    //  x, y are the coordinates of the top-left pixel of a 4x4 block.
    //  ypat_code, xpat_code are the 2-bit codes for the y and x coordinates.
    //  ycode = ypat_code ^ y, xcode = xpat_code ^ x
    //  color position code =  ycode ^ xcode.
    
    int32_t pat = 0, rpat, rcode, cpat;
    pat =static_cast<int32_t>(bayer_info.cfa_pat); 
    if( pat < 4 || pat > 11){
        cout << "fmsst_ir: invalid cfa_pat = " << static_cast<int32_t>(bayer_info.cfa_pat) << endl;
        return;
    }
    pat -= 4;
    switch(pat){
        case 0: rpat = 0x0; cpat = 0x0; break;
        case 1: rpat = 0x0; cpat = 0x3; break; 
        case 2: rpat = 0x1; cpat = 0x0; break; 
        case 3: rpat = 0x1; cpat = 0x3; break; 
        case 4: rpat = 0x2; cpat = 0x0; break;
        case 5: rpat = 0x2; cpat = 0x3; break;
        case 6: rpat = 0x3; cpat = 0x0; break;
        case 7: rpat = 0x3; cpat = 0x3; break;
        default: break;
    }
    rcode = (ro & 0x3) ^ rpat;
    int32_t ri, bi, g1i, g2i;
    int32_t co, dg, u,v, cg, y, r, b, g1, g2;
    for(int k=0; k < bayer_info.w; k+=4){
        // k's lsb 2 bits are always 00b because k starts 0 and increments by 4.
        ri  = k | (rcode ^ (0x0 ^ cpat));
        g1i = k | (rcode ^ (0x1 ^ cpat));
        bi  = k | (rcode ^ (0x2 ^ cpat));
        g2i = k | (rcode ^ (0x3 ^ cpat));

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
        dst[ri ] = r ; 
        dst[g1i] = g1; 
        dst[bi ] = b ;
        dst[g2i] = g2;  
    }
}

static void sat(vector<cfapix>& dst, const uint32_t mxv){
    for(size_t k=0; k < dst.size(); ++k){
        if(dst[k] < 0){
            dst[k] = 0;
        }else if(dst[k] > mxv){
            dst[k] = mxv;
        }
    }
}
void MSST::enc_msst(vector<msst>& dst, const vector<cfapix>& src){
    int32_t dp = 0, sp = 0;
    bool sel_ir = (static_cast<int32_t>(bayer_info.cfa_pat) >= 4);
    if(sel_ir){
        for(int32_t h = 0; h < bayer_info.h; h++){
            fmsst_ir(dst, dp, src, sp, bayer_info, proc_info, h);
            sp += bayer_info.w;
        }
    }else{
        for(int32_t h = 0; h < bayer_info.h; h+=2){
            fmsst(dst, dp, src, sp, bayer_info, proc_info);
            sp += bayer_info.w*2;
        }
    }
    cout << "bayer_fdwt, fmsst: dp = " << dp << ", sp = " << sp << endl;
}
void MSST::dec_msst(vector<cfapix>& dst, const vector<msst>& src){
    int32_t dp = 0, sp = 0;
    bool sel_ir = (static_cast<int32_t>(bayer_info.cfa_pat) >= 4);
    if(sel_ir){
        for(int32_t h = 0; h < bayer_info.h; h++){
            imsst_ir(dst, dp, src, sp, bayer_info, proc_info, h);
            dp += bayer_info.w;
        }
    }else{
        for(int32_t h = 0; h < bayer_info.h; h+=2){
            imsst(dst, dp, src, sp, bayer_info, proc_info);
            dp += bayer_info.w*2;
        }
    }
    uint32_t mxv = (1U << bayer_info.bpp) - 1;
    sat(dst, mxv);
    cout << "bayer_idwt, imsst: dp = " << dp << ", sp = " << sp << endl;
}

void fycc(vector<msst>& dst, int32_t& dp, const vector<cfapix>& src, int32_t& sp,
    const BayerInfo& bayer_info){
    //cfa_pat:
    //0: rg1/g2b, 1: g1r/bg2, 2: g2b/rg1, 3: bg2/g1r
    int32_t rl, bl;
    const BayerPhase phase = bayer_phase(bayer_info.cfa_pat);
    if(phase.red_on_odd_row){  // blue line first
        bl = sp;
        rl = (sp+bayer_info.w);
    }else{                       // red line first
        rl = sp;
        bl = (sp+bayer_info.w);
    }

    int32_t r, b, g1, g2, gg1, gg2;
    int32_t g1l,g1r, g2l, g2r;
    int32_t cr, cb, y, dg, cc ;
    if(phase.red_on_odd_col){  // green ahead of R
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
    const BayerPhase phase = bayer_phase(bayer_info.cfa_pat);
    if(phase.red_on_odd_row){
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
    if(phase.red_on_odd_col){ // G1, R / B, G2
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
