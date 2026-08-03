#include "entropy.h"

struct GR{
    uint8_t q;
    uint8_t r;
    GR():q(0), r(0){}
};

static inline int8_t signed_qbit(uint8_t a){
    int8_t b = a & 0xf;
    int8_t s = (a & 0x8)>>3;
    for(int k=4; k < 8; ++k){
         b |= s << k;
    }
    return b;
}
static inline uint8_t fzigzag(int8_t a){
    uint8_t b;
    b = (a < 0)? (-a<<1)-1 : a << 1;
    return b;
}
static inline int8_t izigzag(uint8_t a){
    int8_t b;
    b = ( a &0x1)? -((a+1)>>1) : a >> 1;
    return b;
}
static inline void enc_gs_golombrice(GR& gr,const int8_t s,const uint8_t K){
    uint8_t zs = fzigzag(s);
    gr.q = zs >> K;
    gr.r = zs &((1U << K)-1);
}
static inline void dec_gs_golombrice(int8_t& s, const GR& gr, const uint8_t K){
    uint8_t q, d;
    q = gr.q << K;
    d = q | gr.r;
    s = izigzag(d);
}
static inline void extract_gr(GR& gr, const vector<uint8_t>& bstm, uint32_t& pos, uint32_t& bp, uint8_t K){
    uint8_t d = bstm[pos];
    d >>= bp;
    uint8_t q = 0, r=0;
    uint8_t b;
    // getting q
    for(;;){
        b = (d & 0x1);
        d >>= 1;
        q+= b;
        bp++;
        if(bp >= 8){
            ++pos;
            d = bstm[pos];
            bp= 0;
        }
        if(!b) break;
    }
    // getting r
    for(int k =0; k < K; ++k){
        b = (d & 0x1);
        d >>= 1;
        r |= b << k;
        bp++;
        if(bp >= 8){
            ++pos;
            d = bstm[pos];
            bp= 0;
        }
    }
}

static inline void fillbits(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t td, uint32_t blen){
    uint8_t b = rRow[pos];
    for(int32_t m=0; m < blen; ++m){
        b |= (td&0x1) << bp;
        td >>= 1;
        bp++;
        if(bp >=8){
            rRow[pos] = b;
            b = 0;
            bp=0;
            pos++;
        }
    }
    if(bp!=0){
        rRow[pos] = b;
    }
}
static inline void putgsbits(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t td, uint32_t blen){
    uint32_t mx = (1U << SGPCD)-1;
    if(blen >= SGPCD){
        // blen: 3, 4
        fillbits(rRow, pos, bp, td, SGP);
    }else{
        // blen: 0, 1, 2
        fillbits(rRow, pos, bp, td, blen);
    }
}
void enc_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info, 
    const ProcessInfo& proc_info, const QuantInfo& quant_info )
{
    //--   +----+----+----+----+----+----+----+----+----+
    //--   | 0  | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 0  |
    //--   +----+----+----+----+----+----+----+----+----+
    //--    <- 4 or 5-bit bcnt    -> <- signs of group->
    //--     HH band: 4-bit bcnt, others: 5-bit bcnt
    uint32_t my, mdg, mcg, mco;
    uint32_t ybcnt, dgbcnt, cgbcnt, cobcnt;
    uint32_t ybcnt_cl, dgbcnt_cl, cgbcnt_cl, cobcnt_cl;
    uint32_t yscnt, dgscnt, cgscnt, coscnt;
    uint32_t sy, sdg, scg, sco;
    uint32_t sy2, sdg2, scg2, sco2;
    uint32_t sycd, sdgcd, scgcd, scocd;
    int32_t k=0;
    uint32_t pos, bp;
    uint32_t tsz=0;
    uint32_t mx = (1U<<SGPCD)-1;
    int32_t  hwidth = msst_info.w/2;
    int32_t  tw;
    uint32_t bcw;
    #ifdef _DUMP_
    FILE* fp  = fopen("tmp/enc_ent.txt", "w");
    FILE* fp2 = fopen("tmp/enc_in.txt", "w");
    FILE* fp3 = fopen("tmp/qbit.txt", "w");
    fprintf(fp3,"## [sy, sdg, scg, sco]\n");
    #endif

    int32_t step = proc_info.ngrp;
    for(int32_t r=0; r < msst_info.h; ++r){
        vector<uint8_t>& rRow = bitv[r];
        pos = 4;
        bp  = 0;
        //for(int32_t c=0; c < msst_info.w; c+=step, k+=step){
        for(int32_t l=proc_info.dwt_lv+1; l > 0; --l){
            int32_t sl = l > proc_info.dwt_lv ? proc_info.dwt_lv : l; 
            tw = msst_info.w >> sl;
            for(int32_t c=0; c < tw; c+=step, k+=step){
                ybcnt=0 ; yscnt = 0;
                dgbcnt=0; dgscnt= 0;
                cgbcnt=0; cgscnt= 0;
                cobcnt=0; coscnt= 0;

                my = msstv[k].mY   ; 
                mdg= msstv[k].mYdDg; 
                mco= msstv[k].mCrCo; 
                mcg= msstv[k].mCbCg; 
                sy = int(msstv[k].sY)   ; 
                sdg= int(msstv[k].sYdDg); 
                sco= int(msstv[k].sCrCo); 
                scg= int(msstv[k].sCbCg); 
                for(int32_t s=k+1, ss=1; ss < step; ++ss, ++s){
                   my |= msstv[s].mY   ; //msstv[k].mY | msstv[k+1].mY | msstv[k+2].mY | msstv[k+3].mY;
                   mdg|= msstv[s].mYdDg; // msstv[k].mYdDg| msstv[k+1].mYdDg| msstv[k+2].mYdDg| msstv[k+3].mYdDg;
                   mco|= msstv[s].mCrCo; // msstv[k].mCrCo| msstv[k+1].mCrCo| msstv[k+2].mCrCo| msstv[k+3].mCrCo;
                   mcg|= msstv[s].mCbCg; // msstv[k].mCbCg| msstv[k+1].mCbCg| msstv[k+2].mCbCg| msstv[k+3].mCbCg;
                   sy |= int(msstv[s].sY)   << s ; //int(msstv[k].sY) | (int(msstv[k+1].sY) << 1) | (int(msstv[k+2].sY)<<2) | (int(msstv[k+3].sY)<<3);
                   sdg|= int(msstv[s].sYdDg)<< s ; //int(msstv[k].sYdDg) | (int(msstv[k+1].sYdDg) << 1) | (int(msstv[k+2].sYdDg)<<2) | (int(msstv[k+3].sYdDg)<<3);
                   sco|= int(msstv[s].sCrCo)<< s ; //int(msstv[k].sCrCo) | (int(msstv[k+1].sCrCo) << 1) | (int(msstv[k+2].sCrCo)<<2) | (int(msstv[k+3].sCrCo)<<3);
                   scg|= int(msstv[s].sCbCg)<< s ; //int(msstv[k].sCbCg) | (int(msstv[k+1].sCbCg) << 1) | (int(msstv[k+2].sCbCg)<<2) | (int(msstv[k+3].sCbCg)<<3);
                }
                sy2 = sy; sdg2 = sdg; scg2 = scg; sco2 = sco;
                #ifdef _DUMP_
                fprintf(fp2,"[k%4d, mY: %6d, %6d, %6d, %6d, sY: %1d, %1d, %1d, %1d]\n",k, msstv[k].mY, msstv[k+1].mY, msstv[k+2].mY, msstv[k+3].mY, msstv[k].sY, msstv[k+1].sY, msstv[k+2].sY, msstv[k+3].sY);
                ssy = signed_qbit(sy); ssdg = signed_qbit(sdg); sscg = signed_qbit(scg); ssco = signed_qbit(sco);
                int8_t ssy, ssdg, sscg, ssco;
                fprintf(fp3, "[k%6d, %2d, %2d, %2d %2d],", k, ssy, ssdg, sscg, ssco);
                #endif

                //if(k==283728){
                //    int a=10;
                //    int myy ; 
                //    myy = msstv[k].mCo;
                //    myy = msstv[k+1].mCo;
                //    myy = msstv[k+2].mCo;
                //    myy = msstv[k+3].mCo;
                //}
                for(int32_t m=0; m < MBC; ++m){
                    if( my == 0) break;
                    my >>= 1;
                    ybcnt++;
                }
                for(int32_t m=0; m < MBC; ++m){
                    if( mdg == 0) break;
                    mdg >>= 1;
                    dgbcnt++;
                }
                for(int32_t m=0; m < MBC; ++m){
                    if( mcg == 0) break;
                    mcg >>= 1;
                    cgbcnt++;
                }
                for(int32_t m=0; m < MBC; ++m){
                    if( mco == 0) break;
                    mco >>= 1;
                    cobcnt++;
                }
                for(int32_t m=0; m < SGP; ++m){
                    if( sy2 == 0) break;
                    sy2 >>= 1;
                    yscnt++;
                }
                for(int32_t m=0; m < SGP; ++m){
                    if( sdg2 == 0) break;
                    sdg2 >>= 1;
                    dgscnt++;
                }
                for(int32_t m=0; m < SGP; ++m){
                    if( scg2 == 0) break;
                    scg2 >>= 1;
                    cgscnt++;
                }
                for(int32_t m=0; m < SGP; ++m){
                    if( sco2 == 0) break;
                    sco2 >>= 1;
                    coscnt++;
                }

                ybcnt = ybcnt >= 1 ? ybcnt-1 : 0;
                dgbcnt= dgbcnt>= 1 ? dgbcnt-1: 0;
                cgbcnt= cgbcnt>= 1 ? cgbcnt-1: 0;
                cobcnt= cobcnt>= 1 ? cobcnt-1: 0;
                ybcnt_cl = ybcnt +1;
                dgbcnt_cl= dgbcnt+1;
                cgbcnt_cl= cgbcnt+1;
                cobcnt_cl= cobcnt+1;

                sycd = yscnt >=mx ? mx : yscnt;
                sdgcd= dgscnt>=mx ? mx : dgscnt;
                scgcd= cgscnt>=mx ? mx : cgscnt;
                scocd= coscnt>=mx ? mx : coscnt;
                //ybitCnt->signs->4 mY->dgbitCnt->signs->4 mDg->
                bcw = l > 1 /* c < hwidth */ ? BCW : BCW-1;
                fillbits(rRow, pos, bp, ybcnt, bcw  );
                fillbits(rRow, pos, bp, sycd , SGPCD);
                for(int32_t s=k; s < k+step; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mY, ybcnt_cl);
                }
                //fillbits(rRow, pos, bp, msstv[k  ].mY, ybcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+1].mY, ybcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+2].mY, ybcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+3].mY, ybcnt_cl);
                putgsbits(rRow, pos, bp, sy  , sycd );

                fillbits(rRow, pos, bp, dgbcnt, bcw );
                fillbits(rRow, pos, bp, sdgcd , SGPCD);
                for(int32_t s=k; s < k+step; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mYdDg, dgbcnt_cl);
                }
                //fillbits(rRow, pos, bp, msstv[k  ].mYdDg, dgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+1].mYdDg, dgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+2].mYdDg, dgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+3].mYdDg, dgbcnt_cl);
                putgsbits(rRow, pos, bp, sdg  , sdgcd);
                
                fillbits(rRow, pos, bp, cgbcnt, bcw );
                fillbits(rRow, pos, bp, scgcd , SGPCD);
                for(int32_t s=k; s < k+step; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mCbCg, cgbcnt_cl);
                }
                //fillbits(rRow, pos, bp, msstv[k  ].mCbCg, cgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+1].mCbCg, cgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+2].mCbCg, cgbcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+3].mCbCg, cgbcnt_cl);
                putgsbits(rRow, pos, bp, scg  , scgcd);

                fillbits(rRow, pos, bp, cobcnt, bcw );
                fillbits(rRow, pos, bp, scocd , SGPCD);
                for(int32_t s=k; s < k+step; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mCrCo, cobcnt_cl);
                }
                //fillbits(rRow, pos, bp, msstv[k  ].mCrCo, cobcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+1].mCrCo, cobcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+2].mCrCo, cobcnt_cl);
                //fillbits(rRow, pos, bp, msstv[k+3].mCrCo, cobcnt_cl);
                putgsbits(rRow, pos, bp, sco  , scocd);
                #ifdef _DUMP_
                fprintf(fp,"[k%6d,pos%4d,yc%2u,dgc%2u,cgc%2d,coc%2d],",k, pos, ybcnt, dgbcnt, cgbcnt, cobcnt);
                #endif
            }
        }
        #ifdef _DUMP_
        fprintf(fp ,"\n");
        fprintf(fp2,"\n");
        fprintf(fp3,"\n");
        #endif
        uint32_t p2=0, bp2=0;
        fillbits(rRow, p2, bp2, pos, 32);
        if(r < 2){
            fprintf(stdout, "row %d size: %d \n", r, getRowCmpSize(rRow));
        }
        tsz += (pos-4); 
    }
    fprintf(stdout, "Total Compressed size: %u\n", tsz);
    #ifdef _DUMP_
    fclose(fp );
    fclose(fp2);
    fclose(fp3);
    #endif
}
static inline void retrv(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t& td, uint32_t blen){
    uint8_t b = rRow[pos];
    b >>= bp;
    td = 0;
    for(int32_t m=0; m < blen; ++m){
        td |= (b & 0x1) << m;
        b >>= 1;
        bp++;
        if(bp >=8){
            ++pos;
            b = rRow[pos];
            bp=0;
        }
    }
}
static inline void getgsbits(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t& td, uint32_t blen){
    uint32_t mx = (1U << SGPCD)-1;
    if(blen >= mx){
        retrv(rRow, pos, bp, td, SGP);
    }else{
        retrv(rRow, pos, bp, td, blen);
    }
}
void dec_entropy(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info,
    const ProcessInfo& proc_info, const QuantInfo& quant_info )
{
    //--   +----+----+----+----+----+----+----+----+----+
    //--   | 0  | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 0  |
    //--   +----+----+----+----+----+----+----+----+----+
    //--    <- 5-bit bcnt         -> <- signs of group->
    uint32_t my, mdg, mcg, mco;
    uint32_t ybcnt, dgbcnt, cgbcnt, cobcnt;
    uint32_t sy, sdg, scg, sco;
    uint32_t sycd, sdgcd, scgcd, scocd;
    int32_t k=0, k2=0;
    uint32_t pos, bp;
    int32_t  hwidth = msst_info.w/2;
    uint32_t bcw;
    int32_t  tw;
    #ifdef _DUMP_
    FILE* fp  = fopen("tmp/dec_ent.txt", "w");
    FILE* fp2 = fopen("tmp/dec_out.txt", "w");
    #endif
    uint32_t rsz=0;

    for(int32_t r=0; r < msst_info.h; ++r){
        vector<uint8_t>& rRow = bitv[r];
        pos = 4;
        bp = 0;
        rsz = getRowCmpSize(rRow)+4;
        for(int32_t c=0; pos < rsz; k+=4, c+=4){
            bcw = c < hwidth? BCW: BCW-1;
            retrv(rRow, pos, bp, ybcnt, bcw  );
            ybcnt++; 
            retrv(rRow, pos, bp, sycd , SGPCD);
            retrv(rRow, pos, bp, msstv[k  ].mY, ybcnt);
            retrv(rRow, pos, bp, msstv[k+1].mY, ybcnt);
            retrv(rRow, pos, bp, msstv[k+2].mY, ybcnt);
            retrv(rRow, pos, bp, msstv[k+3].mY, ybcnt);
            //if(k==283728){
            //    int a=10;
            //    int myy ; 
            //    myy = msstv[k].mY;
            //    myy = msstv[k+1].mY ;
            //    myy = msstv[k+2].mY ;
            //    myy = msstv[k+3].mY;
            //}
            getgsbits(rRow, pos, bp, sy, sycd);
            msstv[k  ].sY = sy & 0x1; sy >>= 1;
            msstv[k+1].sY = sy & 0x1; sy >>= 1;
            msstv[k+2].sY = sy & 0x1; sy >>= 1;
            msstv[k+3].sY = sy & 0x1; 

            retrv(rRow, pos, bp, dgbcnt, bcw  );
            dgbcnt++;
            retrv(rRow, pos, bp, sdgcd , SGPCD);
            retrv(rRow, pos, bp, msstv[k  ].mYdDg, dgbcnt);
            retrv(rRow, pos, bp, msstv[k+1].mYdDg, dgbcnt);
            retrv(rRow, pos, bp, msstv[k+2].mYdDg, dgbcnt);
            retrv(rRow, pos, bp, msstv[k+3].mYdDg, dgbcnt);
            getgsbits(rRow, pos, bp, sdg, sdgcd);
            msstv[k  ].sYdDg = sdg & 0x1; sdg >>= 1;
            msstv[k+1].sYdDg = sdg & 0x1; sdg >>= 1;
            msstv[k+2].sYdDg = sdg & 0x1; sdg >>= 1;
            msstv[k+3].sYdDg = sdg & 0x1; 
            
            retrv(rRow, pos, bp, cgbcnt, bcw  );
            cgbcnt++;
            retrv(rRow, pos, bp, scgcd , SGPCD);
            retrv(rRow, pos, bp, msstv[k  ].mCbCg, cgbcnt);
            retrv(rRow, pos, bp, msstv[k+1].mCbCg, cgbcnt);
            retrv(rRow, pos, bp, msstv[k+2].mCbCg, cgbcnt);
            retrv(rRow, pos, bp, msstv[k+3].mCbCg, cgbcnt);
            getgsbits(rRow, pos, bp, scg, scgcd);
            msstv[k  ].sCbCg = scg & 0x1; scg >>= 1;
            msstv[k+1].sCbCg = scg & 0x1; scg >>= 1;
            msstv[k+2].sCbCg = scg & 0x1; scg >>= 1;
            msstv[k+3].sCbCg = scg & 0x1; 

            retrv(rRow, pos, bp, cobcnt, bcw  );
            cobcnt++;
            retrv(rRow, pos, bp, scocd , SGPCD);
            retrv(rRow, pos, bp, msstv[k  ].mCrCo, cobcnt);
            retrv(rRow, pos, bp, msstv[k+1].mCrCo, cobcnt);
            retrv(rRow, pos, bp, msstv[k+2].mCrCo, cobcnt);
            retrv(rRow, pos, bp, msstv[k+3].mCrCo, cobcnt);
            getgsbits(rRow, pos, bp, sco, scocd);
            msstv[k  ].sCrCo = sco & 0x1; sco >>= 1;
            msstv[k+1].sCrCo = sco & 0x1; sco >>= 1;
            msstv[k+2].sCrCo = sco & 0x1; sco >>= 1;
            msstv[k+3].sCrCo = sco & 0x1; 

            //fprintf(stdout," k: %d, ",k);
            #ifdef _DUMP_
            fprintf(fp,"[k%6d,pos%4d,yc%2d,dgc%2d,cgc%2d,coc%2d],",k, pos, ybcnt, dgbcnt, cgbcnt, cobcnt);
            fprintf(fp2,"[k%4d, mY: %6d, %6d, %6d, %6d, sY: %1d, %1d, %1d, %1d]\n",k, msstv[k].mY, msstv[k+1].mY, msstv[k+2].mY, msstv[k+3].mY, msstv[k].sY, msstv[k+1].sY, msstv[k+2].sY, msstv[k+3].sY);
            #endif
        }
        //fprintf(stdout,"\n");
        #ifdef _DUMP_
        fprintf(fp ,"\n");
        fprintf(fp2,"\n");
        #endif
        if(k-k2 != msst_info.w){
            fprintf(stdout,"Row Size: %d\n", k-k2);
            fprintf(stdout, "<-- width not matching\n");
        }
        k2 = k;
    }
    fprintf(stdout, "final k: %d\n", k);
    #ifdef _DUMP_
    fclose(fp);
    fclose(fp2);
    #endif
}
