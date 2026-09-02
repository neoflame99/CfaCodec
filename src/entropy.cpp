#include "entropy.h"

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
static inline void putgsbits(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t td, uint32_t blen, uint32_t sgp, uint32_t sgpcd){
    uint32_t mx = (1U << sgpcd)-1;
    if(blen >= sgpcd){
        // blen: 3, 4
        fillbits(rRow, pos, bp, td, sgp);
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
    //--    <- 5-bit bcnt         -> <- signs of group->
    //--     HH band: 5-bit bcnt, others: 5-bit bcnt
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
    int32_t  hwidth = msst_info.w/2;
    int32_t  tw;
    uint32_t bcw;
    #ifdef _DUMP_
    FILE* fp  = fopen("tmp/enc_ent.txt", "w");
    FILE* fp2 = fopen("tmp/enc_in.txt", "w");
    FILE* fp3 = fopen("tmp/qbit.txt", "w");
    fprintf(fp3,"## [sy, sdg, scg, sco]\n");
    #endif

    uint32_t ngrp = proc_info.ngrp;
    uint32_t pbcw = proc_info.bcw;
    uint32_t mbc  = proc_info.mbc;
    uint32_t dwtl = proc_info.dwt_lv;
    uint32_t sgpcd = proc_info.sgpcd;
    uint32_t mx = (1U<<sgpcd)-1;

    for(int32_t r=0; r < msst_info.h; ++r){
        vector<uint8_t>& rRow = bitv[r];
        pos = 4;
        bp  = 0;
        //for(int32_t c=0; c < msst_info.w; c+=ngrp, k+=ngrp){
        for(int32_t l=proc_info.dwt_lv+1; l > 0; --l){
            int32_t sl = l > proc_info.dwt_lv ? proc_info.dwt_lv : l; 
            tw = msst_info.w >> sl;
            for(int32_t c=0; c < tw; c+=ngrp, k+=ngrp){
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
                for(int32_t s=k+1, ss=1; ss < ngrp; ++ss, ++s){
                   my |= msstv[s].mY   ; //msstv[k].mY | msstv[k+1].mY | msstv[k+2].mY | msstv[k+3].mY;
                   mdg|= msstv[s].mYdDg; // msstv[k].mYdDg| msstv[k+1].mYdDg| msstv[k+2].mYdDg| msstv[k+3].mYdDg;
                   mco|= msstv[s].mCrCo; // msstv[k].mCrCo| msstv[k+1].mCrCo| msstv[k+2].mCrCo| msstv[k+3].mCrCo;
                   mcg|= msstv[s].mCbCg; // msstv[k].mCbCg| msstv[k+1].mCbCg| msstv[k+2].mCbCg| msstv[k+3].mCbCg;
                   sy |= int(msstv[s].sY)   << ss; //int(msstv[k].sY) | (int(msstv[k+1].sY) << 1) | (int(msstv[k+2].sY)<<2) | (int(msstv[k+3].sY)<<3);
                   sdg|= int(msstv[s].sYdDg)<< ss; //int(msstv[k].sYdDg) | (int(msstv[k+1].sYdDg) << 1) | (int(msstv[k+2].sYdDg)<<2) | (int(msstv[k+3].sYdDg)<<3);
                   sco|= int(msstv[s].sCrCo)<< ss; //int(msstv[k].sCrCo) | (int(msstv[k+1].sCrCo) << 1) | (int(msstv[k+2].sCrCo)<<2) | (int(msstv[k+3].sCrCo)<<3);
                   scg|= int(msstv[s].sCbCg)<< ss; //int(msstv[k].sCbCg) | (int(msstv[k+1].sCbCg) << 1) | (int(msstv[k+2].sCbCg)<<2) | (int(msstv[k+3].sCbCg)<<3);
                }
                sy2 = sy; sdg2 = sdg; scg2 = scg; sco2 = sco;
                #ifdef _DUMP_
                fprintf(fp2,"[k%4d, mY: %6d, %6d, %6d, %6d, sY: %1d, %1d, %1d, %1d]\n",k, msstv[k].mY, msstv[k+1].mY, msstv[k+2].mY, msstv[k+3].mY, msstv[k].sY, msstv[k+1].sY, msstv[k+2].sY, msstv[k+3].sY);
                ssy = signed_qbit(sy); ssdg = signed_qbit(sdg); sscg = signed_qbit(scg); ssco = signed_qbit(sco);
                int8_t ssy, ssdg, sscg, ssco;
                fprintf(fp3, "[k%6d, %2d, %2d, %2d %2d],", k, ssy, ssdg, sscg, ssco);
                #endif

                for(int32_t m=0; m < mbc; ++m){
                    if( my == 0) break;
                    my >>= 1;
                    ybcnt++;
                }
                for(int32_t m=0; m < mbc; ++m){
                    if( mdg == 0) break;
                    mdg >>= 1;
                    dgbcnt++;
                }
                for(int32_t m=0; m < mbc; ++m){
                    if( mcg == 0) break;
                    mcg >>= 1;
                    cgbcnt++;
                }
                for(int32_t m=0; m < mbc; ++m){
                    if( mco == 0) break;
                    mco >>= 1;
                    cobcnt++;
                }
                for(int32_t m=0; m < ngrp; ++m){
                    if( sy2 == 0) break;
                    sy2 >>= 1;
                    yscnt++;
                }
                for(int32_t m=0; m < ngrp; ++m){
                    if( sdg2 == 0) break;
                    sdg2 >>= 1;
                    dgscnt++;
                }
                for(int32_t m=0; m < ngrp; ++m){
                    if( scg2 == 0) break;
                    scg2 >>= 1;
                    cgscnt++;
                }
                for(int32_t m=0; m < ngrp; ++m){
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
                //ybitCnt->ngrp mY->dgbitCnt->ngrp mDg->
                //sycd ->ngrp sY->sdg->ngrp sDg->
                bcw = pbcw;
                fillbits(rRow, pos, bp, ybcnt, bcw  );
                for(int32_t s=k; s < k+ngrp; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mY, ybcnt_cl);
                }

                fillbits(rRow, pos, bp, dgbcnt, bcw );
                for(int32_t s=k; s < k+ngrp; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mYdDg, dgbcnt_cl);
                }
                
                fillbits(rRow, pos, bp, cgbcnt, bcw );
                for(int32_t s=k; s < k+ngrp; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mCbCg, cgbcnt_cl);
                }

                fillbits(rRow, pos, bp, cobcnt, bcw );
                for(int32_t s=k; s < k+ngrp; ++s){
                    fillbits(rRow, pos, bp, msstv[s  ].mCrCo, cobcnt_cl);
                }
                fillbits (rRow, pos, bp, sycd , sgpcd);
                putgsbits(rRow, pos, bp, sy   , sycd , ngrp, sgpcd);
                fillbits (rRow, pos, bp, sdgcd, sgpcd);
                putgsbits(rRow, pos, bp, sdg  , sdgcd, ngrp, sgpcd);
                fillbits (rRow, pos, bp, scgcd, sgpcd);
                putgsbits(rRow, pos, bp, scg  , scgcd, ngrp, sgpcd);
                fillbits (rRow, pos, bp, scocd, sgpcd);
                putgsbits(rRow, pos, bp, sco  , scocd, ngrp, sgpcd);
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
        if(r < 2 || r == 268){
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
static inline void getgsbits(vector<uint8_t>& rRow, uint32_t& pos, uint32_t& bp, uint32_t& td, uint32_t blen, uint32_t sgp, uint32_t sgpcd){
    uint32_t mx = (1U << sgpcd)-1;
    if(blen >= mx){
        retrv(rRow, pos, bp, td, sgp);
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

    uint32_t ngrp = proc_info.ngrp;
    uint32_t pbcw = proc_info.bcw;
    uint32_t mbc  = proc_info.mbc;
    uint32_t dwtl = proc_info.dwt_lv;
    uint32_t sgpcd = proc_info.sgpcd;

    for(int32_t r=0; r < msst_info.h; ++r){
        vector<uint8_t>& rRow = bitv[r];
        pos = 4;
        bp = 0;
        rsz = getRowCmpSize(rRow)+4;
        for(int32_t c=0; pos < rsz; k+=ngrp, c+=ngrp){
            //bcw = c < hwidth? pbcw : pbcw-1;
            bcw = pbcw ;
            retrv(rRow, pos, bp, ybcnt, bcw  );
            ybcnt++; 
            for(int32_t s=k; s < k+ngrp; ++s){
                retrv(rRow, pos, bp, msstv[s  ].mY, ybcnt);
            }
            retrv(rRow, pos, bp, dgbcnt, bcw  );
            dgbcnt++;
            for(int32_t s=k; s < k+ngrp; ++s){
                retrv(rRow, pos, bp, msstv[s  ].mYdDg, dgbcnt);
            }
            retrv(rRow, pos, bp, cgbcnt, bcw  );
            cgbcnt++;
            for(int32_t s=k; s < k+ngrp; ++s){
                retrv(rRow, pos, bp, msstv[s  ].mCbCg, cgbcnt);
            }   
            retrv(rRow, pos, bp, cobcnt, bcw  );
            cobcnt++;
            for(int32_t s=k; s < k+ngrp; ++s){
                retrv(rRow, pos, bp, msstv[s  ].mCrCo, cobcnt);
            }

            retrv(rRow, pos, bp, sycd , sgpcd);
            getgsbits(rRow, pos, bp, sy, sycd, ngrp, sgpcd);
            for(int32_t s=k; s < k+ngrp; ++s){
                msstv[s].sY = sy & 0x1; sy >>= 1;
            }
            retrv(rRow, pos, bp, sdgcd , sgpcd);
            getgsbits(rRow, pos, bp, sdg, sdgcd, ngrp, sgpcd);
            for(int32_t s=k; s < k+ngrp; ++s){
                msstv[s].sYdDg = sdg & 0x1; sdg >>= 1;
            }
            retrv(rRow, pos, bp, scgcd , sgpcd);
            getgsbits(rRow, pos, bp, scg, scgcd, ngrp, sgpcd);
            for(int32_t s=k; s < k+ngrp; ++s){
                msstv[s].sCbCg = scg & 0x1; scg >>= 1;
            }
            retrv(rRow, pos, bp, scocd , sgpcd);
            getgsbits(rRow, pos, bp, sco, scocd, ngrp, sgpcd);
            for(int32_t s=k; s < k+ngrp; ++s){
                msstv[s].sCrCo = sco & 0x1; sco >>= 1;
            }


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
            fprintf(stdout,"%d Row Size: %d\n", r, k-k2);
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

struct GR{
    uint32_t q;
    uint32_t r;
    GR():q(0), r(0){}
};

static inline void Conv_golombrice(GR& gr, int32_t s,const uint8_t K){
    gr.q = s >> K;
    gr.r = s &((1U << K)-1);
}
static inline void Revr_golombrice(int32_t& s, const GR& gr, const uint8_t K){
    s = gr.q << K;
    s = s | gr.r;
}
static inline void enc_gmrice(int32_t dat, vector<uint8_t>& bstm, uint32_t& pos, uint32_t& bp, uint8_t K){
    GR gr;
    Conv_golombrice(gr, dat, K);
    uint8_t q = 0, r=0;
    uint8_t b;
    uint8_t d = bstm[pos];
    // put q
    for(int n=0; n < gr.q; ++n){
        d |= 0x1 << bp; 
        bp++;
        if(bp >= 8){
            bstm[pos] = d;
            ++pos;
            bp= 0;
            d = 0;
        }
    }
    // put 0
    bp++;
    if(bp >= 8){
        bstm[pos] = d;
        ++pos;
        bp= 0;
        d = 0;
    }
    // put r
    for(int k =0; k < K; ++k){
        d |= (r & 0x1) << k;
        r >>=1;
        bp++;
        if(bp >= 8){
            bstm[pos] = d;
            ++pos;
            bp= 0;
            d = 0;
        }
    }
    if(bp!=0){
        bstm[pos] = d;
    }
}
static inline void dec_gmrice(int32_t& dat, const vector<uint8_t>& bstm, uint32_t& pos, uint32_t& bp, uint8_t K){
    GR gr;
    uint8_t d = bstm[pos];
    d >>= bp;
    uint32_t q = 0, r=0;
    uint32_t b;
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
    gr.q = q;
    gr.r = r;
    Revr_golombrice(dat, gr, K);
}

void enc_entropy_gr(vector<vector<uint8_t>>& bitv, vector<msstSm>& msstv, const MsstInfo& msst_info, 
    const ProcessInfo& proc_info, const QuantInfo& quant_info )
{
    //--   +----+----+----+----+----+----+----+----+----+
    //--   | 0  | 1  | 2  | 3  | 4  | 5  | 6  | 7  | 0  |
    //--   +----+----+----+----+----+----+----+----+----+
    //--    <- 5-bit bcnt         -> <- signs of group->
    //--     HH band: 5-bit bcnt, others: 5-bit bcnt
    uint32_t my, mdg, mcg, mco;
    uint32_t yscnt, dgscnt, cgscnt, coscnt;
    uint32_t sy, sdg, scg, sco;
    uint32_t sy2, sdg2, scg2, sco2;
    uint32_t sycd, sdgcd, scgcd, scocd;
    int32_t k=0;
    uint32_t pos, bp;
    uint32_t tsz=0;
    int32_t  hwidth = msst_info.w/2;
    int32_t  tw;
    #ifdef _DUMP2_
    FILE* fp  = fopen("tmp/enc_ent.txt", "w");
    FILE* fp2 = fopen("tmp/enc_in.txt", "w");
    FILE* fp3 = fopen("tmp/qbit.txt", "w");
    fprintf(fp3,"## [sy, sdg, scg, sco]\n");
    #endif

    uint32_t ngrp = proc_info.ngrp;
    uint32_t dwtl = proc_info.dwt_lv;
    uint32_t sgpcd = proc_info.sgpcd;
    uint32_t mx = (1U<<sgpcd)-1;

    for(int32_t r=0; r < msst_info.h; ++r){
        vector<uint8_t>& rRow = bitv[r];
        pos = 4;
        bp  = 0;
        //for(int32_t c=0; c < msst_info.w; c+=ngrp, k+=ngrp){
        for(int32_t l=proc_info.dwt_lv+1; l > 0; --l){
            int32_t sl = l > proc_info.dwt_lv ? proc_info.dwt_lv : l; 
            tw = msst_info.w >> sl;
            for(int32_t c=0; c < tw; c+=ngrp, k+=ngrp){

                //--------------------------------------------------------
                // ngrp mY -> ngrp mDg -> ngrp mCg -> ngrp mCo ->
                // signs mY -> signs mDg -> signs mCg -> signs mCo
                //--------------------------------------------------------

                //-- find optimum K for Golomb-Rice coding
                uint32_t mnpos = 0xFFFFFFFF;
                uint32_t mnbp  = 0xFFFFFFFF;
                uint32_t K = 1;
                for(int32_t t=4; t <=16; ++t){
                    vector<uint8_t> Kvec(ngrp*4*2, 0);
                    uint32_t pos2 = 0;
                    uint32_t bp2 = 0;
                    for(int32_t n=0, m=k; n < ngrp; ++n, ++m){
                        enc_gmrice(msstv[m].mY   , Kvec, pos2, bp2, t);
                        enc_gmrice(msstv[m].mYdDg, Kvec, pos2, bp2, t);
                        enc_gmrice(msstv[m].mCbCg, Kvec, pos2, bp2, t);
                        enc_gmrice(msstv[m].mCrCo, Kvec, pos2, bp2, t);
                    }
                    if((pos2 < mnpos ) || 
                       (pos2 == mnpos && bp2 < mnbp)){
                        mnpos = pos2;
                        mnbp  = bp2;
                        K = t;
                    }
                }
                //-- encode with optimum K
                fillbits(rRow, pos, bp, K-1, 4); // 4-bit K(1~16 -> 0~15)
                for(int32_t n=0, m=k; n < ngrp; ++n, ++m){
                    enc_gmrice(msstv[m].mY   , rRow, pos, bp, K);
                    enc_gmrice(msstv[m].mYdDg, rRow, pos, bp, K);
                    enc_gmrice(msstv[m].mCbCg, rRow, pos, bp, K);
                    enc_gmrice(msstv[m].mCrCo, rRow, pos, bp, K);
                }

                assert(ngrp >= NSGRP && ngrp % NSGRP == 0);
                int32_t iter = ngrp / NSGRP;
                for(int32_t i=0, m=k; i < iter; ++i, m += NSGRP){
                    yscnt = 0;
                    dgscnt= 0;
                    cgscnt= 0;
                    coscnt= 0;
                    sy = int(msstv[m].sY)   ; 
                    sdg= int(msstv[m].sYdDg); 
                    sco= int(msstv[m].sCrCo); 
                    scg= int(msstv[m].sCbCg); 
                    for(int32_t s=m+1, ss=1; ss < NSGRP; ++ss, ++s){
                       sy |= int(msstv[s].sY)   << ss; 
                       sdg|= int(msstv[s].sYdDg)<< ss; 
                       sco|= int(msstv[s].sCrCo)<< ss; 
                       scg|= int(msstv[s].sCbCg)<< ss; 
                    }
                    sy2 = sy; sdg2 = sdg; scg2 = scg; sco2 = sco;

                    for(int32_t s=0; s < NSGRP; ++s){
                        yscnt  += int(sy2 > 0);
                        dgscnt += int(sdg2 > 0);
                        cgscnt += int(scg2 > 0);
                        coscnt += int(sco2 > 0);
                        sy2  >>= 1;
                        sdg2 >>= 1;
                        scg2 >>= 1;
                        sco2 >>= 1;
                    }
                    sycd = yscnt >=mx ? mx : yscnt;
                    sdgcd= dgscnt>=mx ? mx : dgscnt;
                    scgcd= cgscnt>=mx ? mx : cgscnt;
                    scocd= coscnt>=mx ? mx : coscnt;

                    fillbits(rRow, pos, bp, sycd , sgpcd);
                    putgsbits(rRow, pos, bp, sy  , sycd, NSGRP, sgpcd);

                    fillbits(rRow, pos, bp, sdgcd , sgpcd);
                    putgsbits(rRow, pos, bp, sdg  , sdgcd, NSGRP, sgpcd);

                    fillbits(rRow, pos, bp, scgcd , sgpcd);
                    putgsbits(rRow, pos, bp, scg  , scgcd, NSGRP, sgpcd);

                    fillbits(rRow, pos, bp, scocd , sgpcd);
                    putgsbits(rRow, pos, bp, sco  , scocd, NSGRP, sgpcd);
                }
            }
        }
        #ifdef _DUMP2_
        fprintf(fp ,"\n");
        fprintf(fp2,"\n");
        fprintf(fp3,"\n");
        #endif
        uint32_t p2=0, bp2=0;
        fillbits(rRow, p2, bp2, pos, 32);
        if(r < 2 || r == 268){
            fprintf(stdout, "row %d size: %d \n", r, getRowCmpSize(rRow));
        }
        tsz += (pos-4); 
    }
    fprintf(stdout, "Total Compressed size: %u\n", tsz);
    #ifdef _DUMP2_
    fclose(fp );
    fclose(fp2);
    fclose(fp3);
    #endif
}