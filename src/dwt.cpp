#include "dwt.h"

void fdwt53(vector<msst>& dst, int32_t S, int32_t N ){
    // in-place 5/3 DWT transform
    int32_t hN = N/2;
    int32_t hN_r = hN + (N & 0x1);
    vector<msst> d(hN_r, msst());
    vector<msst> s(hN_r, msst());

    int32_t m = 0, k = S;
    for(; m < hN-1; m++, k+=2){
        d[m] = dst[k+1] - ((dst[k]+dst[k+2])>>1);
    }
    d[m] = dst[k+1] - ((dst[k]+dst[k])>>1); // using symmetric extension for dst[N] = dst[N-2]
   #ifdef DEBUG_DWT
    if(S==0){
        fprintf(stderr, "msst Y of dst[%4d]: %6d ,%6d ,%6d ,%6d \n", S, dst[S].Y, dst[S+1].Y, dst[S+2].Y, dst[S+3].Y);
        fprintf(stderr, "msst Y of dn [%4d]: %6d ,%6d ,%6d ,%6d \n", 0, d[0].Y, d[1].Y, d[2].Y, d[3].Y);
    }
   #endif 
    m = 0; k = S;
    s[m] = dst[k] + ((d[0]+d[0]+2)>>2);  // using symmetric extension 
    m++; k+=2;
    for(; m < hN; m++, k+=2){
        s[m] = dst[k] + ((d[m-1]+d[m]+2)>>2);
    }
   #ifdef DEBUG_DWT
    if(S==0){
        fprintf(stderr, "msst Y of sn [%4d]: %6d ,%6d ,%6d ,%6d \n", 0, s[0].Y, s[1].Y, s[2].Y, s[3].Y);
    }
   #endif

    k=S; m = S+hN;
    for(int32_t n=0; n < hN; k++, m++, n++){
        dst[k] = s[n];
        dst[m] = d[n];
    }
    m = S+hN;
   #ifdef DEBUG_DWT
    if(S==0){
        fprintf(stderr, "msst[%4d]: %6d ,%6d ,%6d ,%6d \n", S, dst[S].Y, dst[S+1].Y, dst[S+2].Y, dst[S+3].Y);
        fprintf(stderr, "msst[%4d]: %6d ,%6d ,%6d ,%6d \n", m, dst[m].Y, dst[m+1].Y, dst[m+2].Y, dst[m+3].Y);
    }
   #endif
}
void idwt53(vector<msst>& dst, int32_t S, int32_t N ){
    // in-place 5/3 inverse DWT transform
    int32_t hN = N/2;
    int32_t hN_r = hN + (N & 0x1);
    vector<msst> d(hN_r, msst());
    vector<msst> s(hN_r, msst());
    int32_t m, k;

    k=S; m = S+hN;
    for(int32_t n=0; n < hN; k++, m++, n++){
        s[n] = dst[k];
        d[n] = dst[m];
    }

    m = 0; k = S;
    dst[k] = s[m] - ((d[0]+d[0]+2)>>2);  // using symmetric extension
    m++; k+=2;
    for(; m < hN; m++, k+=2){
        dst[k] = s[m] - ((d[m-1]+d[m]+2)>>2);
    }

    m= 0; k = S;
    for(; m < hN-1; m++, k+=2){
        dst[k+1] = d[m] + ((dst[k]+dst[k+2])>>1);
    }
    dst[k+1] = d[m] + ((dst[k]+dst[k])>>1); // using symmetric extension for dst[N] = dst[N-2]
   #ifdef DEBUG_DWT
    if(S==0){
        fprintf(stderr, "idwt Y of dst[%4d]: %6d, %6d, %6d, %6d \n", 0, dst[0].Y, dst[1].Y, dst[2].Y, dst[3].Y);
    }
   #endif
}
