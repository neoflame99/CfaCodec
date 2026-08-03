#include <stdio.h>
#include <iostream>
#include "common.h"
#include "msst.h"
#include "dwt.h"
#include "imgloader.h"
#include "loadparam.h"
#include "bayer_dwt.h"
#include "monitor.h"
#include "entropy.h"


using namespace std;
int main(int args, char* argv[]){


    if(args < 2){
        cout << "Usage: " << argv[0] << " param_file" << endl;
        return -1;
    }

    string param_file = argv[1];
    BayerImgInfo bayer_img_info;
    bool loadparam_result = loadBayerImgInfo(bayer_img_info, param_file);
    if(!loadparam_result){
        cout << "Failed to load parameter file: " << param_file << endl;
        return -1;
    }

    vector<cfapix> cfaimg;
    vector<cfapix> cfaproc(bayer_img_info.w * bayer_img_info.h, 0);
    vector<msstSm> msstv_enc(bayer_img_info.w * bayer_img_info.h /4, msst());
    vector<msstSm> msstv_dec(bayer_img_info.w * bayer_img_info.h /4, msst());
    loadbayerimg(cfaimg, bayer_img_info.filename, bayer_img_info.bpp);
    if(cfaimg.empty()){
        cout << "Failed to load Bayer image: " << bayer_img_info.filename << endl;
        return -1;
    }else{
        cout << "Successfully loaded Bayer image: " << bayer_img_info.filename << endl;
        cout << "Image Total Size: " << bayer_img_info.w * bayer_img_info.h << endl;
        cout << "cfaimg's size: " << cfaimg.size() << endl;
    }

    constexpr int32_t dwt_l= DWTL;
    constexpr int32_t ngrp = 8;
    ProcessInfo process_info{1.0f, dwt_l, ngrp};
    BayerInfo bayer_info{bayer_img_info.w, bayer_img_info.h, 0, bayer_img_info.bpp};
    MsstInfo msst_info{bayer_img_info.w/2, bayer_img_info.h/2};
    QuantInfo quant_info;
    quant_info.Qp = 5;
    quant_info.Rp = 0;
    int32_t NB = 1U << process_info.dwt_lv; 
    for(int k=0; k < NB; ++k){
        if(k==0){
            quant_info.Gb[k] = quant_info.Qp-4; //-k; 
        }else{
            quant_info.Gb[k] = quant_info.Qp;
        }
        quant_info.Pb[k] = 0;
    }
    quant_info.getTb();

    vector<vector<uint8_t>> entp;
    for(int32_t k=0; k < msst_info.h; ++k){
        vector<uint8_t> rRow(msst_info.w*4*2, 0);
        entp.push_back(rRow);
    }

    savebayertxt(cfaimg, "cfaimg.txt", bayer_info);

    bayer_fdwt(msstv_enc, cfaimg, bayer_info, process_info, quant_info);
    enc_entropy(entp, msstv_enc, msst_info, process_info, quant_info);
    dumpmsst(msstv_enc, "msstv_enc.txt", msst_info);
    //--
    dec_entropy(entp, msstv_dec, msst_info, process_info, quant_info);
    bayer_idwt(cfaproc, msstv_dec, bayer_info, process_info, quant_info);

    dumpmsst(msstv_dec, "msstv_dec.txt", msst_info);
    savebayertxt(cfaproc, "cfaproc.txt", bayer_info);

    size_t sz = cfaimg.size();
    int32_t cnt=0;
    for(size_t i=0; i < sz; i++){
        if(cfaimg[i] != cfaproc[i]){
            cout << "Mismatch at index " << i << ": cfaimg = " << cfaimg[i] << ", cfaproc = " << cfaproc[i] << endl;
            cnt++;
            if(cnt >= 10) break;
        }
    }


    return 0;
}