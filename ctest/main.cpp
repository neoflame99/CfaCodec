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
    SaveInfo save_info;

    float g  = 1.f; //1.0f; // for NLT
    uint32_t dwt_l= DWTL;
    uint32_t ngrp = 8;
    uint32_t bcw  = BCW;
    uint32_t mbc  = MBC;
    uint32_t sgpcd= SGPCD;
    bool sel_ycc  = false;
    ProcessInfo process_info{ g, dwt_l, ngrp, bcw, mbc, sgpcd, sel_ycc };
    bool loadparam_result = loadBayerImgInfo(bayer_img_info, save_info, process_info, param_file);
    if(!loadparam_result){
        cout << "Failed to load parameter file: " << param_file << endl;
        return -1;
    }

    process_info.mbc = bayer_img_info.bpp + process_info.dwt_lv;

    vector<cfapix> cfaorg ;
    vector<cfapix> cfaproc(bayer_img_info.w * bayer_img_info.h, 0);
    vector<cfapix> cfaimg (bayer_img_info.w * bayer_img_info.h, 0);
    vector<msstSm> msstv_enc(bayer_img_info.w * bayer_img_info.h /4, msst());
    vector<msstSm> msstv_dec(bayer_img_info.w * bayer_img_info.h /4, msst());
    loadbayerimg(cfaorg, bayer_img_info.filename, bayer_img_info.bpp);
    if(cfaorg.empty()){
        cout << "Failed to load Bayer image: " << bayer_img_info.filename << endl;
        return -1;
    }else{
        cout << "Successfully loaded Bayer image: " << bayer_img_info.filename << endl;
        cout << "Image Total Size: " << bayer_img_info.w * bayer_img_info.h << endl;
        cout << "cfaimg's size: " << cfaorg.size() << endl;
        size_t cfaorg_size = cfaorg.size();
        for(size_t i=0; i < cfaorg_size; i++){
            cfaimg[i] = cfaorg[i];
        }
    }

    BayerInfo bayer_info{bayer_img_info.w, bayer_img_info.h, 0, bayer_img_info.bpp};
    MsstInfo msst_info{bayer_img_info.w/2, bayer_img_info.h/2};
    QuantInfo quant_info;
    quant_info.Qp = 5;
    quant_info.Rp = 0;
    int32_t NB = 1U << process_info.dwt_lv; 
    for(int k=0; k < NB; ++k){
        if(k==0){
            quant_info.Gb[k] = quant_info.Qp;//-4; //-k; 
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


    //=========== Encoding ======================//
    bayer_fdwt(msstv_enc, cfaimg, bayer_info, process_info, quant_info);
    enc_entropy(entp, msstv_enc, msst_info, process_info, quant_info);
    //===========================================//

    //--------Dump Enc Side----------------------//
    if(save_info.save_cfaimg){
        savebayertxt(cfaorg, "cfaimg.txt", bayer_info);
    }
    if(save_info.save_msstv_enc){
        dumpmsst(msstv_enc, "msstv_enc.txt", msst_info);
    }
    //-------------------------------------------//
    
    write_bitstream(entp, save_info.bitstream_filename);

    //=========== Decoding ======================//
    dec_entropy(entp, msstv_dec, msst_info, process_info, quant_info);
    bayer_idwt(cfaproc, msstv_dec, bayer_info, process_info, quant_info);
    //===========================================//

    //--------Dump Dec Side----------------------//
    if(save_info.save_msstv_dec){
        dumpmsst(msstv_dec, "msstv_dec.txt", msst_info);
    }
    if(save_info.save_cfaproc){
        savebayertxt(cfaproc, "cfaproc.txt", bayer_info);
    }
    //-------------------------------------------//

    size_t sz = cfaorg.size();
    int32_t cnt=0;
    for(size_t i=0; i < sz; i++){
        if(cfaorg[i] != cfaproc[i]){
            cout << "Mismatch at index " << i << ": cfaorg = " << cfaorg[i] << ", cfaproc = " << cfaproc[i] << endl;
            cnt++;
            if(cnt >= 10) break;
        }
    }

    return 0;
}
