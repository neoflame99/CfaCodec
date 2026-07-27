#include "monitor.h"


void savebayertxt(const vector<cfapix> &cfaimg, const string& filename, const BayerInfo& bayer_info){
    FILE *fp;
    fp = fopen(filename.c_str(), "w");
    if(fp == nullptr){
        return;
    }
    size_t numpixels = cfaimg.size();
    if(bayer_info.w * bayer_info.h != numpixels){
        fprintf(stderr, "savebayertxt: image size mismatch, w*h = %d, numpixels = %zu\n", bayer_info.w * bayer_info.h, numpixels);
        fclose(fp);
        return;
    }
    uint32_t i =0;
    for(int32_t m=0; m < bayer_info.h; m++){
        fprintf(fp, "%4d: ", m);
        for(int32_t n=0; n < bayer_info.w; n++){
            fprintf(fp, "%5d ", cfaimg[i]);
            i++;
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void dumpmsst(const vector<msst> &msstv, const string& filename, const MsstInfo& msst_info){
    FILE *fp;
    fp = fopen(filename.c_str(), "w");
    if(fp == nullptr){
        return;
    }
    size_t numpixels = msstv.size();
    if(msst_info.w * msst_info.h != numpixels){
        fprintf(stderr, "dumpmsst: image size mismatch, w*h = %d, numpixels = %zu\n", msst_info.w * msst_info.h, numpixels);
        fclose(fp);
        return;
    }
    uint32_t i =0;
    for(int32_t m=0; m < msst_info.h; m++){
        fprintf(fp, "%4d: ", m);
        for(int32_t n=0; n < msst_info.w; n++){
            fprintf(fp, "(k%6d, %6d ,%6d ,%6d ,%6d) ", i, msstv[i].Y, msstv[i].YdDg, msstv[i].CrCo, msstv[i].CbCg);
            i++;
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}
void dumpmsst(const vector<msstSm> &msstv, const string& filename, const MsstInfo& msst_info){
    size_t len = msstv.size();
    vector<msst> v(len, msst());
    ConvMsstSm2Msst(v, msstv);
    dumpmsst( v, filename, msst_info); 
}