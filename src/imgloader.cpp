#include "imgloader.h"

void loadbayerimg(vector<cfapix> &cfaimg, const string filename, const int bpp){
    FILE *fp;
    fp = fopen(filename.c_str(), "rb");
    if(fp ==nullptr){
        return;
    }
    fseek(fp, 0, SEEK_END);
    size_t filesize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    uint8_t *imbuf = (uint8_t *)malloc(filesize);
    size_t readsize = fread(imbuf, 1, filesize, fp);
    if(readsize != filesize){
        free(imbuf);
        fclose(fp);
        return ;
    }
    size_t numpixels = filesize / bpp;
    int32_t *imbuf32 = (int32_t *)malloc(numpixels * sizeof(int32_t));
    if(imbuf32 == nullptr){
        free(imbuf);
        fclose(fp);
        return ;
    }
    size_t pp = 0, bsz=0, pidx=0;
    int32_t dat = 0;
    while( pp < filesize ){
        //process for bayer image
        dat |= ((int32_t)imbuf[pp]) << (8*bsz);
        if( (++bsz) >= bpp){
            imbuf32[pidx++] = dat;
            bsz = 0;
            dat = 0;
        }
        pp++; 
    }
    for(size_t i=0; i < numpixels; i++){
        cfaimg.push_back(imbuf32[i]);
    }
    free(imbuf);
    fclose(fp);
}
