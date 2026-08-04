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
void write_bitstream(const vector<vector<uint8_t>>& entp, const string& filename){
    FILE *fp;
    fp = fopen(filename.c_str(), "wb");
    if(fp == nullptr){
        fprintf(stderr, "write_bitstream: file open error: %s\n", filename.c_str());
        return;
    }
    size_t h = entp.size();
    size_t w = 0;
    for(size_t m=0; m < h; m++){
        w = entp[m].size();
        for(size_t n=0; n < w; n++){
            fwrite(&entp[m][n], sizeof(uint8_t), 1, fp);
        }
    }
    fclose(fp);
}

void read_bitstream(vector<vector<uint8_t>>& entp, const string& filename){
    FILE *fp;
    fp = fopen(filename.c_str(), "rb");
    if(fp == nullptr){
        fprintf(stderr, "read_bitstream: file open error: %s\n", filename.c_str());
        return;
    }
    size_t h = entp.size();
    size_t w = 0;
    uint32_t rd_sz = 0;
    while(1){
        vector<uint8_t> tenp;
        uint8_t temp[4];
        fread(&rd_sz, sizeof(uint32_t), 1, fp);
        temp[3] = (rd_sz >> 24) & 0xFF;
        temp[2] = (rd_sz >> 16) & 0xFF;
        temp[1] = (rd_sz >> 8 ) & 0xFF;
        temp[0] = (rd_sz      ) & 0xFF;
        tenp.push_back(temp[0]);
        tenp.push_back(temp[1]);
        tenp.push_back(temp[2]);
        tenp.push_back(temp[3]);
        for(size_t n=0; n < rd_sz; n++){
            uint8_t byte;
            fread(&byte, sizeof(uint8_t), 1, fp);
            tenp.push_back(byte);
        }
        entp.push_back(tenp);

        if(feof(fp)){
            break;
        }
    }
    fclose(fp);
}
