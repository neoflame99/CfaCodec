#include "monitor.h"

namespace {

static void write_le16(FILE* fp, uint16_t value){
    uint8_t bytes[2] = {
        static_cast<uint8_t>(value & 0xFF),
        static_cast<uint8_t>((value >> 8) & 0xFF)
    };
    fwrite(bytes, sizeof(uint8_t), 2, fp);
}

static void write_le32(FILE* fp, uint32_t value){
    uint8_t bytes[4] = {
        static_cast<uint8_t>( value        & 0xFF),
        static_cast<uint8_t>((value >> 8 ) & 0xFF),
        static_cast<uint8_t>((value >> 16) & 0xFF),
        static_cast<uint8_t>((value >> 24) & 0xFF)
    };
    fwrite(bytes, sizeof(uint8_t), 4, fp);
}

static uint8_t cfapix_upper8(cfapix value, int32_t bpp){
    if(value < 0){
        value = 0;
    }

    const uint32_t raw = static_cast<uint32_t>(value);
    if(bpp > 8){
        return static_cast<uint8_t>(raw >> (bpp - 8));
    }
    if(bpp > 0){
        return static_cast<uint8_t>(raw << (8 - bpp));
    }
    return 0;
}

static int bayer_channel(int32_t pattern, int32_t y, int32_t x){
    const bool odd_y = (y & 1) != 0;
    const bool odd_x = (x & 1) != 0;

    switch(pattern){
        case 0: // rggb
            if(!odd_y && !odd_x) return 0;
            if( odd_y &&  odd_x) return 2;
            return 1;
        case 1: // grbg
            if(!odd_y &&  odd_x) return 0;
            if( odd_y && !odd_x) return 2;
            return 1;
        case 2: // gbrg
            if( odd_y && !odd_x) return 0;
            if(!odd_y &&  odd_x) return 2;
            return 1;
        case 3: // bggr
            if( odd_y &&  odd_x) return 0;
            if(!odd_y && !odd_x) return 2;
            return 1;
        default:
            return 1;
    }
}

} // namespace

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

void savebayerbmp(const vector<cfapix> &cfaimg, const string& filename, const BayerInfo& bayer_info){
    const size_t numpixels = cfaimg.size();
    const size_t expected_pixels = static_cast<size_t>(bayer_info.w) * static_cast<size_t>(bayer_info.h);
    if(expected_pixels != numpixels){
        fprintf(stderr, "savebayerbmp: image size mismatch, w*h = %zu, numpixels = %zu\n", expected_pixels, numpixels);
        return;
    }
    if(bayer_info.w <= 0 || bayer_info.h <= 0){
        fprintf(stderr, "savebayerbmp: invalid image size, w = %d, h = %d\n", bayer_info.w, bayer_info.h);
        return;
    }

    FILE *fp = fopen(filename.c_str(), "wb");
    if(fp == nullptr){
        fprintf(stderr, "savebayerbmp: file open error: %s\n", filename.c_str());
        return;
    }

    const uint32_t width = static_cast<uint32_t>(bayer_info.w);
    const uint32_t height = static_cast<uint32_t>(bayer_info.h);
    const uint32_t row_bytes = width * 3;
    const uint32_t padding = (4 - (row_bytes % 4)) % 4;
    const uint32_t stride = row_bytes + padding;
    const uint32_t pixel_bytes = stride * height;
    const uint32_t header_bytes = 54;
    const uint32_t file_bytes = header_bytes + pixel_bytes;

    fwrite("BM", sizeof(char), 2, fp);
    write_le32(fp, file_bytes);
    write_le16(fp, 0);
    write_le16(fp, 0);
    write_le32(fp, header_bytes);

    write_le32(fp, 40);
    write_le32(fp, width);
    write_le32(fp, height);
    write_le16(fp, 1);
    write_le16(fp, 24);
    write_le32(fp, 0);
    write_le32(fp, pixel_bytes);
    write_le32(fp, 2835);
    write_le32(fp, 2835);
    write_le32(fp, 0);
    write_le32(fp, 0);

    const uint8_t pad[3] = {0, 0, 0};
    for(int32_t y = bayer_info.h - 1; y >= 0; --y){
        for(int32_t x = 0; x < bayer_info.w; ++x){
            const size_t index = static_cast<size_t>(y) * static_cast<size_t>(bayer_info.w) + static_cast<size_t>(x);
            const uint8_t value = cfapix_upper8(cfaimg[index], bayer_info.bpp);
            uint8_t bgr[3] = {0, 0, 0};
            const int channel = bayer_channel(bayer_info.cfa_pat, y, x);
            if(channel == 0){
                bgr[2] = value;
            }else if(channel == 1){
                bgr[1] = value;
            }else{
                bgr[0] = value;
            }
            fwrite(bgr, sizeof(uint8_t), 3, fp);
        }
        fwrite(pad, sizeof(uint8_t), padding, fp);
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
    size_t w;
    for(size_t m=0; m < h; m++){
        w = entp[m][0];
        w|= entp[m][1]<< 8;
        w|= entp[m][2]<<16;
        w|= entp[m][3]<<24;
        for(size_t n=0; n < w+4; n++){
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
static void write_raw8(FILE* fp, const vector<cfapix>& cfa){
    size_t n = cfa.size();
    vector<uint8_t> buf(n);
    for(size_t i=0; i < n; ++i){
        buf[i] = (uint8_t)(cfa[i] & 0xFF);
    }
    fwrite(buf.data(), sizeof(uint8_t), buf.size(), fp);
}

static void write_raw10_csi2(FILE* fp, const vector<cfapix>& cfa){
    // 4 pixels (10-bit each) packed into 5 bytes
    size_t n = cfa.size();
    size_t ngrp = n/4;
    vector<uint8_t> buf(ngrp*5);
    size_t bi = 0;
    for(size_t g=0; g < ngrp; ++g){
        uint32_t p0 = (uint32_t)cfa[g*4  ] & 0x3FF;
        uint32_t p1 = (uint32_t)cfa[g*4+1] & 0x3FF;
        uint32_t p2 = (uint32_t)cfa[g*4+2] & 0x3FF;
        uint32_t p3 = (uint32_t)cfa[g*4+3] & 0x3FF;
        buf[bi++] = (uint8_t)(p0 >> 2);
        buf[bi++] = (uint8_t)(p1 >> 2);
        buf[bi++] = (uint8_t)(p2 >> 2);
        buf[bi++] = (uint8_t)(p3 >> 2);
        buf[bi++] = (uint8_t)((p0 & 0x3) | ((p1 & 0x3)<<2) | ((p2 & 0x3)<<4) | ((p3 & 0x3)<<6));
    }
    fwrite(buf.data(), sizeof(uint8_t), buf.size(), fp);
    if(n % 4 != 0){
        fprintf(stderr, "write_dec_cfa: RAW10 packing needs a pixel count multiple of 4, %zu leftover pixel(s) dropped\n", n % 4);
    }
}

static void write_raw12_csi2(FILE* fp, const vector<cfapix>& cfa){
    // 2 pixels (12-bit each) packed into 3 bytes
    size_t n = cfa.size();
    size_t ngrp = n/2;
    vector<uint8_t> buf(ngrp*3);
    size_t bi = 0;
    for(size_t g=0; g < ngrp; ++g){
        uint32_t p0 = (uint32_t)cfa[g*2  ] & 0xFFF;
        uint32_t p1 = (uint32_t)cfa[g*2+1] & 0xFFF;
        buf[bi++] = (uint8_t)(p0 >> 4);
        buf[bi++] = (uint8_t)(p1 >> 4);
        buf[bi++] = (uint8_t)((p0 & 0xF) | ((p1 & 0xF)<<4));
    }
    fwrite(buf.data(), sizeof(uint8_t), buf.size(), fp);
    if(n % 2 != 0){
        fprintf(stderr, "write_dec_cfa: RAW12 packing needs a pixel count multiple of 2, %zu leftover pixel(s) dropped\n", n % 2);
    }
}

static void write_raw14_csi2(FILE* fp, const vector<cfapix>& cfa){
    // 4 pixels (14-bit each) packed into 7 bytes
    size_t n = cfa.size();
    size_t ngrp = n/4;
    vector<uint8_t> buf(ngrp*7);
    size_t bi = 0;
    for(size_t g=0; g < ngrp; ++g){
        uint32_t p0 = (uint32_t)cfa[g*4  ] & 0x3FFF;
        uint32_t p1 = (uint32_t)cfa[g*4+1] & 0x3FFF;
        uint32_t p2 = (uint32_t)cfa[g*4+2] & 0x3FFF;
        uint32_t p3 = (uint32_t)cfa[g*4+3] & 0x3FFF;
        uint32_t low = (p0 & 0x3F) | ((p1 & 0x3F)<<6) | ((p2 & 0x3F)<<12) | ((p3 & 0x3F)<<18);
        buf[bi++] = (uint8_t)(p0 >> 6);
        buf[bi++] = (uint8_t)(p1 >> 6);
        buf[bi++] = (uint8_t)(p2 >> 6);
        buf[bi++] = (uint8_t)(p3 >> 6);
        buf[bi++] = (uint8_t)(low        & 0xFF);
        buf[bi++] = (uint8_t)((low >> 8) & 0xFF);
        buf[bi++] = (uint8_t)((low >> 16)& 0xFF);
    }
    fwrite(buf.data(), sizeof(uint8_t), buf.size(), fp);
    if(n % 4 != 0){
        fprintf(stderr, "write_dec_cfa: RAW14 packing needs a pixel count multiple of 4, %zu leftover pixel(s) dropped\n", n % 4);
    }
}

static void write_compact(FILE* fp, const vector<cfapix>& cfa, int32_t bpp){
    // continuous bpp-bit bitstream, LSB-first, no per-pixel padding/expansion
    size_t n = cfa.size();
    size_t total_bits = n * (size_t)bpp;
    vector<uint8_t> buf((total_bits+7)/8, 0);
    size_t pos = 0, bp = 0;
    uint32_t mask = (bpp >= 32) ? 0xFFFFFFFFU : ((1U << bpp) - 1);
    //fprintf(stdout, "bayer_val[888962] = %d\n", cfa[888962]);
    for(size_t i=0; i < n; ++i){
        uint32_t v = (uint32_t)cfa[i] & mask;
        //if(i == 888962){
        //    fprintf(stdout, "write_compact: i=%zu, v=%d\n", i, v);
        //}
        for(int32_t b=0; b < bpp; ++b){
            buf[pos] |= (uint8_t)(((v >> b) & 0x1) << bp);
            bp++;
            if(bp >= 8){
                bp = 0;
                pos++;
            }
        }
    }
    uint32_t d = buf[888962*2+1];
    d <<= 8;
    d |= buf[888962*2];

    //fprintf(stdout, "bayer_val[888962] = %d\n", d);
    fwrite(buf.data(), sizeof(uint8_t), buf.size(), fp);
}

void write_dec_cfa(const string& fname, const vector<cfapix>& cfa, const BayerInfo& bayer_info, const bool csi2_style){
    size_t numpixels = cfa.size();
    if((size_t)(bayer_info.w * bayer_info.h) != numpixels){
        fprintf(stderr, "write_dec_cfa: image size mismatch, w*h = %d, numpixels = %zu\n", bayer_info.w * bayer_info.h, numpixels);
        return;
    }

    FILE *fp = fopen(fname.c_str(), "wb");
    if(fp == nullptr){
        fprintf(stderr, "write_dec_cfa: file open error: %s\n", fname.c_str());
        return;
    }

    if(csi2_style){
        switch(bayer_info.bpp){
            case 8:
                write_raw8(fp, cfa);
                break;
            case 10:
                write_raw10_csi2(fp, cfa);
                break;
            case 12:
                write_raw12_csi2(fp, cfa);
                break;
            case 14:
                write_raw14_csi2(fp, cfa);
                break;
            case 16:
                write_compact(fp, cfa, 16); // CSI2 RAW16 has no sub-byte packing
                break;
            default:
                fprintf(stderr, "write_dec_cfa: CSI2 packing not defined for bpp=%d, falling back to compact storage\n", bayer_info.bpp);
                write_compact(fp, cfa, bayer_info.bpp);
                break;
        }
    }else{
        write_compact(fp, cfa, bayer_info.bpp);
    }

    fclose(fp);
}

void write_report(const string& fname,const BayerInfo& bayer_info, 
        ProcessInfo& proc_info, SaveInfo& save_info, 
        QuantInfo& quant_info, size_t tsz, const string& InputBayerFname)
{
    FILE *fp = fopen(fname.c_str(), "w");
    if(fp == nullptr){
        fprintf(stderr, "write_report: file open error: %s\n", fname.c_str());
        return ;
    }
    fprintf(fp,"CFA Source: %s\n", InputBayerFname.c_str());
    fprintf(fp, "Bayer Info:\n");
    fprintf(fp, "    width  : %5d\n", bayer_info.w);
    fprintf(fp, "    height : %5d\n", bayer_info.h);
    fprintf(fp, "    pattern: %5d //0: rggb, 1: grbg, 2: gbrg, 3: bggr\n", bayer_info.cfa_pat);
    fprintf(fp, "    bpp    : %5d\n", bayer_info.bpp);
    fprintf(fp, "Save Info:\n");
    fprintf(fp, "    bitstream size : %zu bytes\n", tsz);
    fprintf(fp, "    bitstream_fname: %s\n", save_info.bitstream_fname.c_str());
    fprintf(fp, "    dec_cfa_fname  : %s\n", save_info.dec_cfa_fname.c_str());
    fprintf(fp, "    dec_cfa_csi2_style: %d\n", save_info.dec_cfa_csi2_style);
    fprintf(fp, "Process Info:\n");
    fprintf(fp, "    g       : %5.2f\n", proc_info.g);
    fprintf(fp, "    dwt_lv  : %5u\n"  , proc_info.dwt_lv);
    fprintf(fp, "    ngrp    : %5u\n"  , proc_info.ngrp);
    fprintf(fp, "    mbc     : %5u\n"  , proc_info.mbc);
    fprintf(fp, "    sgpcd   : %5u\n"  , proc_info.sgpcd);
    fprintf(fp, "    sel_ycc : %5u\n"  , proc_info.sel_ycc);
    fprintf(fp, "Quant Info:\n");
    fprintf(fp, "    Qp   : %5d\n"  , quant_info.Qp);
    fprintf(fp, "    Rp   : %5d\n"  , quant_info.Rp);
    for(int k=0; k <= proc_info.dwt_lv; ++k){
        fprintf(fp, "    Pb[%d]: %5d\n"  ,k, quant_info.Pb[k]);
        fprintf(fp, "    Gb[%d]: %5d\n"  ,k, quant_info.Gb[k]);
        fprintf(fp, "    Tb[%d]: %5d <-\n" ,k, quant_info.Tb[k]);
    }

    fclose(fp);
} 
