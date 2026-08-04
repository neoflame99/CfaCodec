#include "loadparam.h"

#include <fstream>
#include <sstream>


static std::string trim(const std::string& value)
{
    const std::string whitespace = " \t\r\n";
    const size_t begin = value.find_first_not_of(whitespace);
    if(begin == std::string::npos){
        return "";
    }

    const size_t end = value.find_last_not_of(whitespace);
    return value.substr(begin, end - begin + 1);
}

static bool parseInt(const std::string& value, int32_t& out)
{
    std::istringstream iss(value);
    int32_t parsed = 0;
    iss >> parsed;
    if(iss.fail()){
        return false;
    }

    iss >> std::ws;
    if(!iss.eof()){
        return false;
    }

    out = parsed;
    return true;
}
static bool parseUInt(const std::string& value, uint32_t& out)
{
    std::istringstream iss(value);
    uint32_t parsed = 0;
    iss >> parsed;
    if(iss.fail()){
        return false;
    }

    iss >> std::ws;
    if(!iss.eof()){
        return false;
    }

    out = parsed;
    return true;
}

bool loadBayerImgInfo(BayerImgInfo& bayer_img_info, SaveInfo& save_info, ProcessInfo& proc_info, const std::string& param_file)
{
    std::ifstream ifs(param_file);
    if(!ifs.is_open()){
        return false;
    }

    BayerImgInfo parsed_info;
    bool has_filename = false;
    bool has_width = false;
    bool has_height = false;
    bool has_bpp = false;

    std::string line;
    while(std::getline(ifs, line)){
        const size_t comment_pos = line.find("##");
        if(comment_pos != std::string::npos){
            line = line.substr(0, comment_pos);
        }

        line = trim(line);
        if(line.empty()){
            continue;
        }

        const size_t equal_pos = line.find('=');
        if(equal_pos == std::string::npos){
            continue;
        }

        const std::string key = trim(line.substr(0, equal_pos));
        const std::string value = trim(line.substr(equal_pos + 1));

        save_info.bitstream_filename = "bitstream.bin";
        save_info.save_cfaimg = false;
        save_info.save_cfaproc = false;
        save_info.save_msstv_enc = false;
        save_info.save_msstv_dec = false;
        if(key == "BayerImg"){
            parsed_info.filename = value;
            has_filename = !value.empty();
        }else if(key == "width"){
            has_width = parseInt(value, parsed_info.w);
        }else if(key == "height"){
            has_height = parseInt(value, parsed_info.h);
        }else if(key == "bpp"){
            has_bpp = parseInt(value, parsed_info.bpp);
        }else if(key == "bitstream_filename"){
            save_info.bitstream_filename = value;
        }else if(key == "save_cfaimg"){
            save_info.save_cfaimg = (value == "true" || value == "1");
        }else if(key == "save_cfaproc"){
            save_info.save_cfaproc = (value == "true" || value == "1");
        }else if(key == "save_msstv_enc"){
            save_info.save_msstv_enc = (value == "true" || value == "1");
        }else if(key == "save_msstv_dec"){
            save_info.save_msstv_dec = (value == "true" || value == "1");
        }else if(key == "g"){
            float g_value;
            std::istringstream iss(value);
            iss >> g_value;
            if(!iss.fail() && iss.eof()){
                proc_info.g = g_value;
            }
        }else if(key == "dwt_lv"){
            uint32_t dwt_lv;
            if(parseUInt(value, dwt_lv)){
                proc_info.dwt_lv = dwt_lv;
            }
        }else if(key == "ngrp"){
            uint32_t ngrp;
            if(parseUInt(value, ngrp)){
                proc_info.ngrp = ngrp;
            }
        //}else if(key == "bcw"){
        //    uint32_t bcw;
        //    if(parseUInt(value, bcw)){
        //        proc_info.bcw = bcw;
        //    }
        //}else if(key == "mbc"){
        //    uint32_t mbc;
        //    if(parseUInt(value, mbc)){
        //        proc_info.mbc = mbc;
        //    }
        }else if(key == "sgpcd"){
            uint32_t sgpcd;
            if(parseUInt(value, sgpcd)){
                proc_info.sgpcd = sgpcd;
            }
        }else if(key == "sel_ycc"){
            proc_info.sel_ycc = (value == "true" || value == "1");
        }
    }

    if(!(has_filename && has_width && has_height && has_bpp)){
        return false;
    }

    bayer_img_info = parsed_info;
    return true;
}
