#include "imgloader.h"

namespace {

bool is_supported_bpp(int bpp)
{
    return bpp == 10 || bpp == 12 || bpp == 14 || bpp == 16;
}

bool get_csi2_group_size(int bpp, size_t& pixels_per_group, size_t& bytes_per_group)
{
    switch(bpp){
        case 10: pixels_per_group = 4; bytes_per_group = 5; return true;
        case 12: pixels_per_group = 2; bytes_per_group = 3; return true;
        case 14: pixels_per_group = 4; bytes_per_group = 7; return true;
        case 16: pixels_per_group = 1; bytes_per_group = 2; return true;
        default: return false;
    }
}

bool get_expected_size(size_t pixels, int bpp, bool csi2_style, size_t& bytes)
{
    if(csi2_style){
        size_t pixels_per_group = 0;
        size_t bytes_per_group = 0;
        if(!get_csi2_group_size(bpp, pixels_per_group, bytes_per_group) ||
           pixels % pixels_per_group != 0){
            return false;
        }
        bytes = (pixels / pixels_per_group) * bytes_per_group;
        return true;
    }

    bytes = (pixels * static_cast<size_t>(bpp) + 7) / 8;
    return true;
}

bool trim_to_expected_payload(vector<uint8_t>& input, size_t expected_bytes,
                              bool excess_bytes_as_offset, const char* packing_name)
{
    if(input.size() < expected_bytes){
        fprintf(stderr, "loadbayerimg: input size mismatch for %s: expected at least %zu bytes, got %zu bytes\n",
                packing_name, expected_bytes, input.size());
        return false;
    }

    if(input.size() == expected_bytes){
        return true;
    }

    const size_t extra_bytes = input.size() - expected_bytes;
    const size_t src_offset = excess_bytes_as_offset ? extra_bytes : 0;
    vector<uint8_t> payload(input.begin() + static_cast<std::ptrdiff_t>(src_offset),
                            input.begin() + static_cast<std::ptrdiff_t>(src_offset + expected_bytes));
    input.swap(payload);
    fprintf(stderr, "loadbayerimg: input has %zu extra bytes; %s\n",
            extra_bytes,
            excess_bytes_as_offset ? "using them as a front offset" : "discarding them from the tail");
    return true;
}

bool unpack_compact(vector<cfapix>& pixels, const vector<uint8_t>& input, int bpp, size_t count)
{
    const uint32_t mask = (1U << bpp) - 1U;
    pixels.reserve(count);
    size_t bit_pos = 0;
    for(size_t i = 0; i < count; ++i){
        uint32_t value = 0;
        for(int bit = 0; bit < bpp; ++bit, ++bit_pos){
            value |= static_cast<uint32_t>((input[bit_pos / 8] >> (bit_pos % 8)) & 0x1) << bit;
        }
        pixels.push_back(static_cast<cfapix>(value & mask));
    }
    return true;
}

bool unpack_csi2(vector<cfapix>& pixels, const vector<uint8_t>& input, int bpp)
{
    size_t pixels_per_group = 0;
    size_t bytes_per_group = 0;
    if(!get_csi2_group_size(bpp, pixels_per_group, bytes_per_group) ||
       input.size() % bytes_per_group != 0){
        return false;
    }

    pixels.reserve((input.size() / bytes_per_group) * pixels_per_group);
    for(size_t offset = 0; offset < input.size(); offset += bytes_per_group){
        switch(bpp){
            case 10: {
                const uint8_t low = input[offset + 4];
                for(size_t i = 0; i < 4; ++i){
                    pixels.push_back((static_cast<cfapix>(input[offset + i]) << 2) |
                                     ((low >> (i * 2)) & 0x3));
                }
                break;
            }
            case 12: {
                const uint8_t low = input[offset + 2];
                pixels.push_back((static_cast<cfapix>(input[offset]) << 4) | (low & 0xF));
                pixels.push_back((static_cast<cfapix>(input[offset + 1]) << 4) | ((low >> 4) & 0xF));
                break;
            }
            case 14: {
                const uint32_t low = static_cast<uint32_t>(input[offset + 4]) |
                                     (static_cast<uint32_t>(input[offset + 5]) << 8) |
                                     (static_cast<uint32_t>(input[offset + 6]) << 16);
                for(size_t i = 0; i < 4; ++i){
                    pixels.push_back((static_cast<cfapix>(input[offset + i]) << 6) |
                                     ((low >> (i * 6)) & 0x3F));
                }
                break;
            }
            case 16:
                pixels.push_back(static_cast<cfapix>(input[offset]) |
                                 (static_cast<cfapix>(input[offset + 1]) << 8));
                break;
            default:
                return false;
        }
    }
    return true;
}

} // namespace

bool loadbayerimg(vector<cfapix>& cfaimg, const string& filename, int bpp,
                  bool csi2_style, size_t expected_pixels, bool excess_bytes_as_offset)
{
    cfaimg.clear();
    if(!is_supported_bpp(bpp)){
        fprintf(stderr, "loadbayerimg: unsupported bpp=%d (supported: 10, 12, 14, 16)\n", bpp);
        return false;
    }

    FILE* fp = fopen(filename.c_str(), "rb");
    if(fp == nullptr){
        fprintf(stderr, "loadbayerimg: file open error: %s\n", filename.c_str());
        return false;
    }
    if(fseek(fp, 0, SEEK_END) != 0){
        fclose(fp);
        return false;
    }
    const long file_length = ftell(fp);
    if(file_length < 0 || fseek(fp, 0, SEEK_SET) != 0){
        fclose(fp);
        return false;
    }

    vector<uint8_t> input(static_cast<size_t>(file_length));
    if(!input.empty() && fread(input.data(), 1, input.size(), fp) != input.size()){
        fclose(fp);
        return false;
    }
    fclose(fp);

    size_t pixel_count = expected_pixels;
    size_t expected_bytes = 0;
    if(expected_pixels != 0){
        if(!get_expected_size(expected_pixels, bpp, csi2_style, expected_bytes)){
            fprintf(stderr, "loadbayerimg: CSI-2 RAW%d needs a pixel count that is a multiple of its packing group\n", bpp);
            return false;
        }
        if(!trim_to_expected_payload(input, expected_bytes, excess_bytes_as_offset,
                                     csi2_style ? "CSI-2" : "compact")){
            return false;
        }
    }else if(csi2_style){
        // Without expected_pixels there is no known payload size, so there is
        // no reliable way to distinguish a leading offset from valid payload.
        size_t pixels_per_group = 0;
        size_t bytes_per_group = 0;
        get_csi2_group_size(bpp, pixels_per_group, bytes_per_group);
        if(input.size() % bytes_per_group != 0){
            fprintf(stderr, "loadbayerimg: CSI-2 RAW%d input has an incomplete pixel group\n", bpp);
            return false;
        }
        pixel_count = (input.size() / bytes_per_group) * pixels_per_group;
    }else{
        // Automatic compact-size inference treats all readable bytes as payload.
        pixel_count = (input.size() * 8) / static_cast<size_t>(bpp);
    }

    vector<cfapix> decoded;
    const bool decoded_ok = csi2_style
        ? unpack_csi2(decoded, input, bpp)
        : unpack_compact(decoded, input, bpp, pixel_count);
    if(!decoded_ok || decoded.size() != pixel_count){
        fprintf(stderr, "loadbayerimg: failed to decode %s Bayer RAW%d input\n",
                csi2_style ? "CSI-2" : "compact", bpp);
        return false;
    }

    cfaimg.swap(decoded);
    return true;
}
