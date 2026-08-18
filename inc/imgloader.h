
#include <stdio.h>
#include <iostream>
#include <string>
#include "common.h"

#ifndef _IMGLOADER_H_
#define _IMGLOADER_H_

using namespace std;

// Reads a Bayer RAW image encoded either as a continuous compact bitstream or
// with the CSI-2 RAW10/12/14 byte layout.  expected_pixels, when non-zero,
// validates the expected payload size.  If the file is larger than that payload,
// excess_bytes_as_offset selects whether the extra bytes are skipped at the
// front (true) or discarded from the tail (false).
bool loadbayerimg(vector<cfapix>& cfaimg, const string& filename, int bpp,
                  bool csi2_style = false, size_t expected_pixels = 0,
                  bool excess_bytes_as_offset = true);

#endif //_IMGLOADER_H_
