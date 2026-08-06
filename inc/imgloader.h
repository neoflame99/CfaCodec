
#include <stdio.h>
#include <iostream>
#include <string>
#include "common.h"

#ifndef _IMGLOADER_H_
#define _IMGLOADER_H_

using namespace std;

// Reads a Bayer RAW image encoded either as a continuous compact bitstream or
// with the CSI-2 RAW10/12/14 byte layout.  expected_pixels, when non-zero,
// also validates that the input has exactly the expected size.
bool loadbayerimg(vector<cfapix>& cfaimg, const string& filename, int bpp,
                  bool csi2_style = false, size_t expected_pixels = 0);

#endif //_IMGLOADER_H_
