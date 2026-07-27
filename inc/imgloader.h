
#include <stdio.h>
#include <iostream>
#include <string>
#include "common.h"

#ifndef _IMGLOADER_H_
#define _IMGLOADER_H_

using namespace std;

void loadbayerimg(vector<cfapix>& cfaimg, const string filename, const int bpp = 2);

#endif //_IMGLOADER_H_
