#include "common.h"
#include <vector>

#ifndef _MONITOR_H_
#define _MONITOR_H_

using namespace std;

void savebayertxt(const vector<cfapix> &cfaimg, const string& filename, const BayerInfo& bayer_info);
void dumpmsst(const vector<msst> &msstv, const string& filename, const MsstInfo& msst_info);
void dumpmsst(const vector<msstSm> &msstv, const string& filename, const MsstInfo& msst_info);
void write_bitstream(const vector<vector<uint8_t>>& entp, const string& filename);
void read_bitstream ( vector<vector<uint8_t>>& entp, const string& filename);
#endif //_MONITOR_H_