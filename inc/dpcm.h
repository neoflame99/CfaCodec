
#ifndef _PREDICT_H_
#define _PREDICT_H_

#include <cstdio>
#include <iostream>
#include "common.h"

void dpcm_mixc(cfa_fmt& ocfa, const cfa_fmt& icfa);
void dpcm_monoc(cfa_fmt& ocfa, const cfa_fmt& icfa);

void idpcm_mixc(cfa_fmt& ocfa, const cfa_fmt& icfa);
void idpcm_monoc(cfa_fmt& ocfa, const cfa_fmt& icfa);
#endif