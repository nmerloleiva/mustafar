// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <chrono>
#include <stdio.h>
#include <tchar.h>

#define CREATE_DB 0

#define DOUBLE_PRECISION 0
#if DOUBLE_PRECISION
typedef double REAL;
#else
typedef float REAL;
#endif

#define LAMBDA_DEN 9e+4 // 9e+5 gives epsilon out of range 

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
