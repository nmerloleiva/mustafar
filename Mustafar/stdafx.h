// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

#include "targetver.h"

#include <chrono>
#include <stdio.h>
#include <tchar.h>
# include <conio.h>
# include <malloc.h>

typedef double REAL;
#define REAL_EPSILON DBL_EPSILON

#define experimentalDisturbance 0.001 //cambiarla para varios calculos

#define LAMBDA_DEN 9e+4 // 9e+5 gives epsilon out of range 

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
