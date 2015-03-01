// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once

// This is to stop visual studio from giving a heap load of warnings because it does not like fopen and such calls.
// Just commet the next 3 lines out, if you want to see.
#ifdef _WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif


//#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include <string>

// TODO: reference additional headers your program requires here
#include "constants.h"
#include "FastEnsemble.h"
#include "FastIon.h"
#include "forces.h"
#include "integrator.h"
#include "BisectionRootFinder.h"
#include "Trap.h"