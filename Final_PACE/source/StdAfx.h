// StdAfx.h : include file for standard system include files,
//  or project specific include files that are used frequently, but
//      are changed infrequently
// (c) Alexander Alexeev, 2006 
//

#ifndef _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES 1
#endif
//#define _CRT_SECURE_NO_DEPRECATE
#pragma warning(disable : 4996)

#if !defined(AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_)
#define AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if !defined(SLASH)
#if defined(_WIN32)
#define SLASH "\\"
#else 
#define SLASH "/"
#endif //_WIN32
#endif //SLASH

#if !defined(_WIN32)
#define _OPEN_MP
#endif //_WIN32

#if defined(_WIN32)
#define _OPERATION_SYSTEM  "WIN32"
#else 
#define _OPERATION_SYSTEM  "UNIX"
#endif //_WIN32


// TODO: reference additional headers your program requires here
#include "ConstDef.h"


#include <omp.h>
#ifdef LB_USE_OPENMP
#include <omp.h>
#endif //LB_USE_OPENMP

#include <stdlib.h>
#include <math.h>
//#include <iostream>
#if defined(_WIN32)
#include <direct.h>
#endif //_WIN32
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>

typedef int BOOL;
#ifndef TRUE
#define FALSE	0
#define TRUE	1
#endif

#define MOD(A, B)		(((A) % (B) + (B)) % (B))
#define MODFAST(A, B)	(((A) >= (B)) ? ((A) - (B)) : (((A) < 0) ? ((A) + (B)) : (A)))
#define MAX(A, B)		(((A) > (B)) ? (A) : (B))
#define MIN(A, B)		(((A) < (B)) ? (A) : (B))

#ifdef _DEBUG
//#define LB_DEBUG
//#define LS_DEBUG
//#define TR_DEBUG
#endif //_DEBUG

#define VERSION		6.00


#if !defined(__EXTERN_VARIABLES__)
#define __EXTERN_VARIABLES__

class clConstants;
class clVariables;
class clVariablesLS;
class clVariablesLB;
class clVariablesTr;
class clMethod;
class clGrid;
class clProblem;
class clOutput;
class clDisplay;

extern clConstants c;
extern clVariables v;
extern clVariablesLS vls;
extern clVariablesLB vlb;
extern clVariablesTr vtr;
extern clMethod m;
extern clGrid g;
extern clProblem p;
extern clOutput o;
extern clDisplay d;

#endif // !defined(__EXTERN_VARIABLES__)

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_STDAFX_H__24F045E1_A23A_42BE_BCBB_F3F2F6EE9AF1__INCLUDED_)

#include "clVector.h"
#include "clArray1D.h"
#include "clArray2D.h"
#include "clArray3D.h"

//#include "clConstants.h"
//#include "clMethod.h"
//#include "clVariables.h"
//#include "clProblem.h"
//#include "clOutput.h"
//#include "clDisplay.h"
