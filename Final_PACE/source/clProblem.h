// clProblem.h: interface for the clProblem class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLPROBLEM_H__INCLUDED_)
#define AFX_CLPROBLEM_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clOutput.h"
#include "clDisplay.h"

class clProblem 
{
public:
	clProblem(){};
	virtual ~clProblem(){};
	static void ini(void);
	static void update(void);
	static void	run(void);
	static void debug(void);
	static void start(void);
	static void step(void);
	static void finish(void);
	static void error();
	static void error(char*);
	static clOutput o;
	static clDisplay d;
};

#endif // !defined(AFX_CLPROBLEM_H__INCLUDED_)
