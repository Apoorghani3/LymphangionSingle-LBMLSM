// clProblem.cpp: implementation of the clProblem class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clConstants.h"
//#include "clGrid.h"
#include "clMethod.h"
#include "clVariables.h"
#include "clProblem.h"


clDisplay clProblem::d;
clOutput clProblem::o;

void clProblem::ini(void)
{
	c.ini();
//	g.ini();
	v.ini();
	vls.ini();
	vlb.ini();

#ifdef TR_SOLVER
	vtr.ini();
#endif

	m.ini();
	o.ini();
	d.ini();
}

void clProblem::update(void)
{
}

void clProblem::run(void)
{
	m.run();
}

void clProblem::start(void)
{
	d.start();
	o.start();
}

void clProblem::step(void)
{
	d.step();
//	o.step();
}

void clProblem::debug(void)
{
	d.step();
	o.field();
}

void clProblem::finish(void)
{
	d.finish();
	o.finish();
}

void clProblem::error(void)
{
	error("Program error!");
}

void clProblem::error(char *err)
{
	d.error(err);

	finish();
	exit(0);
}

