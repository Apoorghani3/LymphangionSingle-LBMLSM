// clVariables.cpp: implementation of the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clVariables.h"

void clVariables::ini(void)
{
}

void clVariables::update(void)
{
}

void clVariables::error(char* err)
{
	p.error(err);
}

void clVariablesLS::clD::remove_link(int n1, int n2)
{
	vlb.M.remove_s(n1, n2);

	c.ls.removeBL(n1, n2);
//	c.ls.removedoubleBL();

	c.ls.N.remove_links(n1, n2);

//	c.ls.N.test_links();

//	vlb.M.inifill_S();

};
void clVariablesLS::clD::correct_tables()
{
#ifdef LS_RUPTURE
	c.ls.removeextraBL();
	c.ls.N.test_links();
	vlb.M.inifill_S();
	c.ls.updateO();
	vls.P.force_update();
#endif //LS_RUPTURE
};


