// ArrayBase3D.h: interface for the ArrayBase3D class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAYBASE3D_H__INCLUDED_)
#define AFX_ARRAYBASE3D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "clArrayBase2D.h"
template <class T> class clArrayBase2D;

template <class T> class clArrayBase3D: public clArrayBase2D <T>
{
protected:
	int Size3, Step3;

public:
	clArrayBase3D(){}
	virtual ~clArrayBase3D(){}

	inline int getSize3(void){return Size3;}
	inline int getStep3(void){return Step3;}
	
	inline virtual int getSize(int n)
	{
		if(n == 3)
			return Size3;
		return clArrayBase2D <T>::getSize(n);
	}

	inline virtual int getStep(int n)
	{
		if(n == 3)
			return Step3;
		return clArrayBase2D <T>::getStep(n);
	}
};

#endif // !defined(AFX_ARRAYBASE3D_H__INCLUDED_)
