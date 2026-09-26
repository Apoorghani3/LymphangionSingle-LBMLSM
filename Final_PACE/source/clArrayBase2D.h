// ArrayBase2D.h: interface for the ArrayBase2D class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAYBASE2D_H__INCLUDED_)
#define AFX_ARRAYBASE2D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "clArrayBase1D.h"
template <class T> class clArrayBase1D;

template <class T> class clArrayBase2D: public clArrayBase1D <T> 
{
protected:
	int Size2, Step2;

public:
	clArrayBase2D(){}
	virtual ~clArrayBase2D(){}

	inline int getSize2(void){return Size2;}
	inline int getStep2(void){return Step2;}
	
	inline virtual int getSize(int n)
	{
		if(n == 2)
			return Size2;
		return clArrayBase1D<T>::getSize(n);
	}

	inline virtual int getStep(int n)
	{
		if(n == 2)
			return Step2;
		return clArrayBase1D<T>::getStep(n);
	}

};

#endif // !defined(AFX_ARRAYBASE2D_H__INCLUDED_)
