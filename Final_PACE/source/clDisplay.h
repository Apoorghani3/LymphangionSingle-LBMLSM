// clDisplay.h: interface for the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
#define AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include "clConstants.h"

class clDisplay  
{
	time_t StartTimeSys;
	char cStartTimeSys[26];
	int counter;

public:
	clDisplay(){counter = 0;};
	virtual ~clDisplay(){};
	void error(char* err)
	{
//		cout << err << endl;
	};
	void start(){ShowTime();}
	void finish()
	{
		ShowTime();
#if defined(_WIN32)
		printf("\n");
#endif //_WIN32
	}
	void step()
	{
		if(++counter < c.t.DisplaySignal) 
			return;

		ShowTime();
		counter = 0;
	}

	void field(void);
	void ShowTime(void);
	void ShowdP(void);

	void ini(void);

};

#endif // !defined(AFX_CLDISPLAY_H__904010C9_8CE3_405E_9266_71D4D0D77458__INCLUDED_)
