// clDisplay.cpp: implementation of the clDisplay class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clDisplay.h"
#include "clVariables.h"
#include "clOutput.h"

void clDisplay::ShowTime()
{
	//#ifndef _DEBUG
	time_t CurrTimeSys, ElapsTimeSys, FinishTimeSys, RemTimeSys;
	char cFinishTimeSys[26];

	time(&CurrTimeSys);

	ElapsTimeSys = CurrTimeSys - StartTimeSys;
	FinishTimeSys = StartTimeSys + (time_t)((double)ElapsTimeSys * c.t.StopTime / c.t.Time);
	int ElapsH = 0, ElapsM = 0, ElapsS = 0;
	if(ElapsTimeSys > 0)
	{
		ElapsH = (int)(ElapsTimeSys / (60 * 60));
		ElapsM = (int)((ElapsTimeSys - ElapsH * 60 * 60) / 60);
		ElapsS = (int)((ElapsTimeSys - ElapsH * 60 * 60 - ElapsM * 60));
	}
	int RemH = -1, RemM = 0, RemS = 0;
	if(FinishTimeSys > 0)
	{
		RemTimeSys = FinishTimeSys - CurrTimeSys;
		RemH = (int)(RemTimeSys / (60 * 60));
		RemM = (int)((RemTimeSys - RemH * 60 * 60) / 60);
		RemS = (int)((RemTimeSys - RemH * 60 * 60 - RemM * 60));
	}

#if defined(_WIN32)
	if(FinishTimeSys > 0)
		strcpy(cFinishTimeSys,&(ctime(&FinishTimeSys)[11]));
	else
		sprintf(cFinishTimeSys,"00:00:00");

#ifdef TR_SOLVER
	//		printf("M%9.2f R%6.3f E%6.3f S%.8s F%.8s E%.8s R%.8s P%6.2f%% \r", 
	//printf("M%8.1f R%5.3f Dep%5.4f S%5.3f N%d/%d BN%d O%d/%d E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%% \r", 
	//	vlb.Ro.mass(1), vlb.Ro.avr(1),  vtr.A.avr, vls.D.Smax,
	//	c.ls.nN, c.ls.nNd, c.ls.nBN, c.ls.nO, c.ls.nOp - c.ls.nNd,
	//	ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 
	//	100. * c.t.Time / c.t.StopTime);
	printf("M%8.1f R%5.3f Dep%5.4f S%5.3f N%d/%d BN%d E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%% \r", 
		vlb.Ro.mass(1), vlb.Ro.avr(1),  vtr.C.OutputVal.x, vls.D.Smax,
		c.ls.nN, c.ls.nNd, c.ls.nBN,
		ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 
		100. * c.t.Time / c.t.StopTime);


#else
#ifdef BEADS_MIXING
	//		printf("M%9.2f R%6.3f E%6.3f S%.8s F%.8s E%.8s R%.8s P%6.2f%% \r", 
	printf("M%8.1f R%5.3f Mix%5.4f S%5.3f N%d/%d BN%d O%d/%d E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%% \r", 
		vlb.Ro.mass(1), vlb.Ro.avr(1), vlb.NM.S, vls.D.Smax,
		c.ls.nN, c.ls.nNd, c.ls.nBN, c.ls.nO, c.ls.nOp - c.ls.nNd,
		ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 
		100. * c.t.Time / c.t.StopTime);
#else
		printf("M%8.1f R%5.3f S%5.3f N%d/%d BN%d O%d/%d E%d:%.2d:%.2d R%d:%.2d:%.2d P%5.2f%% \r", 
		vlb.Ro.mass(1), vlb.Ro.avr(1), vls.D.Smax,
		c.ls.nN, c.ls.nNd, c.ls.nBN, c.ls.nO, c.ls.nOp - c.ls.nNd,
		ElapsH, ElapsM, ElapsS, RemH, RemM, RemS, 
		100. * c.t.Time / c.t.StopTime);
#endif // BEADS_MIXING
#endif // TR_SOLVER

#else 

	if(FinishTimeSys > 0)
		strcpy(cFinishTimeSys,ctime(&FinishTimeSys));
	else
		sprintf(cFinishTimeSys,"-1:00:00\n");

	char cTimeSys[26];
	strcpy(cTimeSys,ctime(&CurrTimeSys));

	FILE *stream;
	const char *FileName = UNIX_OUTPUT_FILE;
	
#ifdef TR_SOLVER
	if((stream = fopen(FileName, "a")) != NULL)
	{
		//fprintf(stream, "Mass %9.2f Density %6.3f Strain %6.3f Dep %4.3f Percent %6.2f%%\n Time %.26s Start %.26s Finish %.26s Elapsed %d:%.2d:%.2d\n Remain %d:%.2d:%.2d\n\n",  
		//	vlb.Ro.mass(1), vlb.Ro.avr(1), vls.D.Smax, vtr.A.avr, 100. * c.t.Time / c.t.StopTime,
		//	cTimeSys, cStartTimeSys, cFinishTimeSys, ElapsH, ElapsM, ElapsS, RemH, RemM, RemS);
		//fclose(stream);
		fprintf(stream, "Mass %9.2f Density %6.3f Strain %6.3f Dep %4.3f Percent %6.2f%%\n Time %.26s Start %.26s Finish %.26s Elapsed %d:%.2d:%.2d\n Remain %d:%.2d:%.2d\n\n",  
			vlb.Ro.mass(1), vlb.Ro.avr(1), vls.D.Smax, vtr.C.OutputVal.x, 100. * c.t.Time / c.t.StopTime,
			cTimeSys, cStartTimeSys, cFinishTimeSys, ElapsH, ElapsM, ElapsS, RemH, RemM, RemS);
		fclose(stream);
	}
#else
#ifdef BEADS_MIXING
		if((stream = fopen(FileName, "a")) != NULL)
	{
		fprintf(stream, "Mass %9.2f Density %6.3f Mix %4.3f Strain %6.3f Percent %6.2f%%\n Time %.26s Start %.26s Finish %.26s Elapsed %d:%.2d:%.2d\n Remain %d:%.2d:%.2d\n\n",  
			vlb.Ro.mass(1), vlb.Ro.avr(1), vlb.NM.S, vls.D.Smax, 100. * c.t.Time / c.t.StopTime,
			cTimeSys, cStartTimeSys, cFinishTimeSys, ElapsH, ElapsM, ElapsS, RemH, RemM, RemS);
		fclose(stream);
	}
#else
		if((stream = fopen(FileName, "a")) != NULL)
	{
		fprintf(stream, "Mass %9.2f Density %6.3f Strain %6.3f Percent %6.2f%%\n Time %.26s Start %.26s Finish %.26s Elapsed %d:%.2d:%.2d\n Remain %d:%.2d:%.2d\n\n",  
			vlb.Ro.mass(1), vlb.Ro.avr(1), vls.D.Smax, 100. * c.t.Time / c.t.StopTime,
			cTimeSys, cStartTimeSys, cFinishTimeSys, ElapsH, ElapsM, ElapsS, RemH, RemM, RemS);
		fclose(stream);
	}
#endif // BEADS_MIXING
#endif // TR_SOLVER
#endif //_WIN32
	//#endif //_DEBUG
};
void clDisplay::field()
{
#if !defined(_WIN32)
	time_t CurrTimeSys;
	char cTimeSys[26];

	time(&CurrTimeSys);
	strcpy(cTimeSys,ctime(&CurrTimeSys));

	FILE *stream;
	const char *FileName = UNIX_OUTPUT_FILE;
	if((stream = fopen(FileName, "a")) != NULL)
	{
		fprintf(stream, "Snapshot done #%d : Time %.26s\n", p.o.countRename, cTimeSys);
		fclose(stream);
	}
#endif //_WIN32
}

void clDisplay::ini()
{
	counter = 0;
	time(&StartTimeSys);

#if !defined(_WIN32)
	strcpy(cStartTimeSys,ctime(&StartTimeSys));

	const char *FileName = "out.txt";
	FILE *stream;
	stream = fopen(FileName, "w+");
	if(stream == NULL)
		return;

	fclose(stream);
#else //_WIN32
	strcpy(cStartTimeSys,&(ctime(&StartTimeSys)[11]));
#endif //_WIN32
};
