// clOutput.h: interface for the clOutput class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLOUTPUT_H__INCLUDED_)
#define AFX_CLOUTPUT_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"

class clOutput  
{
	FILE *hCurrOut;
	char *outDir;
	int saveStep, savePeriod;
	int curPeriod;
	double FieldOutTime, periodFieldSaves, fieldOutPeriodStep;
	clArray3Dd LBout, LSout;
	clArray1Dd Timeout;
	int curLine;

	int avrCounter, avrPeriod;

	void RenameOutputFiles();
	
	void MakeDir(int Period)
	{
		char NewDir[80];
		sprintf(NewDir, "%s" SLASH "%d", outDir, Period);
#if defined(_WIN32)
		_mkdir(NewDir);
#else 
		mkdir(NewDir, S_IRWXU | S_IRWXG | S_IRWXO);
#endif //_WIN32
	}

	
	void RenameFile(char* SourceName)
	{
		RenameFile(SourceName, ".txt");
	}
	void RenameFileExt(char* SourceName)
	{
		RenameFile(SourceName, "");
	}
	void RenameFile(char* SourceName, char* SourceExt)
	{
		char DestFileName[40];
		char SourceFileName[40];

		sprintf(DestFileName, "%s" SLASH "%d" SLASH "%s%s", 
			outDir, countRename, SourceName, SourceExt);
		sprintf(SourceFileName, "%s%s", SourceName, SourceExt);
		
		RenameFileAtt(SourceFileName, DestFileName);
	}
	void RenameFileAtt(char* Source, char* Dest)
	{
		int maxAttempt = 100;

		for(int i = 0; i < maxAttempt; i++)
		{
			remove(Dest);
			if(rename(Source, Dest) == NULL)
				return;
			delay(0.1);
			FILE *s = fopen(Source, "r");
			if(s != NULL)
				fclose(s);
			else
				return;
		}
	}

	void delay(double sec)
	{
		time_t StartTime, CurrTime;
		time(&StartTime);
		do
		{
			time(&CurrTime);
		}
		while(difftime(CurrTime, StartTime) < sec);
	}

	void RenameFile(char* SourceName, char* SourceExt, int counter)
	{
		char DestFileName[40];
		char SourceFileName[40];

		sprintf(DestFileName, "%s" SLASH "%d" SLASH "%s%d.%s", 
					outDir, curPeriod, SourceName, counter, SourceExt);
		sprintf(SourceFileName, "%s.%s", SourceName, SourceExt);

		RenameFileAtt(SourceFileName, DestFileName);
	}

	FILE * CopyCurrState(FILE *hFileTmp, char* FileName, char* FileNameTmp)
	{
		FILE *hFileCopy;
		const int BufLength = 1024;
		int ReadLength;
		char Buffer[BufLength];

		if(hFileTmp != NULL)
			fclose(hFileTmp);

		int maxAttempt = 100;
		for(int i = 0; i < maxAttempt; i++)
		{
			hFileCopy = fopen(FileName, "a");
			if(hFileCopy != NULL)
				break;

			delay(0.1);
		}

		if(hFileCopy == NULL)
			return NULL;

		hFileTmp = fopen(FileNameTmp, "r");

		if(hFileTmp == NULL)
			return NULL;

		do
		{
			ReadLength = (int)fread(Buffer, sizeof(char), BufLength, hFileTmp);
			fwrite(Buffer, sizeof(char), ReadLength, hFileCopy);
		}
		while(ReadLength == BufLength);

		fclose(hFileCopy);
		fclose(hFileTmp);
		hFileTmp = fopen(FileNameTmp, "w+");
		return hFileTmp;
	}

	void CopyCurrState(char* FileName)
	{
		FILE *hFile;

		int maxAttempt = 100;
		for(int i = 0; i < maxAttempt; i++)
		{
			hFile = fopen(FileName, "a");
			if(hFile != NULL)
				break;
			printf("cannot open output file %s, waiting", FileName);
			delay(0.1);
		}

		if(hFile == NULL)
		{
			printf("cannot open output file %s, data lost", FileName);
			curLine = 0;
			return;
		}

		int nLBfliuds = LBout.getSize2();
		int nLBfileds = LBout.getSize3();
		int nLSobjects = LSout.getSize2();
		int nLSfields = LSout.getSize3();
		for(int i = 0; i < curLine; i++)
		{
			fprintf(hFile, "%g", Timeout(i));
			for(int j = 1; j < nLBfliuds; j++)
				for(int k = 0; k < nLBfileds; k++)
					fprintf(hFile, "\t%g", LBout(i,j,k));
			for(int j = 0; j < nLSobjects; j++)
				for(int k = 0; k < nLSfields; k++)
					fprintf(hFile, "\t%g", LSout(i,j,k));
			fprintf(hFile, "\n");
		}

		fclose(hFile);
		curLine = 0;
		return;
	}

	BOOL CopyFile(char* NameSrc, char* NameDest)
	{
		FILE *hSrc, *hDest;
		const int BufLength = 1024;
		int ReadLength;
		char Buffer[BufLength];
		int maxAttempt = 60;

		for(int i = 0; i < maxAttempt; i++)
		{
			hSrc = fopen(NameSrc, "r");
			if(hSrc != NULL)
				break;

			delay(0.1);
		}
		if(hSrc == NULL)
			return FALSE;

		for(int i = 0; i < maxAttempt; i++)
		{
			hDest = fopen(NameDest, "w+");
			if(hDest != NULL)
				break;

			delay(0.1);
		}
		if(hSrc == NULL)
		{
			fclose(hSrc);
			return FALSE;
		}

		do
		{
			ReadLength = (int)fread(Buffer, sizeof(char), BufLength, hSrc);
			fwrite(Buffer, sizeof(char), ReadLength, hDest);
		}
		while(ReadLength == BufLength);

		fclose(hSrc);
		fclose(hDest);
		return TRUE;
	}
	BOOL CleanFile(char* Name)
	{
		int maxAttempt = 100;
		for(int i = 0; i < maxAttempt; i++)
		{
			FILE *sOut = fopen(CURR_OUTPUT_FILE, "w+");
			if(sOut != NULL)
			{
				fclose(sOut);
				return TRUE;
			}
			delay(0.1);
		}
		return FALSE;
	}
	//template <class T> void average(T &source, clArray2Dd &dest)
	//{
	//	int i, j;
	//	for(i = 0; i < c.nX; i++)
	//		for(j = 0; j < c.nY; j++)
	//			dest(i,j) += source(i,j);
	//}

	void zeroAverage()
	{
	}

	void saveAverage()
	{
		MakeDir(avrPeriod);
	}

	void saveAverage(clArray2Dd &avr, char* FileName)
	{
		char DestFileName[80];
		sprintf(DestFileName, "%s" SLASH "%d" SLASH "%s", outDir, avrPeriod, FileName);

		avr.div(avrCounter);
		avr.save2file(DestFileName);
	}

	template <class T> void saveBin(T &bin, char* FileName)
	{
		char DestFileName[80];
		sprintf(DestFileName, "%s" SLASH "%s", outDir, FileName);
		bin.savebin(FileName);
	}

public:
	int countRename;
	clOutput()
	{
		hCurrOut = NULL;
		//hCurrOut = fopen(TMP_OUTPUT_FILE, "w+");
	};
	virtual ~clOutput()
	{
		if(hCurrOut != NULL) fclose(hCurrOut);
	};
	void ini();
	void start();

	void step()
	{
		current();
	};

	void finish()
	{
//		current();
		field();
	};

	void current();
	void field();
	void parameters();
	void average();
	void addAverage();
	void bin();
	void saveBin();

};

#endif // !defined(AFX_CLOUTPUT_H__INCLUDED_)
