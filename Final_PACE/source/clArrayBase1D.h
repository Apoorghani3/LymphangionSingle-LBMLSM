// clArrayBase1D.h: interface for the clArrayBase1D class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAYBASE1D_H__INCLUDED_)
#define AFX_ARRAYBASE1D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <stdio.h>
#include <memory.h>

typedef int BOOL;
#ifndef TRUE
#define FALSE 0
#define TRUE 1
#endif

#ifdef _DEBUG
#define ARRAY_DEBUG
#endif //_DEBUG


#ifndef BIN_SAVE_DIR
#define BIN_SAVE_DIR	""
#endif //BIN_SAVE_DIR

template <class T> class clArrayBase1D 
{
protected:
	T *Array, *currEl;
	int Size1, FullSize, Step1;
	int flDelete;

	void ini(void)
	{
		Size1 = 0;
		FullSize = 0;
		Step1 = 0;
		flDelete = 0;
		Array = NULL;
	}

	//memory allocation
	void allocateMemory(void)
	{
		flDelete = 1;
		Array = allocateMemory(FullSize);
	}

	//allocate memory and return pointer
	T* allocateMemory(int size)
	{
		if(size <= 0)
			return NULL;

		T *p = new T [size];
		if(p == NULL)
		{
			printf("Error in memory allocation!");
			exit(0);
		}
		return p;
	}

	//memory deallocation
	void delArray(void)
	{
		if(flDelete == 1 && Array != NULL)
			delete [] Array;
		flDelete = 0;
		currEl = 0;
	}
	
public:
	clArrayBase1D() {}
	virtual ~clArrayBase1D() {}

	T* getArray(void){return Array;}
	int getSize1(void){return Size1;}
	int getStep1(void){return Step1;}
	int getFullSize(void){return FullSize;}
	//fill array by an integer value
	void setInitValue(const int val)
	{
		memset(Array, val, sizeof(T) * FullSize);
	}
	//fill array by a double value
	void fill(const T val)
	{
		T *p = Array;
		for(int i = 0; i < FullSize; i++) *(p++) = val;
	}
	//fill array by zeros
	void zeros()
	{
		setInitValue(0);
	}
	//divide all elements in array
	void div(const T val)
	{
		T *p = Array;
		for(int i = 0; i < FullSize; i++) *(p++) /= val;
	}
	//copy array between objects, array size must be the same
	void copy(clArrayBase1D &source)
	{
#ifdef ARRAY_DEBUG
		if(FullSize != source.getFullSize())
		{
			printf("Error in memory copy between arrays! Source size [%d], targed size [%d]", FullSize, source.getFullSize());
			exit(0);
		}
#endif //ARRAY_DEBUG

		memcpy(Array, source.getArray(), sizeof(T) * FullSize);
	}
	//change pointers between arrays
	void change(clArrayBase1D &source)
	{
		T *p = Array;
		int fl = flDelete;
		Array = source.Array;
		flDelete = source.flDelete;
		source.Array = p;
		source.flDelete = fl;
	}
	//copy object without memory
	void copyPointer(clArrayBase1D &source)
	{
		delArray();
		Array = source.getArray();
		Size1 = source.getSize1();
		Step1 = source.getStep1();
		FullSize = source.getFullSize();
	}
	//return size of array
	virtual int getSize(int n)
	{
		if(n == 1)
			return Size1;
		return NULL;
	}
	//return memory increment step in array
	virtual int getStep(int n)
	{
		if(n == 1)
			return Step1;
		return NULL;
	}

	//save memory dump to file
	BOOL savebin(const char* FileName)
	{
		char Buf[80];
		sprintf(Buf, BIN_SAVE_DIR "%s.bin", FileName);
		return savebin_ext(Buf);
	};

	BOOL savebin_ext(const char* FileName)
	{
		// Open file in binary mode:
#ifndef CREATE_BIN_FILES
		return TRUE;
#endif //CREATE_BIN_FILES
		FILE *stream = fileopen(FileName, "w+b");
		if(stream == NULL)
			return FALSE;

		BOOL ret = savebin(stream);
		fclose(stream);
		if(ret != TRUE)
			printf("\nProblem with output to file: %s\n", FileName);
		return ret;
	};

	BOOL savebin(FILE *stream)
	{
		if(stream != NULL)
		{
#ifdef CREATE_BIN_FILES
			fwrite(Array, sizeof( T ), FullSize, stream);
#endif //CREATE_BIN_FILES
			return TRUE;
		}
		return FALSE;
	};

	//load memory dump from a file
	BOOL loadbin(const char* FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.bin", FileName);
		return loadbin_ext(Buf);
	};
	BOOL loadbin_ext(const char* FileName)
	{
		FILE *stream = fopen(FileName, "r+b" );
		if(stream == NULL)
			return FALSE;

		BOOL ret = loadbin(stream);
		fclose(stream);
		return ret;
	};
	BOOL loadbin(FILE *stream)
	{
		if(stream != NULL)
		{
			fread(Array, sizeof( T ), FullSize, stream);
			return TRUE;
		}
		return FALSE;
	};
	int count(const char * FileName)
	{
		int counter = countbin(FileName);
		if(counter > 0) return counter;
		return counttxt(FileName);
	}
	int counttxt(const char * FileName)
	{
		int counter = 0;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return 0;
	
		while(fscanf(stream, "%lg", &Buf) > 0)
			counter++;

		fclose(stream);
		return counter;
	};
	int countbin(const char * FileName)
	{
		int counter = 0;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.bin", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r+b")) == NULL)
			return 0;
	
		while(fread(&Buf, sizeof(T), 1, stream) > 0)
			counter++;

		fclose(stream);
		return counter;
	};
	void delay(double sec)
	{
		time_t StartTime, CurrTime;
		time(&StartTime);
		do
		{
			time(&CurrTime);
		}
		while(difftime(CurrTime, StartTime) < sec);
	};

	FILE* fileopen(const char* FileName, const char* Mode)
	{
		return fileopen(FileName, Mode, 100);
	};

	FILE* fileopen(const char* FileName, const char* Mode, const int att)
	{
		FILE *stream;
		for(int i = 0; i < att; i++)
		{
			stream = fopen(FileName, Mode);
			if(stream != NULL)
				return stream;
			delay(0.1);
		}
		return NULL;
	};


};

#endif // !defined(AFX_ARRAYBASE1D_H__INCLUDED_)
