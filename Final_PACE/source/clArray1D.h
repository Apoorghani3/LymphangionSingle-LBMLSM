// Array1D.h: interface for the Array1D class.
// 1D array derived from ArrayBase1D
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAY1D_H__INCLUDED_)
#define AFX_ARRAY1D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
//#include "StdAfx.h"
#include "clArrayBase1D.h"

template <class T> class  clArray2D;
template <class T> class  clArrayBase1D;

template <class T> class  clArray1D: public clArrayBase1D <T>
{
public:
	// Construction without memory allocation
	clArray1D(){clArrayBase1D<T>::ini();}
	// Construction with memory allocation
	clArray1D(const int n1)
	{
		ini();
		allocate(n1);
	}
	// Destruction
	virtual ~clArray1D(){clArrayBase1D<T>::delArray();}

	T& operator[](const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif 
		return *(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1);
	}

	T& operator()(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return *(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1);
	}

	T& el(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return *(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1);
	}

	//return element
	T getEl(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return *(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1);
	}
	//set element
	void setEl(const T val, const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		*(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1) = val;
	}
	//add to element
	void addEl(const T val, const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		*(clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1) += val;
	}
	//divide element by val
/*
	void div(const T val, const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		*(Array + n1 * Step1) /= val;
	}
*/
	//return pointer on element
	T * getPointerEl(const int n1)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1);
#endif
		return (clArrayBase1D<T>::Array + n1 * clArrayBase1D<T>::Step1);
	}
	//return refernece on element

	//memory allocation
	void allocate(const int n1)
	{
		clArrayBase1D<T>::delArray();

		this->Size1 = n1;
		this->Step1 = 1;
		this->FullSize = this->Size1;

		this->allocateMemory();
	}

	//memory reallocation
	void reallocate(const int n1)
	{
		T * p = clArrayBase1D<T>::allocateMemory(n1);
		memset(p, 0, sizeof(T) * n1);

		for(int i = 0; i < clArrayBase1D<T>::Size1 && i < n1; i++)
			*(p + 1 * i) = el(i);
		
		clArrayBase1D<T>::delArray();

		this->flDelete = 1;
		this->Array = p;
		this->Size1 = n1;
		this->Step1 = 1;
		this->FullSize = n1;
	}

	//allocates memory and fills by zeros
	void zeros(const int n1)
	{
		allocate(n1);
		zeros();
	}

	void zeros()
	{
		clArrayBase1D <T>::zeros();
	}

	void map1Don2D1(clArray2D<T> &p2D, const int n2)
	{
		setPointer(p2D.getPointerEl(0, n2), 
					p2D.getSize1(),  p2D.getStep1(), p2D.getFullSize());
	}

	void map1Don2D2(clArray2D<T> &p2D, const int n1)
	{
		setPointer(p2D.getPointerEl(n1, 0), 
					p2D.getSize2(), p2D.getStep2(), p2D.getFullSize());
	}

	//int findnonzero(int n)
	//{
	//	for(int i = n; i < clArrayBase1D<T>::Size1; i++)
	//		if(this->el(i) ~= 0)
	//			return i;
	//	return -1;
	//}

	clArray1D<T>& operator +=(const T val)
	{
		for(int i = 0; i < clArrayBase1D<T>::Size1; i++)
			el(i) += val;
		return *this;
	}

	clArray1D<T>& operator -=(const T val)
	{
		for(int i = 0; i < clArrayBase1D<T>::Size1; i++)
			el(i) -= val;
		return *this;
	}

	clArray1D<T>& operator *=(const T val)
	{
		for(int i = 0; i < clArrayBase1D<T>::Size1; i++)
			el(i) *= val;
		return *this;
	}

	clArray1D<T>& operator /=(const T val)
	{
		for(int i = 0; i < clArrayBase1D<T>::Size1; i++)
			el(i) /= val;
		return *this;
	}

	//save array to a file
	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && clArrayBase1D<T>::savebin(FileName);
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return save2file(Buf);
	}

	BOOL save2file(const char * FileName)
	{
		int i;
		FILE *stream;
		//if((stream = fopen(FileName, "w+")) == NULL)
		if((stream = this->fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fprintf(stream, "%g\n", (double)getEl(i)) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
		}

		fclose(stream);
		return TRUE;
	}

	//load array from a file
	BOOL loadtxt(const char * FileName)
	{
		int i;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;

		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fscanf(stream, "%lg", &Buf) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			setEl((T)Buf, i);
		}

		fclose(stream);
		return TRUE;
	}

	BOOL load(const char * FileName)
	{
		return !((clArrayBase1D<T>::loadbin(FileName) == FALSE) && (loadtxt(FileName) == FALSE));
	}
	
	BOOL loadtxtnew(const char * FileName)
	{
		zeros(this->counttxt(FileName));
		return loadtxt(FileName);
	}

	BOOL loadbinnew(const char * FileName)
	{
		zeros(countbin(FileName));
		return loadbin(FileName);
	}

	BOOL loadnew(const char * FileName)
	{
		int counter = this->countbin(FileName);
		if(counter > 0)
		{
			zeros(counter);
			return this->loadbin(FileName);
		}
		return loadtxtnew(FileName);
	}

private:
	void testIndex(const int n1)
	{
		if(0 > n1 || n1 >= clArrayBase1D<T>::Size1)
		{
			printf("Error in 1D array [%d]! Index n1=%d out of range [0...%d]", this->Size1, n1, this->Size1-1);
			exit(0);
		}

		if(n1 < 0 || n1 >= this->FullSize)
		{
			printf("Error in 1D array!");
			exit(0);
		}
	}

	void setPointer(T *p, const int n1, const int s1, const int fSize)
	{
		delArray();

		clArrayBase1D<T>::Array = p;
		clArrayBase1D<T>::flDelete = 0;

		clArrayBase1D<T>::Size1 = n1;
		clArrayBase1D<T>::Step1 = s1;
		clArrayBase1D<T>::FullSize = fSize;
	}

};

typedef  clArray1D <double> clArray1Dd;
typedef  clArray1D <int> clArray1Di;
typedef  clArray1D <signed char> clArray1Dc;
//typedef  clVector2D <int> clVector2Di;

#include "clVector.h"
template <class T> class clArray1Dv2: public clArray1D <clVector2D <T> >
{
public:
	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && savebin(FileName);
	}
	BOOL savetxt_xy(const char * FileName)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		return clArray1Dv2<T>::savetxt(Bufx, Bufy);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		BOOL sx = clArray1Dv2<T>::save2file(Bufx, 0);
		BOOL sy = clArray1Dv2<T>::save2file(Bufy, 1);
		return sx && sy;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray1Dv2<T>::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
			if(fprintf(stream, "%g\n", (double)getEl(i).e(ind)) <= 0)
			{
				fclose(stream);
				return FALSE;
			}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
			if(fprintf(stream, "%g\t%g\n", (double)getEl(i).x, (double)getEl(i).y) <= 0)
			{
				fclose(stream);
				return FALSE;
			}

		fclose(stream);
		return TRUE;
	}
	//load array from a file
	BOOL load(const char * FileName)
	{
		return !((loadbin(FileName) == FALSE) && (loadtxt(FileName) == FALSE));
	}
	BOOL loadtxt(const char * FileName)
	{
		int i;
		double Buf1, Buf2;
		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fscanf(stream, "%lg%lg", &Buf1, &Buf2) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			setEl(clVector2D<T>((T)Buf1, (T)Buf2), i);
		}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt_xy(const char * FileName)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		return clArray1Dv2<T>::loadtxt(Bufx, Bufy);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY)
	{
		BOOL lx = clArray1Dv2<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray1Dv2<T>::loadtxt(FileNameY, 1);
		return lx && ly;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fscanf(stream, "%lg", &Buf) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			el(i).e(ind) = (T)Buf;
		}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxtnew(const char * FileName)
	{
		zeros(counttxt(FileName)/2);
		return loadtxt(FileName);
	}
};

template <class T> class clArray1Dv3: public clArray1D <clVector3D <T> >
{
public:
	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && this->savebin(FileName);
	}
	BOOL savetxt_xyz(const char * FileName)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		sprintf(Bufz, "%sz", FileName);
		return clArray1Dv3<T>::savetxt(Bufx, Bufy, Bufz);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		sprintf(Bufz, "%s.txt", FileNameZ);
		BOOL sx = clArray1Dv3<T>::save2file(Bufx, 0);
		BOOL sy = clArray1Dv3<T>::save2file(Bufy, 1);
		BOOL sz = clArray1Dv3<T>::save2file(Bufz, 2);
		return sx && sy && sz;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray1Dv3<T>::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
			if(fprintf(stream, "%g\n", (double)getEl(i).e(ind)) <= 0)
			{
				fclose(stream);
				return FALSE;
			}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i;
		FILE *stream;
		if((stream = this->fileopen(FileName, "w+")) == NULL)
			return FALSE;
	

		for(i = 0; i < this->Size1; i++)
		{
			clVector3D <T> Vel = this->getEl(i);
			if(fprintf(stream, "%g\t%g\t%g\n", (double)Vel.x, (double)Vel.y, (double)Vel.z) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
		}
		fclose(stream);
		return TRUE;
	}
	//load array from a file
	BOOL load(const char * FileName)
	{
		return !((this->loadbin(FileName) == FALSE) && (loadtxt(FileName) == FALSE));
	}
	BOOL loadtxt(const char * FileName)
	{
		int i;
		double Buf1, Buf2, Buf3;
		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;

		for(i = 0; i < this->Size1; i++)
		{
			if(fscanf(stream, "%lg%lg%lg", &Buf1, &Buf2, &Buf3) <= 0)
			{
			fclose(stream);
				return FALSE;
			}
			this->setEl(clVector3D<T>((T)Buf1, (T)Buf2, (T)Buf3), i);
		}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt_xyz(const char * FileName)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		sprintf(Bufz, "%sz", FileName);
		return clArray1Dv3<T>::loadtxt(Bufx, Bufy, Bufz);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		BOOL lx = clArray1Dv3<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray1Dv3<T>::loadtxt(FileNameY, 1);
		BOOL lz = clArray1Dv3<T>::loadtxt(FileNameZ, 2);
		return lx && ly && lz;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			if(fscanf(stream, "%lg", &Buf) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			this->el(i).e(ind) = (T)Buf;
		}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxtnew(const char * FileName)
	{
		zeros(counttxt(FileName)/3);
		return loadtxt(FileName);
	}
};

template <class T> class clArray1Dts3: public clArray1D <clTensorS3D <T> >
{
public:
	
	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && savebin(FileName);
	}
	BOOL savetxt_xyz(const char * FileName)
	{
		char Bufxx[80],Bufyy[80],Bufzz[80], Bufxy[80],Bufxz[80],Bufyz[80];
		sprintf(Bufxx, "%sxx", FileName);
		sprintf(Bufyy, "%syy", FileName);
		sprintf(Bufzz, "%szz", FileName);
		sprintf(Bufxy, "%sxy", FileName);
		sprintf(Bufxz, "%sxz", FileName);
		sprintf(Bufyz, "%syz", FileName);
		return clArray1Dts3<T>::savetxt(Bufxx, Bufyy, Bufzz, Bufxy, Bufxz, Bufyz);
	}

	BOOL savetxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameZZ, const char * FileNameXY, const char * FileNameXZ, const char * FileNameYZ)
	{
		BOOL sxx = clArray1Dts3<T>::save2file(FileNameXX, 0);
		BOOL syy = clArray1Dts3<T>::save2file(FileNameYY, 1);
		BOOL szz = clArray1Dts3<T>::save2file(FileNameZZ, 2);
		BOOL sxy = clArray1Dts3<T>::save2file(FileNameXY, 3);
		BOOL sxz = clArray1Dts3<T>::save2file(FileNameXZ, 4);
		BOOL syz = clArray1Dts3<T>::save2file(FileNameYZ, 5);
		return sxx && syy && szz && sxy && sxz && syz;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray1Dts3::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		stream = fileopen(FileNameExt, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fprintf(stream, "%g\n", (double)getEl(i).e(ind)) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i;
		FILE *stream;
		stream = fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fprintf(stream, "%g\t%g\t%g\t%g\t%g\t%g\n", 
					(double)getEl(i).xx, (double)getEl(i).yy, (double)getEl(i).zz, 
					(double)getEl(i).xy, (double)getEl(i).xz, (double)getEl(i).yz) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
		}

		fclose(stream);
		return TRUE;
	}
	//load array from a file
	BOOL load(const char * FileName)
	{
		return !((loadbin(FileName) == FALSE) && (loadtxt(FileName) == FALSE));
	}

	BOOL loadtxt_xyz(const char * FileName)
	{
		char Bufxx[80],Bufyy[80],Bufzz[80], Bufxy[80],Bufxz[80],Bufyz[80];
		sprintf(Bufxx, "%sxx", FileName);
		sprintf(Bufyy, "%syy", FileName);
		sprintf(Bufzz, "%szz", FileName);
		sprintf(Bufxy, "%sxy", FileName);
		sprintf(Bufxz, "%sxz", FileName);
		sprintf(Bufyz, "%syz", FileName);
		return clArray1Dts3<T>::loadtxt(Bufxx, Bufyy, Bufzz, Bufxy, Bufxz, Bufyz);
	}

	BOOL loadtxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameZZ, const char * FileNameXY, const char * FileNameXZ, const char * FileNameYZ)
	{
		BOOL sxx = clArray1Dts3<T>::loadtxt(FileNameXX, 0);
		BOOL syy = clArray1Dts3<T>::loadtxt(FileNameYY, 1);
		BOOL szz = clArray1Dts3<T>::loadtxt(FileNameZZ, 2);
		BOOL sxy = clArray1Dts3<T>::loadtxt(FileNameXY, 3);
		BOOL sxz = clArray1Dts3<T>::loadtxt(FileNameXZ, 4);
		BOOL syz = clArray1Dts3<T>::loadtxt(FileNameYZ, 5);
		return sxx && syy && szz && sxy && sxz && syz;
	}
	BOOL loadtxt(const char * FileName)
	{
		int i;
		double Buf1, Buf2, Buf3, Buf4, Buf5, Buf6;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fscanf(stream, "%lg%lg%lg%lg%lg%lg", &Buf1, &Buf2, &Buf3, &Buf4, &Buf5, &Buf6) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			setEl(clTensorS3D<T>((T)Buf1, (T)Buf2, (T)Buf3, (T)Buf4, (T)Buf5, (T)Buf6), i);
		}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase1D<T>::Size1; i++)
		{
			if(fscanf(stream, "%lg", &Buf) <= 0)
			{
				fclose(stream);
				return FALSE;
			}
			el(i).e(ind) = (T)Buf;
		}

		fclose(stream);
		return TRUE;
	}

};

template <class T> class clArray1DLl: public clArray1D <clLinklist <T> >
{
public:
	void clear()
	{
		int i;
		for(i = 0; i < this->Size1; i++)
		{
			this->el(i).~clLinklist();
		}
	}

	clArray1D<int> get_indicies(int i)
	{
		int ind = 0;
		clArray1D<int> index;
		int num = this->el(i).length();
		index.zeros(num);
		clNode<int> *tempt;
		tempt = this->el(i).firstlocal();
		while(tempt != NULL)
		{
			index(ind++) = tempt->getData();
			tempt = tempt->getLink();
		}
		return index;
	}

};

typedef  clArray1DLl <int> clArray1DLli;


typedef  clArray1Dv2 <double> clArray1Dv2d;
typedef  clArray1Dv2 <int> clArray1Dv2i;
typedef  clArray1Dv2 <char> clArray1Dv2c;

typedef  clArray1Dv3 <double> clArray1Dv3d;
typedef  clArray1Dv3 <int> clArray1Dv3i;
typedef  clArray1Dv3 <char> clArray1Dv3c;

typedef  clArray1Dts3 <double> clArray1Dts3d;
typedef  clArray1Dts3 <int> clArray1Dts3i;
typedef  clArray1Dts3 <char> clArray1Dts3c;


#define ARRAY_INCREASE_STEP	4

template <class T> class clArray1Dst: public clArray1D <T>
{
public:
	int nN, nA;
	clArray1Dst()
	{
		ini();
	}

	void ini()
	{
		nN = 0;
		nA = 0;
	}
	void increase(int nIncr)
	{
		clArray1D<T>::reallocate(nA + nIncr);
		nA += nIncr;
	}
	void add(int n)
	{
		if(nN >= nA)
			increase(ARRAY_INCREASE_STEP);
		this->el(nN++) = n;
	}
	void reset()
	{
		nN = 0;
	}
	int find(int n)
	{
		for(int i = 0; i < nN; i++)
			if(this->el(i) == n)
				return i;
		return -1;
	}
//#ifdef ARRAY_DEBUG
//	T& clArray1Dst::operator()(const int n1)
//	{
//		if(n1 >= nN || n1 >= nA || n1 < 0)
//			printf("Error in clArray1Dst, %d from %d tot %d", n1, nN, nA);
//		return clArray1D<T>::operator()(n1);
//	}
//#endif
};

typedef  clArray1Dst <double> clArray1Dstd;
typedef  clArray1Dst <int> clArray1Dsti;
typedef  clArray1Dst <signed char> clArray1Dstc;


#endif // !defined(AFX_ARRAY1D_H__INCLUDED_)
