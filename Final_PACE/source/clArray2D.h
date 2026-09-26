// Array2D.h: interface for the Array2D class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAY2D_H__INCLUDED_)
#define AFX_ARRAY2D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <stdio.h>
#include <memory.h>
#include "clArray1D.h"
#include "clArray3D.h"
#include "clArrayBase2D.h"

template <class T> class clArray1D;
template <class T> class clArray3D;

template <class T> class clArray2D: public clArrayBase2D <T>
{
public:

//////////////////////////////////////////////////////////////////////
// Construction
//////////////////////////////////////////////////////////////////////

	clArray2D()
	{
		clArrayBase2D<T>::ini();
		clArrayBase2D<T>::Size2 = 0;
		clArrayBase2D<T>::Step2 = 0;
	}

	clArray2D(const T *p, const int n1, const int n2)
	{
		ini();
		setPointer(p, n1, n2);
	}

	clArray2D(const int n1, const int n2)
	{
		ini();
		allocate(n1, n2);
	}

//////////////////////////////////////////////////////////////////////
// Destruction
//////////////////////////////////////////////////////////////////////
	virtual ~clArray2D(){clArrayBase2D<T>::delArray();}

	T& operator()(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		return *(clArrayBase2D<T>::Array + clArrayBase2D<T>::Step1 * n1 + clArrayBase2D<T>::Step2 * n2);
	}

	T& operator()(clVector2Di n)
	{
		return operator()(n.x, n.y);
	}

	T& el(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		return *(clArrayBase2D<T>::Array + clArrayBase2D<T>::Step1 * n1 + clArrayBase2D<T>::Step2 * n2);
	}

	T getEl(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		return *(clArrayBase2D<T>::Array + clArrayBase2D<T>::Step1 * n1 + clArrayBase2D<T>::Step2 * n2);
	}

	T *getPointerEl(const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		return (clArrayBase2D<T>::Array + clArrayBase2D<T>::Step1 * n1 + clArrayBase2D<T>::Step2 * n2);
	}


	void setEl(const T val, const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		*(clArrayBase2D<T>::Array + clArrayBase2D<T>::Step1 * n1 + clArrayBase2D<T>::Step2 * n2) = val;
	}

	void addEl(const T val, const int n1, const int n2)
	{
#ifdef ARRAY_DEBUG
		testIndex(n1, n2);
#endif
		*(this->Array + this->Step1 * n1 + this->Step2 * n2) += val;
	}


	void map2Don3D2(clArray3D<T> &p3D, const int n1)
	{
		setPointer(p3D.getPointerEl(n1, 0, 0), 
					p3D.getSize2(),  p3D.getSize3(), p3D.getFullSize());
	}

	void map2Don3D1(clArray3D<T> &p3D, const int n2)
	{
		setPointer(p3D.getPointerEl(0, n2, 0), 
					p3D.getSize1(), p3D.getSize3(), 
					p3D.getSize2() * p3D.getSize3(), 1, p3D.getFullSize());
	}

	void	setPointer(T *p, const int n1, const int n2)
	{
		setPointer(p, n1, n2, n2, 1, n1 * n2);
	}

	void setPointer(const T *p, const int n1, const int n2, const int fSize)
	{
		setPointer(p, n1, n2, n2, 1, fSize);
	}

	void setPointer(T *p, const int n1, const int n2, const int s1, const int s2, const int fSize)
	{
		this->delArray();

		this->Array = p;
		this->flDelete = 0;

		this->Size1 = n1;
		this->Size2 = n2;
		this->Step1 = s1;
		this->Step2 = s2;
		this->FullSize = fSize;
	}

	void allocate(const clVector2Di nn)
	{
		allocate(nn.x, nn.y);
	}

	void allocate(const int n1, const int n2)
	{
		setPointer(NULL, n1, n2);

		this->allocateMemory();
	}

	//memory reallocation
	void reallocate(const int n1, const int n2)
	{
		T * p = clArrayBase2D<T>::allocateMemory(n1 * n2);
		memset(p, 0, sizeof(T) * n1 * n2);

		for(int i = 0; i < clArrayBase2D<T>::Size1 && i < n1; i++)
			for(int j = 0; j < clArrayBase2D<T>::Size2 && j < n2; j++)
				*(p + n2 * i + 1 * j) = el(i, j);
		
		clArrayBase2D<T>::delArray();

		setPointer(p, n1, n2);
		clArrayBase2D<T>::flDelete = 1;
	}

	//allocates memory and fills by zeros
	void zeros(const int n1, const int n2)
	{
		allocate(n1, n2);
		zeros();
	}

	void zeros(const clVector2Di n)
	{
		zeros(n.x, n.y);
	}

	void zeros()
	{
		clArrayBase1D<T>::zeros();
	}

	int findnonzero(int n, int m)
	{
		for(int i = n; i < clArrayBase2D<T>::Size1; i++)
			if(this->el(i, m) != 0)
				return i;
		return -1;
	}

	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && clArrayBase2D<T>::savebin(FileName);
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return save2file(Buf);
	}


	BOOL save2file(const char * FileName)
	{
		int i, j;
		FILE *stream;
		if((stream = clArrayBase2D<T>::fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase2D<T>::Size1; i++)
		{
			for(j = 0; j < clArrayBase2D<T>::Size2; j++)
			{
#ifdef ARRAY_DEBUG
//				if(fprintf(stream, "%18.12e\t", (double)getEl(i, j)) <= 0)
				if(fprintf(stream, "%g\t", (double)getEl(i, j)) <= 0)
#else
				if(fprintf(stream, "%g\t", (double)getEl(i, j)) <= 0)
#endif
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
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
		int i, j;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase2D<T>::Size1; i++)
			for(j = 0; j < clArrayBase2D<T>::Size2; j++)
			{
				if(fscanf(stream, "%lg", &Buf) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				setEl((T)Buf, i, j);
			}

		fclose(stream);
		return TRUE;
	}
	void zerosboundary(void)
	{
		int i, j;
		for(i = 0; i < this->Size1; i++)
		{
			setEl(0, i, 0);
			setEl(0, i, this->Size2-1);
		}
		for(j = 0; j < this->Size2; j++)
		{
			setEl(0, 0, j);
			setEl(0, this->Size1-1, j);
		}

		return;
	}
	BOOL loadtxtnew(const char * FileName, const int s2)
	{
		zeros(this->counttxt(FileName) / s2, s2);
		return loadtxt(FileName);
	}

	BOOL loadbinnew(const char * FileName, const int s2)
	{
		zeros(countbin(FileName) / s2, s2);
		return loadbin(FileName);
	}

	BOOL loadnew(const char * FileName, const int s2)
	{
		int counter = clArrayBase2D<T>::countbin(FileName);
		if(counter > 0)
		{
			zeros(counter / s2, s2);
			return clArrayBase2D<T>::loadbin(FileName);
		}
		return loadtxtnew(FileName, s2);
	}
private:
	void testIndex(const int n1, const int n2)
	{
		if(0 > n1 || n1 >= this->Size1)
		{
			printf("Error in 2D array [%d, %d]! Index n1=%d out of range [0...%d]", this->Size1, this->Size2, n1, this->Size1-1);
			exit(0);
		}

		if(0 > n2 || n2 >= this->Size2)
		{
			printf("Error in 2D array [%d, %d]! Index n2=%d out of range [0...%d]", this->Size1, this->Size2, n2, this->Size2-1);
			exit(0);
		}

		if(this->Step1 * n1 + this->Step2 * n2 < 0 || this->Step1 * n1 + this->Step2 * n2 >= this->FullSize)
		{
			printf("Error in 2D array!");
			exit(0);
		}
	}


};

typedef  clArray2D <double> clArray2Dd;
typedef  clArray2D <int> clArray2Di;
typedef  clArray2D <signed char> clArray2Dc;

#include "clVector.h"
template <class T> class clArray2Dv2: public clArray2D <clVector2D <T> >
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
		return clArray2Dv2<T>::savetxt(Bufx, Bufy);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		BOOL sx = clArray2Dv2<T>::save2file(Bufx, 0);
		BOOL sy = clArray2Dv2<T>::save2file(Bufy, 1);
		return sx && sy;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray2Dv2<T>::save2file(Buf);
	}

	BOOL save2file(const char * FileName)
	{
		int i, j;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				if(fprintf(stream, "%g\t%g\t", 
							(double)getEl(i, j).x, (double)getEl(i, j).y) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				if(fprintf(stream, "%g\t", (double)getEl(i, j).e(ind)) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
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
		int i, j;
		double Buf1, Buf2;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg%lg", &Buf1, &Buf2) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				setEl(clVector2D<T>((T)Buf1, (T)Buf2), i, j);
			}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt_xy(const char * FileName)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		return clArray2Dv2<T>::loadtxt(Bufx, Bufy);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY)
	{
		BOOL lx = clArray2Dv2<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray2Dv2<T>::loadtxt(FileNameY, 1);
		return lx && ly;
	}
	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg", &Buf) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				el(i,j).e(ind) = (T)Buf;
			}

		fclose(stream);
		return TRUE;
	}
};


template <class T> class clArray2Dts2: public clArray2D <clTensorS2D <T> >
{
public:
	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && savebin(FileName);
	}
	BOOL savetxt_xy(const char * FileName)
	{
		char Bufxx[80],Bufyy[80],Bufxy[80];
		sprintf(Bufxx, "%sxx", FileName);
		sprintf(Bufyy, "%syy", FileName);
		sprintf(Bufxy, "%sxy", FileName);
		return clArray2Dts2<T>::savetxt(Bufxx, Bufyy, Bufxy);
	}

	BOOL savetxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameXY)
	{
		char Bufxx[80],Bufyy[80],Bufxy[80];
		sprintf(Bufxx, "%s.txt", FileNameXX);
		sprintf(Bufyy, "%s.txt", FileNameYY);
		sprintf(Bufxy, "%s.txt", FileNameXY);
		BOOL sxx = clArray2Dts2<T>::save2file(Bufxx, 0);
		BOOL syy = clArray2Dts2<T>::save2file(Bufyy, 1);
		BOOL sxy = clArray2Dts2<T>::save2file(Bufxy, 2);
		return sxx && syy && sxy;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray2Dts2<T>::save2file(Buf);
	}

	BOOL save2file(const char * FileName)
	{
		int i, j;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				if(fprintf(stream, "%g\t%g\t%g\t", 
							(double)getEl(i, j).xx, (double)getEl(i, j).yy, (double)getEl(i, j).xy) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				if(fprintf(stream, "%g\t", (double)getEl(i, j).e(ind)) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
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
		int i, j;
		double Buf1, Buf2, Buf3;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg%lg%lg", &Buf1, &Buf2, &Buf3) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				setEl(clTensorS2D<T>((T)Buf1, (T)Buf2, (T)Buf3), i, j);
			}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt_xy(const char * FileName)
	{
		char Bufxx[80],Bufyy[80],Bufxy[80];
		sprintf(Bufxx, "%sxx", FileName);
		sprintf(Bufyy, "%syy", FileName);
		sprintf(Bufxy, "%sxy", FileName);
		return clArray1Dv2<T>::loadtxt(Bufxx, Bufyy, Bufxy);
	}

	BOOL loadtxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameXY)
	{
		BOOL lxx = clArray2Dv2<T>::loadtxt(FileNameXX, 0);
		BOOL lyy = clArray2Dv2<T>::loadtxt(FileNameYY, 1);
		BOOL lxy = clArray2Dv2<T>::loadtxt(FileNameXY, 2);
		return lxx && lyy && lxy;
	}
	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg", &Buf) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				el(i,j).e(ind) = (T)Buf;
			}

		fclose(stream);
		return TRUE;
	}
};

template <class T> class clArray2Dv3: public clArray2D <clVector3D <T> >
{
public:

	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && savebin(FileName);
	}
	BOOL savetxt_xyz(const char * FileName)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		sprintf(Bufz, "%sz", FileName);
		return clArray2Dv3<T>::savetxt(Bufx, Bufy, Bufz);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		sprintf(Bufz, "%s.txt", FileNameZ);
		BOOL sx = clArray2Dv3<T>::save2file(Bufx, 0);
		BOOL sy = clArray2Dv3<T>::save2file(Bufy, 1);
		BOOL sz = clArray2Dv3<T>::save2file(Bufz, 2);
		return sx && sy && sz;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray2Dv3<T>::save2file(Buf);
	}

	BOOL save2file(const char * FileName)
	{
		int i, j;
		FILE *stream;
		if((stream = this->fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				clVector3D <T> Vel = this->getEl(i, j);
				if(fprintf(stream, "%g\t%g\t%g\t", (double)Vel.x, (double)Vel.y, (double)Vel.z) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j;
		FILE *stream;
		if((stream = fileopen(FileName, "w+")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				if(fprintf(stream, "%g\t", (double)getEl(i, j).e(ind)) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
			}
			fprintf(stream, "\n");
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
		int i, j;
		double Buf1, Buf2, Buf3;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg%lg", &Buf1, &Buf2, &Buf3) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				setEl(clVector2D<T>((T)Buf1, (T)Buf2, (T)Buf3), i, j);
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
		return clArray2Dv3<T>::loadtxt(Bufx, Bufy, Bufz);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		BOOL lx = clArray2Dv3<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray2Dv3<T>::loadtxt(FileNameY, 1);
		BOOL lz = clArray2Dv3<T>::loadtxt(FileNameZ, 2);
		return lx && ly && lz;
	}
	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
			{
				if(fscanf(stream, "%lg", &Buf) <= 0)
				{
					fclose(stream);
					return FALSE;
				}
				el(i,j).e(ind) = (T)Buf;
			}

		fclose(stream);
		return TRUE;
	}
};

typedef  clArray2Dv2 <double> clArray2Dv2d;
typedef  clArray2Dv2 <int> clArray2Dv2i;
typedef  clArray2Dv2 <char> clArray2Dv2c;

typedef  clArray2Dts2 <double> clArray2Dts2d;
typedef  clArray2Dts2 <int> clArray2Dts2i;
typedef  clArray2Dts2 <char> clArray2Dts2c;

typedef  clArray2Dv3 <double> clArray2Dv3d;
typedef  clArray2Dv3 <int> clArray2Dv3i;
typedef  clArray2Dv3 <char> clArray2Dv3c;


template <class T> class clArray2DLl: public clArray2D <clLinklist <T> >
{
public:
	void clear()
	{
		int i, j;
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
		{
			this->el(i,j).~clLinklist();
		}
	}

	clArray1D<int> get_indicies(int i, int j)
	{
		int ind = 0;
		clArray1D<int> index;
		int num = this->el(i,j).length();
		index.zeros(num);
		clNode<int> *tempt;
		tempt = this->el(i,j).firstlocal();
		while(tempt != NULL)
		{
			index(ind++) = tempt->getData();
			tempt = tempt->getLink();
		}
		return index;
	}

};

typedef  clArray2DLl <int> clArray2DLli;

#endif // !defined(AFX_ARRAY2D_H__INCLUDED_)
