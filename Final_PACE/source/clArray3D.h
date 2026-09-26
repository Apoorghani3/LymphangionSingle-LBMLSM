// Array3D.h: interface for the Array3D class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAY3D_H__INCLUDED_)
#define AFX_ARRAY3D_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <stdio.h>
#include <memory.h>
#include "clArray1D.h"
#include "clArray2D.h"
#include "clArrayBase3D.h"

template <class T> class clArray1D;
template <class T> class clArray2D;

template <class T> class clArray3D: public clArrayBase3D <T>
{
public:

	clArray3D()
	{
		clArrayBase3D<T>::ini();
		clArrayBase3D<T>::Size3 = 0;
		clArrayBase3D<T>::Step3 = 0;
	}
	clArray3D(const int n1, const int n2, const int n3)
	{
		clArrayBase3D<T>::ini();
		allocate(n1, n2, n3);
	}
	virtual ~clArray3D(){clArrayBase3D<T>::delArray();}


	T& operator()(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		return el(n1, n2, n3);
	}
	T& operator()(clVector3Di n)
	{
		return el(n.x, n.y, n.z);
	}
	T& el(clVector3Di n)
	{
		return el(n.x, n.y, n.z);
	}

	T& el(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		return *(clArrayBase3D<T>::Array + clArrayBase3D<T>::Step1 * n1 + clArrayBase3D<T>::Size3 * n2 + n3);
	}

	T getEl(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		return *(clArrayBase3D<T>::Array + clArrayBase3D<T>::Step1 * n1 + clArrayBase3D<T>::Size3 * n2 + n3);
	}

	T * getPointerEl(const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		return (clArrayBase3D<T>::Array + clArrayBase3D<T>::Step1 * n1 + clArrayBase3D<T>::Size3 * n2 + n3);
	}

	void setEl(const T val, const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		*(clArrayBase3D<T>::Array + clArrayBase3D<T>::Step1 * n1 + clArrayBase3D<T>::Size3 * n2 + n3) = val;
	}

	void addEl(const T val, const int n1, const int n2, const int n3)
	{
#ifdef ARRAY_DEBUG
		TestIndex(n1, n2, n3);
#endif
		*(clArrayBase3D<T>::Array + clArrayBase3D<T>::Step1 * n1 + clArrayBase3D<T>::Size3 * n2 + n3) += val;
	}

	void allocate(const clVector3Di nn)
	{
		allocate(nn.x, nn.y, nn.z);
	}

	void allocate(const int n1, const int n2, const int n3)
	{
		this->delArray();

		this->Size1 = n1;
		this->Size2 = n2;
		this->Size3 = n3;
		this->Step1 = this->Size2 * this->Size3;
		this->Step2 = this->Size3;
		this->Step3 = 1;
		this->FullSize = this->Size1 * this->Size2 * this->Size3;

		this->allocateMemory();
		this->currEl = this->Array;
	}

	//memory reallocation
	void reallocate(const int n1, const int n2, const int n3)
	{
		const int n23 = n2 * n3;
		T * p = allocateMemory(n1 * n23);
		memset(p, 0, sizeof(T) * n1 * n23);

		for(int i = 0; i < this->Size1 && i < n1; i++)
			for(int j = 0; j < this->Size2 && j < n2; j++)
				for(int k = 0; k < this->Size3 && k < n3; k++)
					*(p + n23 * i + n3 * j + 1 * k) = el(i, j, k);
		
		delArray();

		this->Size1 = n1;
		this->Size2 = n2;
		this->Size3 = n3;
		this->Step1 = this->Size2 * this->Size3;
		this->Step2 = this->Size3;
		this->Step3 = 1;
		this->FullSize = this->Size1 * this->Size2 * this->Size3;

		this->Array = p;
		this->flDelete = 1;
	}

	//allocates memory and fills by zeros
	void zeros(const int n1, const int n2, const int n3)
	{
		allocate(n1, n2, n3);
		zeros();
	}
	
	void zeros(const clVector3Di n)
	{
		zeros(n.x, n.y, n.z);
	}
	
	void zeros()
	{
		clArrayBase1D<T>::zeros();
	}

	void copy2Dto3D(const int n1, clArray2D<T> *p2D)
	{
		this->currEl = this->Array + this->Step1 * n1;
		copy2Dto3D(p2D);
	}

	void copy2Dto3D(clArray2D<T> *p2D)
	{
		memcpy(this->currEl, p2D->getArray(), sizeof(double) * p2D->getFullSize());
		this->currEl += this->Step1;
	}

	BOOL save(const char * FileName)
	{
		return savetxt(FileName) && this->savebin(FileName);
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return save2file(Buf);
	}

	BOOL save2file(const char * FileName)
	{
		int i, j, k;
		FILE *stream;
		stream = this->fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
#ifdef ARRAY_DEBUG
//					fprintf(stream, "%18.12e\t", (double)getEl(i, j, k));
					fprintf(stream, "%g\t", (double)getEl(i, j, k));
#else
					fprintf(stream, "%g\t", (double)getEl(i, j, k));
#endif
				}
				fprintf(stream, "\n");
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
		int i, j, k;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < clArrayBase3D<T>::Size1; i++)
			for(j = 0; j < clArrayBase3D<T>::Size2; j++)
				for(k = 0; k < clArrayBase3D<T>::Size3; k++)
				{
					if(fscanf(stream, "%lg", &Buf) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					setEl((T)Buf, i, j, k);
				}

		fclose(stream);
		return TRUE;
	}

private:
	void TestIndex(const int n1, const int n2, const int n3)
	{
		if(0 > n1 || n1 >= this->Size1)
		{
			printf("Error in 3D array [%d, %d, %d]! Index n1=%d out of range [0...%d]", this->Size1, this->Size2, this->Size3, n1, this->Size1-1);
			exit(0);
		}

		if(0 > n2 || n2 >= this->Size2)
		{
			printf("Error in 3D array [%d, %d, %d]! Index n2=%d out of range [0...%d]", this->Size1, this->Size2, this->Size3, n2, this->Size2-1);
			exit(0);
		}

		if(0 > n3 || n3 >= this->Size3)
		{
			printf("Error in 3D array [%d, %d, %d]! Index n3=%d out of range [0...%d]", this->Size1, this->Size2, this->Size3, n3, this->Size3-1);
			exit(0);
		}

		if(this->Step1 * n1 + this->Size3 * n2 + n3 < 0 || this->Step1 * n1 + this->Size3 * n2 + n3 >= this->FullSize)
		{
			printf("Error in 3D array!");
			exit(0);
		}
	}
};



typedef  clArray3D <double> clArray3Dd;
typedef  clArray3D <int> clArray3Di;
typedef  clArray3D <char> clArray3Dc;


#include "clVector.h"
template <class T> class clArray3Dv2: public clArray3D <clVector2D <T> >
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
		return clArray3Dv2<T>::savetxt(Bufx, Bufy);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		BOOL sx = clArray3Dv2<T>::save2file(Bufx, 0);
		BOOL sy = clArray3Dv2<T>::save2file(Bufy, 1);
		return sx && sy;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray3Dv2::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j, k;
		FILE *stream;
		stream = fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t", 
							(double)getEl(i, j, k).e(ind)) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i, j, k;
		FILE *stream;
		stream = fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t%g\t", 
							(double)getEl(i, j, k).x, (double)getEl(i, j, k).y) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
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

	BOOL loadtxt_xy(const char * FileName)
	{
		char Bufx[80],Bufy[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		return clArray3Dv2<T>::loadtxt(Bufx, Bufy);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY)
	{
		BOOL lx = clArray3Dv2<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray3Dv2<T>::loadtxt(FileNameY, 1);
		return lx && ly;
	}

	BOOL loadtxt(const char * FileName)
	{
		int i, j, k;
		double Buf1, Buf2;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg%lg", &Buf1, &Buf2) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					setEl(clVector3D<T>((T)Buf1, (T)Buf2), i, j, k);
				}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j, k;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg", &Buf) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					el(i,j,k).e(ind) = (T)Buf;
				}

		fclose(stream);
		return TRUE;
	}

};

template <class T> class clArray3Dv3: public clArray3D <clVector3D <T> >
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
		return clArray3Dv3<T>::savetxt(Bufx, Bufy, Bufz);
	}

	BOOL savetxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%s.txt", FileNameX);
		sprintf(Bufy, "%s.txt", FileNameY);
		sprintf(Bufz, "%s.txt", FileNameZ);
		BOOL sx = clArray3Dv3<T>::save2file(Bufx, 0);
		BOOL sy = clArray3Dv3<T>::save2file(Bufy, 1);
		BOOL sz = clArray3Dv3<T>::save2file(Bufz, 2);
		return sx && sy && sz;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray3Dv3::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j, k;
		FILE *stream;
		stream = fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t", 
							(double)getEl(i, j, k).e(ind)) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i, j, k;
		FILE *stream;
		stream = this->fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t%g\t%g\t", 
							(double)this->getEl(i, j, k).x, (double)this->getEl(i, j, k).y, (double)this->getEl(i, j, k).z) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
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

	BOOL loadtxt_xyz(const char * FileName)
	{
		char Bufx[80],Bufy[80],Bufz[80];
		sprintf(Bufx, "%sx", FileName);
		sprintf(Bufy, "%sy", FileName);
		sprintf(Bufz, "%sz", FileName);
		return clArray3Dv3<T>::loadtxt(Bufx, Bufy, Bufz);
	}

	BOOL loadtxt(const char * FileNameX, const char * FileNameY, const char * FileNameZ)
	{
		BOOL lx = clArray3Dv3<T>::loadtxt(FileNameX, 0);
		BOOL ly = clArray3Dv3<T>::loadtxt(FileNameY, 1);
		BOOL lz = clArray3Dv3<T>::loadtxt(FileNameZ, 2);
		return lx && ly && lz;
	}

	BOOL loadtxt(const char * FileName)
	{
		int i, j, k;
		double Buf1, Buf2, Buf3;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg%lg%lg", &Buf1, &Buf2, &Buf3) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					setEl(clVector3D<T>((T)Buf1, (T)Buf2, (T)Buf3), i, j, k);
				}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j, k;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg", &Buf) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					el(i,j,k).e(ind) = (T)Buf;
				}

		fclose(stream);
		return TRUE;
	}

};

template <class T> class clArray3Dts3: public clArray3D <clTensorS3D <T> >
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
		return clArray3Dts3<T>::savetxt(Bufxx, Bufyy, Bufzz, Bufxy, Bufxz, Bufyz);
	}

	BOOL savetxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameZZ, const char * FileNameXY, const char * FileNameXZ, const char * FileNameYZ)
	{
		BOOL sxx = clArray3Dts3<T>::save2file(FileNameXX, 0);
		BOOL syy = clArray3Dts3<T>::save2file(FileNameYY, 1);
		BOOL szz = clArray3Dts3<T>::save2file(FileNameZZ, 2);
		BOOL sxy = clArray3Dts3<T>::save2file(FileNameXY, 3);
		BOOL sxz = clArray3Dts3<T>::save2file(FileNameXZ, 4);
		BOOL syz = clArray3Dts3<T>::save2file(FileNameYZ, 5);
		return sxx && syy && szz && sxy && sxz && syz;
	}

	BOOL savetxt(const char * FileName)
	{
		char Buf[80];
		sprintf(Buf, "%s.txt", FileName);
		return clArray3Dts3::save2file(Buf);
	}

	BOOL save2file(const char * FileName, const int ind)
	{
		int i, j, k;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		stream = fileopen(FileNameExt, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t", 
							(double)getEl(i, j, k).e(ind)) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
			}
			fprintf(stream, "\n");
		}

		fclose(stream);
		return TRUE;
	}

	BOOL save2file(const char * FileName)
	{
		int i, j, k;
		FILE *stream;
		stream = this->fileopen(FileName, "w+");
		if(stream == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
		{
			for(j = 0; j < this->Size2; j++)
			{
				for(k = 0; k < this->Size3; k++)
				{
					if(fprintf(stream, "%g\t%g\t%g\t%g\t%g\t%g\t", 
							(double)this->getEl(i, j, k).xx, (double)this->getEl(i, j, k).yy, (double)this->getEl(i, j, k).zz, 
							(double)this->getEl(i, j, k).xy, (double)this->getEl(i, j, k).xz, (double)this->getEl(i, j, k).yz) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
				}
				fprintf(stream, "\n");
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

	BOOL loadtxt_xyz(const char * FileName)
	{
		char Bufxx[80],Bufyy[80],Bufzz[80], Bufxy[80],Bufxz[80],Bufyz[80];
		sprintf(Bufxx, "%sxx", FileName);
		sprintf(Bufyy, "%syy", FileName);
		sprintf(Bufzz, "%szz", FileName);
		sprintf(Bufxy, "%sxy", FileName);
		sprintf(Bufxz, "%sxz", FileName);
		sprintf(Bufyz, "%syz", FileName);
		return clArray3Dts3<T>::loadtxt(Bufxx, Bufyy, Bufzz, Bufxy, Bufxz, Bufyz);
	}

	BOOL loadtxt(const char * FileNameXX, const char * FileNameYY, const char * FileNameZZ, const char * FileNameXY, const char * FileNameXZ, const char * FileNameYZ)
	{
		BOOL sxx = clArray3Dts3<T>::loadtxt(FileNameXX, 0);
		BOOL syy = clArray3Dts3<T>::loadtxt(FileNameYY, 1);
		BOOL szz = clArray3Dts3<T>::loadtxt(FileNameZZ, 2);
		BOOL sxy = clArray3Dts3<T>::loadtxt(FileNameXY, 3);
		BOOL sxz = clArray3Dts3<T>::loadtxt(FileNameXZ, 4);
		BOOL syz = clArray3Dts3<T>::loadtxt(FileNameYZ, 5);
		return sxx && syy && szz && sxy && sxz && syz;
	}
	BOOL loadtxt(const char * FileName)
	{
		int i, j, k;
		double Buf1, Buf2, Buf3, Buf4, Buf5, Buf6;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg%lg%lg%lg%lg%lg", &Buf1, &Buf2, &Buf3, &Buf4, &Buf5, &Buf6) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					setEl(clTensorS3D<T>((T)Buf1, (T)Buf2, (T)Buf3, (T)Buf4, (T)Buf5, (T)Buf6), i, j, k);
				}

		fclose(stream);
		return TRUE;
	}

	BOOL loadtxt(const char * FileName, const int ind)
	{
		int i, j, k;
		double Buf;

		char FileNameExt[80];
		sprintf(FileNameExt, "%s.txt", FileName);

		FILE *stream;
		if((stream = fopen(FileNameExt, "r")) == NULL)
			return FALSE;
	
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					if(fscanf(stream, "%lg", &Buf) <= 0)
					{
						fclose(stream);
						return FALSE;
					}
					el(i,j,k).e(ind) = (T)Buf;
				}

		fclose(stream);
		return TRUE;
	}

};

typedef  clArray3Dv2 <double> clArray3Dv2d;
typedef  clArray3Dv2 <int> clArray3Dv2i;
typedef  clArray3Dv2 <char> clArray3Dv2c;

typedef  clArray3Dv3 <double> clArray3Dv3d;
typedef  clArray3Dv3 <int> clArray3Dv3i;
typedef  clArray3Dv3 <char> clArray3Dv3c;

typedef  clArray3Dts3 <double> clArray3Dts3d;
typedef  clArray3Dts3 <int> clArray3Dts3i;
typedef  clArray3Dts3 <char> clArray3Dts3c;


template <class T> class clArray3DLl: public clArray3D <clLinklist <T> >
{
public:
	void clear()
	{
		int i, j, k;
		for(i = 0; i < this->Size1; i++)
			for(j = 0; j < this->Size2; j++)
				for(k = 0; k < this->Size3; k++)
				{
					this->el(i,j,k).~clLinklist();
				}
	}

	clArray1D<int> get_indicies(int i, int j, int k) 
	{
		int ind = 0;
		clArray1D<int> index;
		int num = this->el(i,j,k).length();
		index.zeros(num);
		clNode<int> *tempt;
		tempt = this->el(i,j,k).firstlocal();
		while(tempt != NULL)
		{
			index(ind++) = tempt->getData();
			tempt = tempt->getLink();
		}
		return index;
	}

};

typedef  clArray3DLl <int> clArray3DLli;


#endif // !defined(AFX_ARRAY3D_H__INCLUDED_)
