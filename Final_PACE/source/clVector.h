// clVector.h: interface for the clVector class.
//
// (c) Alexander Alexeev, 2005 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__INCLUDED_)
#define AFX_VECTOR_H__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
//#include "StdAfx.h"

template <class T> class clNode
{
public:
	clNode(){}
    
	clNode(T theData, clNode<T>* theLink) : data(theData), link(theLink){}
    clNode<T>* getLink( ) const { return link; }
    T getData() const { return data; }
    void setData(const T& theData) { data = theData; }
    void setLink(clNode<T>* pointer) { link = pointer; }
private:
    T data;
    clNode<T> *link;
};

template <class T> class clLinklist
{
public:
	void initializeList();
	//function to initialize first element in list to zero
	//postcondition: first=NULL last=NULL

	bool isEmpty();
	//function to see if list contains elements

	int Length();
	//function to return the number of elements initialized in the list

	T firstItem();
	//function to return the first item in the list

	clNode<T>* firstlocal();

	T lastItem();
	//function to return the last item in the list

	void appendItem(const T& newItem);
	//function to append an item to the list

	void destroyList();

	void printList();

	int length();

	//void get_indicies();

	clLinklist();
	//constructor

	~clLinklist();
	//destructor

protected:
	clNode<T> *first;	//pointer to first node
	clNode<T> *last;	//pointer to last node
	int count;	//stores number of items in the list
};

template <class T> bool clLinklist<T>::isEmpty()
{
	return(first==NULL);
}

template <class T> clLinklist<T>::clLinklist()
{
	first=NULL;
	last=NULL;
	count=0;
}
template <class T> clLinklist<T>::~clLinklist()
{
	clNode<T> *temp;						//pointer to deallocate the 
									//  memory occupied by the node
	while(first != NULL)			//while there are nodes in the list
	{
		temp = first;				//set temp to the current node
		first = first->getLink();	//advance first to the next node
		delete temp;				//deallocate the memory occupied by temp
	}
	last = NULL;	//iniitialize last to NULL
					//  (first has already been set to NULL
					//		by the while loop)
	count = 0;
}


template<class T> T clLinklist<T>::firstItem()
{
	if(!isEmpty())
	{
		return first->getData();			
		//return first->data;
		//NB: return the first node's data member eg: return any object
	}
	else
	{
		return NULL;
	}
}

template<class T> clNode<T>* clLinklist<T>::firstlocal()
{
	if(!isEmpty())
	{
		return first;			
		//return first->data;
		//NB: return the first node's data member eg: return any object
	}
	else
	{
		return NULL;
	}
}

template<class T> T clLinklist<T>::lastItem()
{
	if(!isEmpty())
	{
		return last->getData();			
		//NB: return the first node's data member eg: return any object
	}
	else
	{
		return NULL;
	}
}
template <class T> void clLinklist<T>::printList()
{
	clNode<T> *tempt;
	tempt = first;
	while(tempt != NULL)
	{
		printf("Data %d\n", tempt->getData());
		tempt = tempt->getLink();
	}
}

template <class T> int clLinklist<T>::length()
{
	return count;
}

//template <class T> void clLinklist<T>::get_indicies()
//{
//	////int ind = 0;
//	//::clArray1Di<int> index;
//	//index.zeros(count);
//	////clNode<T> *tempt;
//	////tempt = first;
//	////while(tempt != NULL)
//	////{
//	////	index(ind) = tempt->getData();
//	////	tempt = tempt->getLink();
//	////	ind++;
//	////}
//	//return index;
//}


template<class T> void clLinklist<T>::appendItem(const T& newData)
{
	clNode<T> *newNode;				//pointer to create the new node	
	newNode = new clNode<T>;			//create the new node	
	newNode->setData(newData);	//store new data in the node	
	newNode->setLink(NULL);		//insert new node at the end pointing to NULL	
	if (first == NULL)			//if the list is empty, newNode is
								//  both the first and last node
	{
		first = newNode;
		last = newNode;
		++count;
	}
	else						//the list is not empty, insert
								//  newNode after the last node
	{
		last->setLink(newNode); //insert newNode after last
		last = newNode;			//make last point to the actual last node
		++count;
	}
}

template <class T> class clVector2D
{
public:
	T x, y;

	clVector2D(const T x_, const T y_)
	{
		ini(x_, y_);
	}
	clVector2D()
	{
		ini();
	}

	void ini()
	{
		ini(0, 0);
	};

	void ini(const T x_, const T y_)
	{
		x = x_;
		y = y_;
	};

	T nx()
	{
		return x / r();
	}
	T ny()
	{
		return y / r();
	}
	clVector2D<T> n()
	{
		T rr = r();
		return clVector2D<T>(x / rr, y / rr);
	}
	T& e(const int ind)
	{
		switch(ind)
		{
		case 0:
			return x;
		}
		return y;
	}

	T r2()
	{
		return (x*x + y*y);
	}
	double r()
	{
		return sqrt(r2());
	}
	double abs()
	{
		return r();
	}
	clVector2D<T> sort()
	{
		if(x < y)
			return clVector2D<T>(x, y);
		return clVector2D<T>(y, x);
	}
	operator clVector2D<double> ()
	{
		return clVector2D<double>((double)x, (double)y);
	}
	operator clVector2D<int> ()
	{
		return clVector2D<int>((int)x, (int)y);
	}
	clVector2D<T>& operator= (const clVector2D<T>& src)
	{
		x = src.x;
		y = src.y;
		return *this;
	}
	clVector2D<T>& operator= (const double src)
	{
		x = y = src;
		return *this;
	}
	clVector2D<T>& operator= (const int src)
	{
		x = y = src;
		return *this;
	}
	clVector2D<T>& operator+= (const clVector2D<T>& src)
	{
		x += src.x;
		y += src.y;
		return *this;
	}
	clVector2D<T>& operator+= (const double src)
	{
		x += src;
		y += src;
		return *this;
	}
	clVector2D<T>& operator+= (const int src)
	{
		x += src;
		y += src;
		return *this;
	}
	clVector2D<T>& operator-= (const clVector2D<T>& src)
	{
		x -= src.x;
		y -= src.y;
		return *this;
	}
	clVector2D<T>& operator-= (const double src)
	{
		x -= src;
		y -= src;
		return *this;
	}
	clVector2D<T>& operator-= (const int src)
	{
		x -= src;
		y -= src;
		return *this;
	}
	clVector2D<T>& operator*= (const double src)
	{
		x *= src;
		y *= src;
		return *this;
	}
	clVector2D<T>& operator*= (const int src)
	{
		x *= src;
		y *= src;
		return *this;
	}
	clVector2D<T>& operator/= (const double src)
	{
		x /= src;
		y /= src;
		return *this;
	}
	clVector2D<T>& operator/= (const int src)
	{
		x /= src;
		y /= src;
		return *this;
	}

	clVector2D<T> operator+ (const clVector2D<T>& src)
	{
		return clVector2D<T>(x + src.x, y + src.y);
	}
	clVector2D<T> operator+ (const double src)
	{
		return clVector2D<double>(x + src, y + src);
	}
	clVector2D<T> operator+ (const int src)
	{
		return clVector2D<T>(x + src, y + src);
	}

	clVector2D<T> operator- (const clVector2D<T>& src)
	{
		return clVector2D<T>(x - src.x, y - src.y);
	}
	clVector2D<T> operator- (const double src)
	{
		return clVector2D<double>(x - src, y - src);
	}
	clVector2D<T> operator- (const int src)
	{
		return clVector2D<T>(x - src, y - src);
	}
	clVector2D<T> operator- ()
	{
		return clVector2D<T>(-x, -y);
	}

	T operator* (const clVector2D<T>& src)
	{
		return (x * src.x + y * src.y);
	}
	clVector2D<T> operator* (const double src)
	{
		return clVector2D<double>(x * src, y * src);
	}
	clVector2D<T> operator* (const int src)
	{
		return clVector2D<T>(x * src, y * src);
	}
	T operator^ (const clVector2D<T>& src)
	{
		return (x * src.y - y * src.x);
	}
	clVector2D<T> operator/ (const double src)
	{
		return clVector2D<double>(x / src, y / src);
	}
	clVector2D<T> operator/ (const int src)
	{
		return clVector2D<T>(x / src, y / src);
	}

	clVector2D<T> operator% (const clVector2D<int>& src)
	{
		return clVector2D<int>((int)x % src.x, (int)y % src.y);
	}
	clVector2D<T> operator% (const int src)
	{
		return clVector2D<int>((int)x % src, (int)y % src);
	}
	int operator== (const clVector2D<T>& src)
	{
		return x == src.x && y == src.y;
	}
	int operator== (const double src)
	{
		return x == src && y == src;
	}
	int operator== (const int src)
	{
		return x == src && y == src;
	}
	int operator!= (const clVector2D<T>& src)
	{
		return x != src.x || y != src.y;
	}
	int operator!= (const double src)
	{
		return x != src || y != src;
	}
	int operator!= (const int src)
	{
		return x != src || y != src;
	}

};

template <class T> class clVector3D
{
public:
	T x, y, z;

	clVector3D(const T x_, const T y_, const T z_)
	{
		ini(x_, y_, z_);
	}
	clVector3D()
	{
		ini();
	}

	void ini()
	{
		ini(0, 0, 0);
	};

	void ini(const T x_, const T y_, const T z_)
	{
		x = x_;
		y = y_;
		z = z_;
	};

	T& e(const int ind)
	{
		switch(ind)
		{
		case 0:
			return x;
		case 1:
			return y;
		}
		return z;
	}
	clVector3D<T> n()
	{
		T rr = r();
		if(rr == (T)0)
			return clVector3D<T>(0, 0, 0);
		return clVector3D<T>(x / rr, y / rr, z / rr);
	}
	T nx()
	{
		return x / r();
	}
	T ny()
	{
		return y / r();
	}
	T nz()
	{
		return z / r();
	}
	T x2()
	{
		return x*x;
	}
	T y2()
	{
		return y*y;
	}
	T z2()
	{
		return z*z;
	}

	T r2()
	{
		return (x*x + y*y + z*z);
	}
	double r()
	{
		return sqrt(r2());
	}
	double abs()
	{
		return r();
	}

	clVector3D<T> floor()
	{
		return clVector3D<T>((T)::floor(x), (T)::floor(y), (T)::floor(z));
	}
	clVector3D<T> ceil()
	{
		return clVector3D<T>((T)::ceil(x), (T)::ceil(y), (T)::ceil(z));
	}

	clVector3D<T> round()
	{
		clVector3D<T> ret = floor();
		if(x - ret.x > 0.5) ret.x++;
		if(y - ret.y > 0.5) ret.y++;
		if(z - ret.z > 0.5) ret.z++;
		return ret;
	}

	clVector3D<T> fmod(const clVector3D<T>& src)
	{
		clVector3D<T> ret;
		ret.x = (T)::fmod((double)x, (double)src.x);
		if(ret.x < 0) ret.x += src.x;
		ret.y = (T)::fmod((double)y, (double)src.y);
		if(ret.y < 0) ret.y += src.y;
		ret.z = (T)::fmod((double)z, (double)src.z);
		if(ret.z < 0) ret.z += src.z;
		return ret;
	}

	clVector3D<T> mod1(const clVector3D<T>& src)
	{
		clVector3D<T> ret = clVector3D<T>(x,y,z);
		if(ret.x < 0) ret.x += src.x;
		else if(ret.x >= src.x) ret.x -= src.x;
		if(ret.y < 0) ret.y += src.y;
		else if(ret.y >= src.y) ret.y -= src.y;
		if(ret.z < 0) ret.z += src.z;
		else if(ret.z >= src.z) ret.z -= src.z;
		return ret;
	}

	clVector3D<T> max(const clVector3D<T>& src)
	{
		clVector3D<T> ret;
		ret.x = (x > src.x ? x : src.x);
		ret.y = (y > src.y ? y : src.y);
		ret.z = (z > src.z ? z : src.z);
		return ret;
	}

	clVector3D<T> min(const clVector3D<T>& src)
	{
		clVector3D<T> ret;
		ret.x = (x < src.x ? x : src.x);
		ret.y = (y < src.y ? y : src.y);
		ret.z = (z < src.z ? z : src.z);
		return ret;
	}
	clVector3D<T> sort()
	{
		T x1 = x, y1 = y, z1 = z;

		if(y < x && y < z)
		{
			x1 = y;
			y1 = x;
		}
		else if(z < x && z < y)
		{
			x1 = z;
			z1 = x;
		}

		if(y1 > z1)
			return clVector3D<T>(x1, z1, y1);

		return clVector3D<T>(x1, y1, z1);
	}
	int compare(clVector3D<T> s)
	{
		int count = 0;
		if(x == s.x || x == s.y || x == s.z)
			count += 1;
		if(y == s.x || y == s.y || y == s.z)
			count += 2;
		if(z == s.x || z == s.y || z == s.z)
			count += 4;
		return count;
	}

	T volume(clVector3D<T> v1, clVector3D<T> v2)
	{
		return fabs(*this * (v1 ^ v2) / 6.);
	}


	operator clVector3D<double> ()
	{
		return clVector3D<double>((double)x, (double)y, (double)z);
	}
	operator clVector3D<int> ()
	{
		return clVector3D<int>((int)x, (int)y, (int)z);
	}
	clVector3D<T>& operator= (const clVector3D<T>& src)
	{
		x = src.x;
		y = src.y;
		z = src.z;
		return *this;
	}
	clVector3D<T>& operator= (const double src)
	{
		x = src;
		y = src;
		z = src;
		return *this;
	}

	clVector3D<T>& operator+= (const clVector3D<T>& src)
	{
		x += src.x;
		y += src.y;
		z += src.z;
		return *this;
	}

	clVector3D<T>& operator+= (const double src)
	{
		x += src;
		y += src;
		z += src;
		return *this;
	}

	clVector3D<T>& operator-= (const int src)
	{
		x -= src;
		y -= src;
		z -= src;
		return *this;
	}

	clVector3D<T>& operator+= (const int src)
	{
		x += src;
		y += src;
		z += src;
		return *this;
	}

	clVector3D<T>& operator-= (const clVector3D<T>& src)
	{
		x -= src.x;
		y -= src.y;
		z -= src.z;
		return *this;
	}
	clVector3D<T>& operator-= (const double src)
	{
		x -= src;
		y -= src;
		z -= src;
		return *this;
	}
	clVector3D<T>& operator*= (const double src)
	{
		x *= src;
		y *= src;
		z *= src;
		return *this;
	}
	clVector3D<T>& operator/= (const double src)
	{
		x /= src;
		y /= src;
		z /= src;
		return *this;
	}
	clVector3D<T> operator^= (const clVector3D<T>& src)
	{
		T xx = y * src.z - z * src.y, yy = z * src.x - x * src.z, zz = x * src.y - y * src.x;
		x = xx;
		y = yy;
		z = zz;
		return *this;
	}

	clVector3D<T> operator%= (const clVector3D<int>& src)
	{
		x = (((x) >= (src.x)) ? ((x) - (src.x)) : (((x) < 0) ? ((x) + (src.x)) : (x)));
		y = (((y) >= (src.y)) ? ((y) - (src.y)) : (((y) < 0) ? ((y) + (src.y)) : (y)));
		z = (((z) >= (src.z)) ? ((z) - (src.z)) : (((z) < 0) ? ((z) + (src.z)) : (z)));
		return *this;
	}

	clVector3D<T> operator&= (const clVector3D<int>& src)
	{
		x = (((x) >= (src.x)) ? (src.x) : (((x) < 0) ? (0) : (x)));
		y = (((y) >= (src.y)) ? (src.y) : (((y) < 0) ? (0) : (y)));
		z = (((z) >= (src.z)) ? (src.z) : (((z) < 0) ? (0) : (z)));
		return *this;
	}

	clVector3D<T> operator+ (const clVector3D<T>& src)
	{
		return clVector3D<T>(x + src.x, y + src.y, z + src.z);
	}
	clVector3D<T> operator+ (const double src)
	{
		return clVector3D<double>(x + src, y + src, z + src);
	}
	clVector3D<T> operator+ (const int src)
	{
		return clVector3D<T>(x + src, y + src, z + src);
	}
	clVector3D<T> operator- ()
	{
		return clVector3D<T>(-x, -y, -z);
	}

	clVector3D<T> operator- (const clVector3D<T>& src)
	{
		return clVector3D<T>(x - src.x, y - src.y, z - src.z);
	}
	clVector3D<T> operator- (const double src)
	{
		return clVector3D<T>(x - src, y - src, z - src);
	}
	clVector3D<T> operator- (const int src)
	{
		return clVector3D<T>(x - src, y - src, z - src);
	}

	double operator* (const clVector3D<T>& src)
	{
		return (x * src.x + y * src.y + z * src.z);
	}
	clVector3D<T> operator* (const double src)
	{
		return clVector3D<double>(x * src, y * src, z * src);
	}
	clVector3D<T> operator* (const int src)
	{
		return clVector3D<T>(x * src, y * src, z * src);
	}
	clVector3D<T> operator| (const clVector3D<T>& src)
	{
		return clVector3D<T>(x * src.x, y * src.y, z * src.z);
	}

	clVector3D<T> operator^ (const clVector3D<T>& src)
	{
		return clVector3D<T>(y * src.z - z * src.y, 
							z * src.x - x * src.z, x * src.y - y * src.x);
	}
	clVector3D<T> operator/ (const int src)
	{
		return clVector3D<T>(x / src, y / src, z / src);
	}
	clVector3D<T> operator/ (const double src)
	{
		return clVector3D<double>(x / src, y / src, z / src);
	}
	clVector3D<T> operator% (const clVector3D<int>& src)
	{
		return clVector3D<int>((int)x % src.x, (int)y % src.y, (int)z % src.z);
	}
	clVector3D<T> operator/ (const clVector3D<T>& src)
	{
		return clVector3D<T>(x / src.x, y / src.y, z / src.z);
	}
	clVector3D<T> operator% (const int src)
	{
		return clVector3D<int>((int)x % src, (int)y % src, (int)z % src);
	}

	int operator== (const clVector3D<T>& src)
	{
		return x == src.x && y == src.y && z == src.z;
	}
	int operator== (const double src)
	{
		return x == src && y == src && z == src;
	}
	int operator!= (const clVector3D<T>& src)
	{
		return x != src.x || y != src.y || z != src.z;
	}
	int operator!= (const double src)
	{
		return x != src || y != src || z != src;
	}
	int operator >= (const clVector3D<T>& src)
	{
		return x >= src.x && y >= src.y && z >= src.z;
	}
	int operator > (const clVector3D<T>& src)
	{
		return x > src.x && y > src.y && z > src.z;
	}
	int operator <= (const clVector3D<T>& src)
	{
		return x <= src.x && y <= src.y && z <= src.z;
	}
	int operator < (const clVector3D<T>& src)
	{
		return x < src.x && y < src.y && z < src.z;
	}


};

template <class T> class clTensorS2D
{
public:
	T xx, yy, xy;

	clTensorS2D(const T xx_, const T yy_, const T xy_)
	{
		ini(xx_, yy_, xy_);
	}
	clTensorS2D()
	{
		ini();
	}
	clTensorS2D(const clVector2D<T>& src)
	{
		vect(src);
	}
	clTensorS2D(const clVector2D<T>& src1, const clVector2D<T>& src2)
	{
		vect(src1, src2);
	}
	void ini()
	{
		ini(0, 0, 0);
	}

	T& e(const int ind)
	{
		switch(ind)
		{
		case 0:
			return xx;
		case 1:
			return yy;
		case 2:
			return xy;
		}
		return xy;
	}

	void unit()
	{
		ini(1, 1, 0);
	}

	void ini(const T xx_, const T yy_, const T xy_)
	{
		xx = xx_;
		yy = yy_;
		xy = xy_;
	}
	clTensorS2D<T> diag()
	{
		return clTensorS2D<T>(xx, yy, 0);
	}
	clTensorS2D<T> trace()
	{
		return diag();
	}
	T tr_sum()
	{
		return xx + yy;
	}
	T det()
	{
		return (xx * yy - xy * xy);
	}
	clTensorS2D<T> inv()
	{
		T d = det();
		if(d != (T)0)
			return clTensorS2D<T>(yy / d, xx / d, -xy / d);
		else
			return clTensorS2D<T>(0, 0, 0);
	}
	clTensorS2D<T>& vect(const clVector2D<T>& src)
	{
		xx = src.x * src.x;
		yy = src.y * src.y;
		xy = src.x * src.y;
		return *this;
	}
	clTensorS2D<T>& vect(const clVector2D<T>& src1, const clVector2D<T>& src2)
	{
		xx = 2. * src1.x * src2.x;
		yy = 2. * src1.y * src2.y;
		xy = src1.x * src2.y + src2.x * src1.y;
		return *this;
	}


	operator clTensorS2D<double> ()
	{
		return clTensorS2D<double>((double)xx, (double)yy, (double)xy);
	}
	operator clTensorS2D<int> ()
	{
		return clTensorS2D<int>((int)xx, (int)yy, (int)xy);
	}
	clTensorS2D<T>& operator= (const clTensorS2D<T>& src)
	{
		xx = src.xx;
		yy = src.yy;
		xy = src.xy;
		return *this;
	}
	clTensorS2D<T>& operator= (const double src)
	{
		xx = yy = xy = src;
		return *this;
	}
	clTensorS2D<T>& operator= (const int src)
	{
		xx = yy = xy = src;
		return *this;
	}
	clTensorS2D<T>& operator+= (const clTensorS2D<T>& src)
	{
		xx += src.xx;
		yy += src.yy;
		xy += src.xy;
		return *this;
	}
	clTensorS2D<T>& operator+= (const double src)
	{
		xx += src;
		yy += src;
		xy += src;
		return *this;
	}
	clTensorS2D<T>& operator+= (const int src)
	{
		xx += src;
		yy += src;
		xy += src;
		return *this;
	}
	clTensorS2D<T>& operator-= (const clTensorS2D<T>& src)
	{
		xx -= src.xx;
		yy -= src.yy;
		xy -= src.xy;
		return *this;
	}
	clTensorS2D<T>& operator-= (const double src)
	{
		xx -= src;
		yy -= src;
		xy -= src;
		return *this;
	}
	clTensorS2D<T>& operator-= (const int src)
	{
		xx -= src;
		yy -= src;
		xy -= src;
		return *this;
	}
	clTensorS2D<T>& operator*= (const double src)
	{
		xx *= src;
		yy *= src;
		xy *= src;
		return *this;
	}
	clTensorS2D<T>& operator*= (const int src)
	{
		xx *= src;
		yy *= src;
		xy *= src;
		return *this;
	}
	clTensorS2D<T>& operator/= (const double src)
	{
		xx /= src;
		yy /= src;
		xy /= src;
		return *this;
	}
	clTensorS2D<T>& operator/= (const int src)
	{
		xx /= src;
		yy /= src;
		xy /= src;
		return *this;
	}

	clTensorS2D<T> operator+ (const clTensorS2D<T>& src)
	{
		return clTensorS2D<T>(xx + src.xx, yy + src.yy, xy + src.xy);
	}
	clTensorS2D<T> operator+ (const double src)
	{
		return clTensorS2D<double>(xx + src, yy + src, xy + src);
	}
	clTensorS2D<T> operator+ (const int src)
	{
		return clTensorS2D<T>(xx + src, yy + src, xy + src);
	}

	clTensorS2D<T> operator- (const clTensorS2D<T>& src)
	{
		return clTensorS2D<T>(xx - src.xx, yy - src.yy, xy - src.xy);
	}
	clTensorS2D<T> operator- (const double src)
	{
		return clTensorS2D<double>(xx - src, yy - src, xy - src);
	}
	clTensorS2D<T> operator- (const int src)
	{
		return clTensorS2D<T>(xx - src, yy - src, xy - src);
	}

	T operator* (const clTensorS2D<T>& src)
	{
		return (xx * src.xx + yy * src.yy + 2. * xy * src.xy);
	}
	clVector2D<T> operator* (const clVector2D<T>& src)
	{
		return clVector2D<T>(xx * src.x + xy * src.y, xy * src.x + yy * src.y);
	}
	clTensorS2D<T> operator^ (const clTensorS2D<T>& src)
	{
		return clTensorS2D<T>(xx * src.xx, yy * src.yy, xy * src.xy);
	}
	clTensorS2D<T> operator* (const double src)
	{
		return clTensorS2D<double>(xx * src, yy * src, xy * src);
	}
	clTensorS2D<T> operator* (const int src)
	{
		return clTensorS2D<T>(xx * src, yy * src, xy * src);
	}
	clTensorS2D<T> operator/ (const double src)
	{
		return clTensorS2D<double>(xx / src, yy / src, xy / src);
	}
	clTensorS2D<T> operator/ (const int src)
	{
		return clTensorS2D<T>(xx / src, yy / src, xy / src);
	}

	int operator== (const clTensorS2D<T>& src)
	{
		return xx == src.xx && yy == src.yy && xy == src.xy;
	}
	int operator== (const double src)
	{
		return xx == src && yy == src && xy == src;
	}
	int operator== (const int src)
	{
		return xx == src && yy == src && xy == src;
	}
	int operator!= (const clTensorS2D<T>& src)
	{
		return xx != src.xx || yy != src.yy || xy != src.xy;
	}
	int operator!= (const double src)
	{
		return xx != src || yy != src || xy != src;
	}
	int operator!= (const int src)
	{
		return xx != src || yy != src || xy != src;
	}
};
template <class T> class clTensorS3D
{
public:
	T xx, yy, zz, xy, xz, yz;

	clTensorS3D(const T xx_, const T yy_, const T zz_, const T xy_, const T xz_, const T yz_)
	{
		ini(xx_, yy_, zz_, xy_, xz_, yz_);
	}
	clTensorS3D()
	{
		ini();
	}
	clTensorS3D(const clVector3D<T>& src)
	{
		vect(src);
	}
	clTensorS3D(const clVector3D<T>& src1, const clVector3D<T>& src2)
	{
		vect(src1, src2);
	}
	void ini()
	{
		ini(0, 0, 0, 0, 0, 0);
	}

	T& e(const int ind)
	{
		switch(ind)
		{
		case 0:
			return xx;
		case 1:
			return yy;
		case 2:
			return zz;
		case 3:
			return xy;
		case 4:
			return xz;
		case 5:
			return yz;
		}
		return yz;
	}

	void unit()
	{
		ini(1, 1, 1, 0, 0, 0);
	}

	void ini(const T xx_, const T yy_, const T zz_, const T xy_, const T xz_, const T yz_)
	{
		xx = xx_;
		yy = yy_;
		zz = zz_;
		xy = xy_;
		xz = xz_;
		yz = yz_;
	}
	clTensorS3D<T> diag()
	{
		return clTensorS3D<T>(xx, yy, zz, 0, 0, 0);
	}
	clTensorS3D<T> trace()
	{
		return diag();
	}
	T tr_sum()
	{
		return xx + yy + zz;
	}
	T det()
	{
		return (xx*yy*zz + (T)2*xy*xz*yz - xz*xz*yy - xy*xy*zz - yz*yz*xx);
	}
	clTensorS3D<T> inv()
	{
		T d = det();
		if(d == (T)0)
			return clTensorS3D<T>(0, 0, 0, 0, 0, 0);

		T axx = yy * zz - yz * yz;
		T axy = -(xy * zz - yz * xz);
		T axz = xy * yz - yy * xz;
		T ayy = xx * zz - xz * xz;
		T ayz = -(xx * yz - xy * xz);
		T azz = xx * yy - xy * xy;
		return clTensorS3D<T>(axx / d, ayy / d, azz / d, axy / d, axz / d, ayz / d);
	}
	clTensorS3D<T>& vect(const clVector3D<T>& src)
	{
		xx = src.x * src.x;
		yy = src.y * src.y;
		zz = src.z * src.z;
		xy = src.x * src.y;
		xz = src.x * src.z;
		yz = src.y * src.z;
		return *this;
	}
	clTensorS3D<T>& vect(const clVector3D<T>& src1, const clVector3D<T>& src2)
	{
		xx = 2. * src1.x * src2.x;
		yy = 2. * src1.y * src2.y;
		zz = 2. * src1.z * src2.z;
		xy = src1.x * src2.y + src2.x * src1.y;
		xz = src1.x * src2.z + src2.x * src1.z;
		yz = src1.y * src2.z + src2.y * src1.z;
		return *this;
	}


	operator clTensorS3D<double> ()
	{
		return clTensorS3D<double>((double)xx, (double)yy, (double)zz, (double)xy, (double)xz, (double)yz);
	}
	operator clTensorS3D<int> ()
	{
		return clTensorS3D<int>((int)xx, (int)yy, (int)zz, (int)xy, (int)xz, (int)yz);
	}
	clTensorS3D<T>& operator= (const clTensorS3D<T>& src)
	{
		xx = src.xx;
		yy = src.yy;
		zz = src.zz;
		xy = src.xy;
		xz = src.xz;
		yz = src.yz;
		return *this;
	}
	clTensorS3D<T>& operator= (const double src)
	{
		xx = yy = zz = xy = xz = yz = (T)src;
		return *this;
	}
	clTensorS3D<T>& operator= (const int src)
	{
		xx = yy = zz = xy = xz = yz = (T)src;
		return *this;
	}
	clTensorS3D<T>& operator+= (const clTensorS3D<T>& src)
	{
		xx += src.xx;
		yy += src.yy;
		zz += src.zz;
		xy += src.xy;
		xz += src.xz;
		yz += src.yz;
		return *this;
	}
	clTensorS3D<T>& operator+= (const double src)
	{
		xx += src;
		yy += src;
		zz += src;
		xy += src;
		xz += src;
		yz += src;
		return *this;
	}
	clTensorS3D<T>& operator+= (const int src)
	{
		xx += src;
		yy += src;
		zz += src;
		xy += src;
		xz += src;
		yz += src;
		return *this;
	}
	clTensorS3D<T>& operator-= (const clTensorS3D<T>& src)
	{
		xx -= src.xx;
		yy -= src.yy;
		zz -= src.zz;
		xy -= src.xy;
		xz -= src.xz;
		yz -= src.yz;
		return *this;
	}
	clTensorS3D<T>& operator-= (const double src)
	{
		xx -= src;
		yy -= src;
		zz -= src;
		xy -= src;
		xz -= src;
		yz -= src;
		return *this;
	}
	clTensorS3D<T>& operator-= (const int src)
	{
		xx -= src;
		yy -= src;
		zz -= src;
		xy -= src;
		xz -= src;
		yz -= src;
		return *this;
	}
	clTensorS3D<T>& operator*= (const double src)
	{
		xx *= src;
		yy *= src;
		zz *= src;
		xy *= src;
		xz *= src;
		yz *= src;
		return *this;
	}
	clTensorS3D<T>& operator*= (const int src)
	{
		xx *= src;
		yy *= src;
		zz *= src;
		xy *= src;
		xz *= src;
		yz *= src;
		return *this;
	}
	clTensorS3D<T>& operator/= (const double src)
	{
		xx /= src;
		yy /= src;
		zz /= src;
		xy /= src;
		xz /= src;
		yz /= src;
		return *this;
	}
	clTensorS3D<T>& operator/= (const int src)
	{
		xx /= src;
		yy /= src;
		zz /= src;
		xy /= src;
		xz /= src;
		yz /= src;
		return *this;
	}

	clTensorS3D<T> operator+ (const clTensorS3D<T>& src)
	{
		return clTensorS3D<T>(xx + src.xx, yy + src.yy, zz + src.zz, xy + src.xy, xz + src.xz, yz + src.yz);
	}
	clTensorS3D<T> operator+ (const double src)
	{
		return clTensorS3D<double>(xx + src, yy + src, zz + src, xy + src, xz + src, yz + src);
	}
	clTensorS3D<T> operator+ (const int src)
	{
		return clTensorS3D<T>(xx + src, yy + src, zz + src, xy + src, xz + src, yz + src);
	}

	clTensorS3D<T> operator- (const clTensorS3D<T>& src)
	{
		return clTensorS3D<T>(xx - src.xx, yy - src.yy, zz - src.zz, xy - src.xy, xz - src.xz, yz - src.yz);
	}
	clTensorS3D<T> operator- (const double src)
	{
		return clTensorS3D<double>(xx - src, yy - src, zz - src, xy - src, xz - src, yz - src);
	}
	clTensorS3D<T> operator- (const int src)
	{
		return clTensorS3D<T>(xx - src, yy - src, zz - src, xy - src, xz - src, yz - src);
	}

	T operator* (const clTensorS3D<T>& src)
	{
		return (xx * src.xx + yy * src.yy + zz * src.zz + 2. * xy * src.xy + 2. * xz * src.xz + 2. * yz * src.yz);
	}
	clVector3D<T> operator* (const clVector3D<T>& src)
	{
		return clVector3D<T>(xx * src.x + xy * src.y + xz * src.z, xy * src.x + yy * src.y + yz * src.z, xz * src.x + yz * src.y + zz * src.z);
	}
	clTensorS3D<T> operator^ (const clTensorS3D<T>& src)
	{
		return clTensorS3D<T>(xx * src.xx, yy * src.yy, zz * src.zz, xy * src.xy, xz * src.xz, yz * src.yz);
	}
	clTensorS3D<T> operator* (const double src)
	{
		return clTensorS3D<double>(xx * src, yy * src, zz * src, xy * src, xz * src, yz * src);
	}
	clTensorS3D<T> operator* (const int src)
	{
		return clTensorS3D<T>(xx * src, yy * src, zz * src, xy * src, xz * src, yz * src);
	}
	clTensorS3D<T> operator/ (const double src)
	{
		return clTensorS3D<double>(xx / src, yy / src, zz / src, xy / src, xz / src, yz / src);
	}
	clTensorS3D<T> operator/ (const int src)
	{
		return clTensorS3D<T>(xx / src, yy / src, zz / src, xy / src, xz / src, yz / src);
	}

	int operator== (const clTensorS3D<T>& src)
	{
		return xx == src.xx && yy == src.yy && zz == src.zz && xy == src.xy && xz == src.xz && yz == src.yz;
	}
	int operator== (const double src)
	{
		return xx == src && yy == src && zz == src && xy == src && xz == src && yz == src;
	}
	int operator== (const int src)
	{
		return xx == src && yy == src && zz == src && xy == src && xz == src && yz == src;
	}
	int operator!= (const clTensorS3D<T>& src)
	{
		return xx != src.xx || yy != src.yy || zz != src.zz || xy != src.xy || xz != src.xz || yz != src.yz;
	}
	int operator!= (const double src)
	{
		return xx != src || yy != src || zz != src || xy != src || xz != src || yz != src;
	}
	int operator!= (const int src)
	{
		return xx != src || yy != src || zz != src || xy != src || xz != src || yz != src;
	}
};
typedef  clVector2D <double> clVector2Dd;
typedef  clVector2D <int> clVector2Di;
typedef  clVector2D <char> clVector2Dc;
typedef  clVector3D <double> clVector3Dd;
typedef  clVector3D <int> clVector3Di;
typedef  clVector3D <char> clVector3Dc;

typedef  clTensorS2D <double> clTensorS2Dd;
typedef  clTensorS2D <int> clTensorS2Di;
typedef  clTensorS2D <char> clTensorS2Dc;
typedef  clTensorS3D <double> clTensorS3Dd;
typedef  clTensorS3D <int> clTensorS3Di;
typedef  clTensorS3D <char> clTensorS3Dc;

#endif 
