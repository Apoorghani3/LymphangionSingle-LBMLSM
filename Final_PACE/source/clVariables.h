// clVariables.h: interface for the clVariables class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLVARIABLES_H__INCLUDED_)
#define AFX_CLVARIABLES_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
//#include "clVariable.h"
#include "clProblem.h"

class clVariables  
{
public:
	clVariables(){};
	virtual ~clVariables(){};

	void ini(void);
	void update(void);
	void error(char *);

};
/////////////////////////////////////////////////////////////////////////
class clVariablesLS
{
public:
	clVariablesLS(){};
	virtual ~clVariablesLS(){};
	class clForce
	{
	public:
		double f, p;
		clForce()
		{
			//f = 0.;
			//p = 0.;
			ini(0., 0.);
		};
		clForce(double f0, double p0)
		{
			ini(f0, p0);
		}
		void ini(double f0, double p0)
		{
			f = f0;
			p = p0;
		};
		clForce operator* (const double src){return clForce(f * src, p * src);};
		clForce operator+ (const clForce src){return clForce(f + src.f, p + src.p);};
		clForce operator- (const clForce src){return clForce(f - src.f, p - src.p);};
	};
	class clForceV
	{
	public:
		double p;
		clVector3Dd f;
		clForceV()
		{
			ini(clVector3Dd(0., 0., 0.), 0.);
		};
		clForceV(clVector3Dd f0, double p0)
		{
			ini(f0, p0);
		};
		void ini(clVector3Dd f0, double p0)
		{
			f = f0;
			p = p0;
		};
		clForceV operator* (const double src){return clForceV(f * src, p * src);};
		clForceV operator+ (const clForceV src){return clForceV(f + src.f, p + src.p);};
		clForceV operator- (const clForceV src){return clForceV(f - src.f, p - src.p);};
	};

	void ini()
	{
		C.ini();
		V.ini();
		F.ini();

		c.ls.L0.ini();
		c.ls.iniAA0();
		c.ls.iniP();

		//		Ta.ini();
		//Tp.ini();
		P.ini();
		D.ini();

		R.ini();
		avr.ini();

		//		I.update();

		//C.ini2();
		//		V.ini2();
		//		F.ini2();

	};

	void update()
	{
		C.update();
		V.update();
		F.update();
		//		Ta.update();

		D.update();

		P.update();
		//Tp.update();

	};


	void save2file()
	{
#ifdef LS_SOLVER
		C.save2file();
		F.save2file();
		V.save2file();
		P.save2file();
		//Tp.save2file();
		D.save2file();
		c.ls.BLFb.savetxt("lsblfb");
#endif //LS_SOLVER
	}

	class clD
	{
	public:
		clArray2Dd Symax, Symin, Scur, Sm1;
		clArray1Dd S;
		clArray1Di Sn;
		double Smax;
		int step_counter;
		void update()
		{
#ifndef LS_RUPTURE
			//double Symax = c.ls.Symax0, Symin = c.ls.Symin0;
			//if(Symax == 0.)	Symax = 1.;
			//if(Symin == 0.)	Symin = 1.;
			//double smax = vls.C.Smax/Symax, smin = -vls.C.Smin/Symin;
			double smax = vls.C.Smax, smin = -vls.C.Smin;
			if(smax > smin)
				Smax = smax;
			else
				Smax = smin;
			return;
#endif //LS_RUPTURE

			//if(step_counter++ < 5)
			//	return;
			//step_counter = 0;

			Sm1.copy(Scur);
			Scur.zeros();
			S.zeros();
			double s_cur;
			double f1 = LS_CRIT_F1 * c.t.lsSteps;
			int nmax0, nmax1, jmax;
			Smax = 0.;
			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n0 = c.ls.BN(i);
				if(c.ls.DO(c.ls.O(n0)) == FALSE || c.ls.N.NL(n0) == 0)
					continue;

				for(int j = 0; j < c.ls.nL; j++)
				{
					int n1 = c.ls.N(n0,j);
					if(n1 < 0)
						continue;
					s_cur = strain(n0, j);
					//Scur(n0,j) = s_cur;
					//s_cur = (s_cur + f1 * Sm1(n0,j)) / (1 + f1);
					Scur(n0,j) = (s_cur + f1 * Sm1(n0,j)) / (1. + f1);
					if(s_cur > Smax)
					{
						nmax0 = n0;
						nmax1 = n1;
						jmax = j;
						Smax = s_cur;
					}
					if(s_cur > S(n0))
					{
						S(n0) = s_cur;
						Sn(n0) = n1;
					}
				}

			}
			double s_max;
			BOOL flCorrect = FALSE;
			while(1)
			{
				s_max = 0.;
				for(int i = 0; i < c.ls.nBN; i++)
				{
					int n0 = c.ls.BN(i);
					if(S(n0) > s_max && S(Sn(n0)) == S(n0))
					{
						s_max = S(n0);
						nmax0 = n0;
						nmax1 = Sn(n0);
					}
				}
				if(s_max > 1.)
				{
					for(int j = 0; j < c.ls.nL; j++)
						if(c.ls.N(nmax0,j) >= 0)
							S(c.ls.N(nmax0,j)) = -1.;

					for(int j = 0; j < c.ls.nL; j++)
						if(c.ls.N(nmax1,j) >= 0)
							S(c.ls.N(nmax1,j)) = -1.;

					S(nmax0) = -1.;
					S(nmax1) = -1.;
					remove_link(nmax0, nmax1);
					flCorrect = TRUE;
					//					break;
				}
				else
					break;
			}

			if(flCorrect == TRUE)
				correct_tables();

			//if(e_max_tot > 1.)
			//	remove_link(nmax0, nmax1);
			//return;

			//if(e_max_tot < 1.)
			//	return;
			//if(jmax == 6 || jmax == 13)
			//{
			//	remove_link(nmax0, nmax1);
			//	return;
			//}
			//int nv0, nv1;
			//if(c.ls.N(nmax0,6) >= 0)
			//	nv0 = c.ls.N(nmax0,6);
			//else
			//	nv0 = c.ls.N(nmax0,13);

			//if(c.ls.N(nmax1,6) >= 0)
			//	nv1 = c.ls.N(nmax1,6);
			//else
			//	nv1 = c.ls.N(nmax1,13);

			//remove_link(nmax0, nmax1);
			//remove_link(nv0, nv1);
			//remove_link(nmax0, nv1);
			//remove_link(nmax1, nv0);

			//for(int j = 0; j < c.ls.nL; j++)
			//{
			//	int n = c.ls.N(nmax0,j);
			//	if(n < 0)
			//		continue;
			//	remove_link(nmax0, n);
			//}
		};
		double strain(int n0, int j)
		{
			double S = vls.C.S(n0,j);

			if(S > 0.)
				return S / Symax(n0,j);
			if(S < 0.)
				return S / Symin(n0,j);

			return 0.;
		};
		void remove_link(int n1, int n2);
		void correct_tables();
		void save2file(){};
		void ini()
		{
#ifndef LS_RUPTURE
			return;
#endif //LS_RUPTURE
			step_counter = 0;
			Scur.zeros(c.ls.nN, c.ls.nL);
			Sm1.zeros(c.ls.nN, c.ls.nL);
			S.zeros(c.ls.nN);
			Sn.zeros(c.ls.nN);

			Symax.zeros(c.ls.nN, c.ls.nL);
			Symin.zeros(c.ls.nN, c.ls.nL);
			if(Symax.load("load" SLASH "lsdsymax") == FALSE)
			{
				for(int i = 0; i < c.ls.nN; i++)
				{
					if(c.ls.DO(c.ls.O(i)) == FALSE)
						continue;
					for(int j = 0; j < c.ls.nL; j++)
						if(c.ls.N(i,j) >= 0)
						{
							Symax(i,j) = c.ls.Symax0;
							if(c.ls.O(i) == 1)
								Symax(i,j) = c.ls.Symax1;
							//Modified
							if(c.ls.O(i) >= 2)
							//Modified
								Symax(i,j) = c.ls.Symax2;
						}
				}
				setnoise_max();
			}
			if(Symin.load("load" SLASH "lsdsymin") == FALSE)
			{
				for(int i = 0; i < c.ls.nN; i++)
				{
					if(c.ls.DO(c.ls.O(i)) == FALSE)
						continue;
					for(int j = 0; j < c.ls.nL; j++)
						if(c.ls.N(i,j) >= 0)
						{
							Symin(i,j) = c.ls.Symin0;
							if(c.ls.O(i) == 1)
								Symin(i,j) = c.ls.Symin1;
							//Modified
							if(c.ls.O(i) >= 2)
							//Modified
								Symin(i,j) = c.ls.Symin2;
						}
				}
				setnoise_min();
			}
			for(int i = 0; i < c.ls.nN; i++)
			{
				for(int j = 0; j < c.ls.nL; j++)
				{
					if(c.ls.N(i,j) < 0)
						continue;
					if(Symax(i,j) == 0.)
						Symax(i,j) = 1.;
					if(Symin(i,j) == 0.)
						Symin(i,j) = -1.;
				}
			}
			Symax.savetxt("lsdsymax");
			Symin.savetxt("lsdsymin");
		};
		void setnoise_max()
		{
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.DO(c.ls.O(i)) == FALSE)
					continue;
				for(int j = 0; j < c.ls.nL; j++)
					if(c.ls.N(i,j) >= 0 && i > c.ls.N(i,j))
					{
						Symax(i,j) = Symax(i,j) * (1. + LS_SY_NOISE * (c.t.random_normal() - 0.5));
						Symax(c.ls.N(i,j),c.ls.N.find(c.ls.N(i,j),i)) = Symax(i,j);
					}
			}
		};
		void setnoise_min()
		{
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.DO(c.ls.O(i)) == FALSE)
					continue;
				for(int j = 0; j < c.ls.nL; j++)
					if(c.ls.N(i,j) >= 0 && i > c.ls.N(i,j))
					{
						Symin(i,j) = Symin(i,j) * (1. + LS_SY_NOISE * (c.t.random_normal() - 0.5));
						Symin(c.ls.N(i,j),c.ls.N.find(c.ls.N(i,j),i)) = Symin(i,j);
					}
			}
		};
	} D;

	class clC: public clArray1Dv3d
	{
	public:
		clArray1Dv3d c0;
		clArray1Dv3i p;
		clArray1Dv3d d;
		clArray2Dd L, S;
		clVector3Dd Cmin,Cmax;		

		double Smax, Smin;
		int MinMaxBool;
		void ini()
		{
			zeros(c.ls.nN);
			c0.zeros(c.ls.nN);
			p.zeros(c.ls.nN);
			d.zeros(c.ls.nN);
			L.zeros(c.ls.nN, c.ls.nL);
			S.zeros(c.ls.nN, c.ls.nL);
			if(load("load" SLASH "lsc") == FALSE)
			{
				if(loadtxt_xyz("load" SLASH "lsc") == FALSE)
				{
					//v.error("Cannot initialze Cx in LS");
				}
			}

			if(c.t.flContinue == FALSE)
			{
#ifdef WING_PITCHING
				c0.copy(*this);
				double sinamp = -OSC_AMPL / WING_LENGTH;
				double cosamp = sqrt(1. - sinamp * sinamp);
				clVector3Dd Czero1 = el(c.ls.nN-1) - clVector3Dd(-0.5, 0., 0.), 
							Czero0 = el(0) - clVector3Dd(-0.5, 0., 0.), dCzero = Czero1 - Czero0;
				double Yzero = Czero1.y - dCzero.y / dCzero.x * Czero1.x;
				for(int i = 0; i < c.ls.nN; i++)
				{
					clVector3Dd C0 = el(i) - clVector3Dd(-0.5, Yzero, 0.), C1;
					C1.x = cosamp * C0.x + sinamp * C0.y;
					C1.y = -sinamp * C0.x + cosamp * C0.y;
					C1.z = C0.z;
					clVector3Dd dC = C1 - C0;
					el(i) += dC;
				}
#else
				for(int i = 0; i < c.ls.nN; i++)
				{
					if(c.ls.O(i) == 0)
						el(i) += c.ls.Cini0;
					if(c.ls.O(i) == 1)
						el(i) += c.ls.Cini1;
					//Modified
					if(c.ls.O(i) >= 2)
					//Modified
						el(i) += c.ls.Cini2;
				}
				c0.copy(*this);
#endif //WING_ROTATION
			}
			else
			{
				if(c0.loadtxt("load" SLASH "lsc") == TRUE)
				{
					for(int i = 0; i < c.ls.nN; i++)
					{
						if(c.ls.O(i) == 0)
							c0(i) += c.ls.Cini0;
						if(c.ls.O(i) == 1)
							c0(i) += c.ls.Cini1;
						//Modified
						if(c.ls.O(i) >= 2)
						//Modified
							c0(i) += c.ls.Cini2;
					}
				}
				else
					c0.copy(*this);
			}

			// Valve extent: computed once here from the loaded/offset geometry.
			// valveMinX = lowest X among nodes of the first valve object (lso == 1).
			// valveMaxX = lowest X among nodes of the last valve object (lso == max(lso)).
			{
				int maxO = 0;
				for(int i = 0; i < c.ls.nN; i++)
					if(c.ls.O(i) > maxO) maxO = c.ls.O(i);

				double minXFirst = 1.e300, minXLast = 1.e300;
				for(int i = 0; i < c.ls.nN; i++)
				{
					if(c.ls.O(i) == 1 && el(i).x < minXFirst)
						minXFirst = el(i).x;
					if(c.ls.O(i) == maxO && el(i).x < minXLast)
						minXLast = el(i).x;
				}

				if(minXFirst < 1.e300) c.ls.valveMinX = minXFirst;
				if(minXLast < 1.e300) c.ls.valveMaxX = minXLast;

				// Peristaltic wavelength: 2 x the intervalve (lymphangion) distance.
				// With the current 2-valve layout, valveMaxX-valveMinX IS the intervalve
				// distance (the span between the only two valves).
				c.ls.waveLength = 2.0 * (c.ls.valveMaxX - c.ls.valveMinX);

				printf("Valve extent: X in [%.2f, %.2f] (lso=1 min-X to lso=%d min-X), waveLength=%.2f\n",
					c.ls.valveMinX, c.ls.valveMaxX, maxO, c.ls.waveLength);
			}

			iniL();

			shift();
			save2file();
		}

		void FindMinMax()
		{
			double xmin = 100., ymin = 100., zmin = 100.;
			double xmax = 0., ymax = 0., zmax = 0.;
			MinMaxBool = 0;
			double xbound = (double)c.lb.nX - 0.5;
			for(int i = 0; i < c.ls.nN; i++)
			{
				clVector3Dd Pos = vls.C(i);
				//if(Pos.x > xmax)
				//	xmax = Pos.x;
				//if(Pos.x < xmin)
				//	xmin = Pos.x;
				//if(Pos.y > ymax)
				//	ymax = Pos.y;
				//if(Pos.y < ymin)
				//	ymin = Pos.y;
				if(Pos.z > zmax)
					zmax = Pos.z;
				if(Pos.z < zmin)
					zmin = Pos.z;
			}

			//if(xmax > xbound)
			//{
			//	xmax -= (double)c.lb.nX;
			//	MinMaxBool = 1;
			//}

			//if(xmin < -0.5)
			//{
			//	xmin += (double)c.lb.nX;
			//	MinMaxBool = 1;
			//}

			Cmin = clVector3Dd(xmin,ymin,zmin);
			Cmax = clVector3Dd(xmax,ymax,zmax);
		}

		void zerosLS()
		{
			L.zeros();
			S.zeros();
		}
		void iniL()
		{
			for(int i = 0; i < c.ls.nN; i++)
			{
				for(int j = 0; j < c.ls.nL; j++)
				{
					int n = c.ls.N(i,j);
					if(n < 0)
						continue;
					clVector3Dd d = vls.C(i) - vls.C(n);
					L(i,j) = d.abs();
				}
			}
		}
		void shift()
		{
			for(int i = 0; i < c.ls.nN; i++)
			{
				shiftx(i);
				shifty(i);
				shiftz(i);
			}
		}
		void shiftx(int i)
		{
			if(c.ls.PO(c.ls.O(i)).x == FALSE)
				return;

			double xi = el(i).x + 0.5;
			int cxi = (int)xi;
			int pxi = (cxi - cxi % c.nX) / c.nX;
			if(xi < 0.) pxi--;
			p(i).x = pxi;
		}

		void shifty(int i)
		{
			if(c.ls.PO(c.ls.O(i)).y == FALSE)
				return;

			double yi = el(i).y + 0.5;
			int cyi = (int)yi;
			int pyi = (cyi - cyi % c.nY) / c.nY;
			if(yi < 0) pyi--;
			p(i).y = pyi;
		}

		void shiftz(int i)
		{
			if(c.ls.PO(c.ls.O(i)).z == FALSE)
				return;

			double zi = el(i).z + 0.5;
			int czi = (int)zi;
			int pzi = (czi - czi % c.nZ) / c.nZ;
			if(zi < 0) pzi--;
			p(i).z = pzi;
		}

		clVector3Dd s(int n)
		{
			return el(n) - (clVector3Dd)(p(n) | c.nn);
			//return clVector3Dd(sx(n), sy(n), sz(n));
		}

		clVector3Dd s(int n, int np)
		{
			return el(n) - (clVector3Dd)(p(np) | c.nn);
			//return clVector3Dd(sx(n,np), sy(n,np), sz(n,np));
		}

		//void shiftx(int i)
		//{
		//	if(c.ls.PO(c.ls.O(i)).x == FALSE)
		//		return;

		//	double xi = el(i).x;
		//	int cxi = (int)(xi + 0.5);
		//	int pxi = (cxi - cxi % c.nX) / c.nX;
		//	if(cxi < 0) pxi--;
		//	p(i).x = pxi;
		//}

		//void shifty(int i)
		//{
		//	if(c.ls.PO(c.ls.O(i)).y == FALSE)
		//		return;

		//	double yi = el(i).y;
		//	int cyi = (int)(yi + 0.5);
		//	int pyi = (cyi - cyi % c.nY) / c.nY;
		//	if(cyi < 0) pyi--;
		//	p(i).y = pyi;
		//}

		//void shiftz(int i)
		//{
		//	if(c.ls.PO(c.ls.O(i)).z == FALSE)
		//		return;

		//	double zi = el(i).z;
		//	int czi = (int)(zi + 0.5);
		//	int pzi = (czi - czi % c.nZ) / c.nZ;
		//	if(czi < 0) pzi--;
		//	p(i).z = pzi;
		//}


		double sx(int i)
		{
			return sx(i, i);
		}

		double sy(int i)
		{
			return sy(i, i);
		}

		double sz(int i)
		{
			return sz(i, i);
		}

		double sx(int ix, int ip)
		{
			return el(ix).x - (double)(p(ip).x * c.nX);
		}

		double sy(int iy, int ip)
		{
			return el(iy).y - (double)(p(ip).y * c.nY);
		}

		double sz(int iz, int ip)
		{
			return el(iz).z - (double)(p(ip).z * c.nZ);
		}

		void ini2()
		{
			for(int i = 0; i < c.ls.nN; i++)
			{
				el(i) = (el(i) - vls.avr.CG(0)) * 1.01 + vls.avr.CG(0);
				//el(i).y *= 1.01;
				//el(i).z *= 1.01;
			}
		}

		void update()
		{
		}

		void save2file()
		{
#ifdef LS_DEBUG
			savetxt_xyz("lsc");
			p.savetxt("lscp");
#endif //LS_DEBUG
			save("lsc");
			L.savetxt("lscl");
			S.savetxt("lscs");
		};
	} C;

	class clV: public clArray1Dv3d
	{
	public:
		double max;
		clArray1Dv3d avr;
		void ini()
		{
			max = 0.;

			zeros(c.ls.nN);
			avr.zeros(c.ls.nN);
			if(load("load" SLASH "lsv"))
				if(loadtxt_xyz("load" SLASH "lsv") == FALSE)
				{
					//						v.error("Cannot initialze Cx in LS");
				}

#ifdef INI_SHEAR
				inishear();
#endif //INI_SHEAR
#ifdef INI_POISEUILLE
				//iniPoiseuille();		// Comment this out if don't want LS to be initialized with Poiseuille
#endif //INI_POISEUILLE

				avr.copy(*this);
				save2file();
				//			for(int i = 0; i < c.ls.nN; i++)
				//			{
				////				el(i).y = LS_INI_VELOCITY;
				//			}
		}
		void ini2()
		{
			//clVector3Dd Omega = clVector3Dd(0., 0., 0.001);
			for(int i = 0; i < c.ls.nN; i++)
			{
				//clVector3Dd R = vls.C(i) - vls.avr.CG(0);
				//R.z = 0.;
				//double R2 = R.r2();
				//clVector3Dd V = Omega ^ R;
				//clVector3Dd	OmegaAbs = (R ^ V) / R.r();
				//double OmegaMag = OmegaAbs.r();
				//clVector3Dd	OmegaNew = (R ^ V) / R2;
				//el(i) = V;
				//el(i).y = -(vls.C(i).x - vls.avr.CG(0).x) * 0.0001;
				//el(i).x = (vls.C(i).y - vls.avr.CG(0).y) * 0.0001;
				//el(i) = clVector3Dd(0.01, 0., 0.);
				//el(i) = (vls.C(i) - vls.avr.CG(0)) * 0.0001;
				el(i) = (vls.C(i) - vls.avr.CG(0)) * LS_INI_VELOCITY;
				//el(i).x += LS_INI_VELOCITY;
			}
			avr.copy(*this);
		}
		void inishear()
		{
			double Cyavr = 0., count = 0.;
			for(int i = 0; i < c.ls.nN; i++)
			{
				Cyavr += vls.C(i).y;
				count++;
			}
			double Ux = c.lb.UY0.x + (c.lb.UY1.x - c.lb.UY0.x) * (Cyavr / count + 0.5) / (double)c.lb.nY;
			double Uy = 0.;

			for(int i = 0; i < c.ls.nN; i++)
			{
				//Ux = c.lb.UX0x + (c.lb.UX1x - c.lb.UX0x) * (vls.C.y(i) + 0.5) / (double)c.lb.nY;
				el(i).x = Ux;
				el(i).y = Uy;
			}
		}

		void iniPoiseuille()
		{
			double Cyavr = 0., count = 0.;
			for(int i = 0; i < c.ls.nN; i++)
			{
				Cyavr += vls.C(i).y;
				count++;
			}
			double a = 0.5 * (double)c.lb.nY;
			double yy = (Cyavr / count + 0.5) - a;
			double Ux = c.lb.F(0).x * (a * a - yy * yy) / 2. * 6.;
			double Uy = 0.;

			for(int i = 0; i < c.ls.nN; i++)
			{
				//yy = (vls.C.y(i) + 0.5) - a;
				//Ux = c.lb.Fx[0] * (a * a - yy * yy) / 2. * 6.;
				el(i).x = Ux;
				el(i).y = Uy;
			}
		}

		void update()
		{
		}

		void save2file()
		{
			save("lsv");
#ifdef LS_DEBUG
			savetxt_xyz("lsv");
#endif //LS_DEBUG
		};
	} V;

	class clF: public clArray1Dv3d
	{
		//		clArray1Dv3d bs;
	public:
		clArray1Dv3d b, c, e, s, B;	//B -- magnetic field
		clArray1Dv3d ba, ca, ea, sa;
		clArray1Dd p;
		clArray1Dd Wb, Wc, We;
		void ini()
		{
			clArray1Dv3d::zeros(::c.ls.nN);
			b.zeros(::c.ls.nN);
			c.zeros(::c.ls.nN);
			e.zeros(::c.ls.nN);
			s.zeros(::c.ls.nN);
			ba.zeros(::c.ls.nN);
			ca.zeros(::c.ls.nN);
			ea.zeros(::c.ls.nN);
			sa.zeros(::c.ls.nN);
//			b_m1.zeros(::c.ls.nN);
//			b_old.zeros(::c.ls.nN);
			p.zeros(::c.ls.nN);
			Wb.zeros(::c.ls.nN);
			Wc.zeros(::c.ls.nN);
			We.zeros(::c.ls.nN);
			//			bs.zeros(::c.ls.nN);

			if(load("load" SLASH "lsf") == FALSE)
				if(loadtxt_xyz("load" SLASH "lsf") == FALSE)
				{
					//						v.error("Cannot initialze F in LS");
				}

				save2file();
			Bini();
		}
		//void ini2()
		//{
		//	for(int i = 0; i < ::c.ls.nN; i++)
		//		b(i) = (vls.C(i) - vls.avr.CG(0)) * LS_INI_FORCE;
		//}

		void update()
		{
		}
		void Bini()
		{
		}
		void updateE()
		{
			//double osc = ::c.t.oscillator1();
			//e.zeros();
			//clVector3Dd Fe = clVector3Dd(::c.t.oscillator1sin(), ::c.t.oscillator1cos(), -1.);
			//for(int i = 0; i < ::c.ls.nN; i++)
			//{
			//	if(::c.ls.O(i) < 13)
			//		e(i) = Fe / N0_NUMBER * PILLAR2CAPSULE;
			//}
			//return;

			e.zeros();
			double dx = 0.5;
			double dy = 0.5;
			double dz = 0.5;
			int dn;
			const int nd = N_DISK;			//number of disc
			const int nrd = N_ROWS_DISK;	// number of rows of discs
			const int ncols = (int)(N_DISK/N_ROWS_DISK);
			clArray1Dv3d Cdisc;			//coordinates of discs
			Cdisc.zeros(nd);
			//for(int i = 0; i < nd; i++)
			//{
			//	Cdisc(i).x = 14.5;
			//	Cdisc(i).y = 17.5 + 30. * i;

			//	Cdisc(i).z = 0;
			//}
			for(int i = 0; i < nrd; i++)
			{
				Cdisc(i).x = CX0_DISK;
				Cdisc(i).y = CY0_DISK + SPACING * i;

				Cdisc(i).z = 0.5 * (PILLAR_H1 + PILLAR_H2);
			}
			if(ncols == 2){
				for(int i = nrd; i < nd; i++)
				{
					Cdisc(i).x = CX0_DISK + SPACING;
					Cdisc(i).y = CY0_DISK - 0.5 * SPACING + SPACING * (i - nrd);

					Cdisc(i).z = 0.5 * (PILLAR_H1 + PILLAR_H2);
				}
			}
			else if(ncols > 2){
				dn = nrd-1;
				for(int id = 1; id < ncols; id++)
				{
					for(int i = 0; i < nrd; i++)
					{
						dn = dn + 1;
						Cdisc(dn).x = CX0_DISK + SPACING * id;
						Cdisc(dn).y = CY0_DISK + SPACING * i;

						Cdisc(dn).z = 0.5 * (PILLAR_H1 + PILLAR_H2);
					}
				}
			}

			//need superimposition to calculate B
			const int ndisc = nd * 9;	//# of disc, including superimposed discs
			clArray1Dv3d Cdisc1;
			Cdisc1.zeros(ndisc);
			//Cdisc1.copy(Cdisc);
			int n = 0;
			for(int i = 0; i < nd; i++)
				for(int j = 1; j < 4; j++)
					for(int k = 1; k < 4; k++)
					{
						Cdisc1(n) = Cdisc(i) + clVector3Dd(::c.nX, 0., 0.) * (j-2) + clVector3Dd(0., ::c.nY, 0.) * (k-2);
						n++;
					}
			//Cdisc1.save2file("cdisc");
			vls.avr.calcCG();
			for(int i = 0; i < N_BEAD; i++)
			{
				//clVector3Dd Ccap = clVector3Dd(vls.avr.CG(i) + clVector3Dd(0.5, 0.5, 0.5)).fmod(::c.lb.nn) - clVector3Dd(0.5, 0.5, 0.5); //Coords of beads
				//clVector3Dd Ccap = vls.C.s(i);
				//Ccap = transformR(clVector3Dd(Ccap.x - Cdisc(0).x, Ccap.y - Cdisc(0).y, Ccap.z)) + Cdisc(0);

				//fix vertical position
				//if (vls.avr.CG(i).z < 1.003 * (R0_NUMBER - 0.5))
				//{
				//	double dCGz = vls.avr.CG(i).z - 1.001 * (R0_NUMBER - 0.5);
				//	for(int j = 0; j < ::c.ls.nN; j++)
				//		if(::c.ls.O(j) == i)
				//			vls.C(j).z -= dCGz;
				//}
				clVector3Dd Ccap = clVector3Dd(vls.avr.CG(i));
				clVector3Dd B;	//magnetic field
				B = calculateB(Ccap, Cdisc1);
				B = transform(B) + clVector3Dd(B_BIAS * ::c.t.oscillator1cos(), B_BIAS * ::c.t.oscillator1sin(), 0.); //total B
				//clVector3Dd Ctemp = Ccap + clVector3Dd(::c.t.oscillator1cos()*dx, -::c.t.oscillator1sin()*dx, 0.);
				clVector3Dd Ctemp = Ccap + clVector3Dd(dx, 0., 0.);
				clVector3Dd B_dx = calculateB(Ctemp, Cdisc1);
				B_dx = transform(B_dx);
				//Ctemp = Ccap + clVector3Dd(-::c.t.oscillator1cos()*dx, ::c.t.oscillator1sin()*dx, 0.);
				Ctemp = Ccap + clVector3Dd(-dx, 0., 0.);
				clVector3Dd B_dx_ = calculateB(Ctemp, Cdisc1);
				B_dx_ = transform(B_dx_);
				// gradient in the x direction
				clVector3Dd G_Bx = (B_dx - B_dx_) / (2.*dx);

				//Ctemp = Ccap + clVector3Dd(::c.t.oscillator1sin()*dy, ::c.t.oscillator1cos()*dy, 0.);
				Ctemp = Ccap + clVector3Dd(0., dy, 0.);
				clVector3Dd B_dy = calculateB(Ctemp, Cdisc1);
				B_dy = transform(B_dy);
				//Ctemp = Ccap + clVector3Dd(-::c.t.oscillator1sin()*dy, -::c.t.oscillator1cos()*dy, 0.);
				Ctemp = Ccap + clVector3Dd(0., -dy, 0.);
				clVector3Dd B_dy_ = calculateB(Ctemp, Cdisc1);
				B_dy_ = transform(B_dy_);
				// gradient in the y direction
				clVector3Dd G_By = (B_dy - B_dy_) / (2.*dy);

				clVector3Dd B_dz = calculateB(Ccap + clVector3Dd(0., 0., dz), Cdisc1);
				clVector3Dd B_dz_ = calculateB(Ccap + clVector3Dd(0., 0., -dz), Cdisc1);
				B_dz = transform(B_dz);
				B_dz_ = transform(B_dz_);
				// gradient in the z direction
				clVector3Dd G_Bz = (B_dz - B_dz_) / (2.*dz);
				//magnetic force
				clVector3Dd Fp;
				Fp.x = B.x*G_Bx.x + B.y*G_By.x +B.z*G_Bz.x;
				Fp.y = B.x*G_Bx.y + B.y*G_By.y +B.z*G_Bz.y;
				Fp.z = B.x*G_Bx.z + B.y*G_By.z +B.z*G_Bz.z;

#ifdef BEADS_MIXING
				for(int j = 0; j < ::c.ls.nN; j++)
				{
					if(::c.ls.O(j) == i)
						e(j) = Fp / N0_NUMBER * MAGNETIC_FORCE;
				}				
#else
				if( B.abs() * MU0_MS > 0.1)
					Fp = Fp * MSP * MU0_MS / B.abs();
				else if( B.abs() * MU0_MS < 0.01)
					Fp = Fp * MU0_MS * MU0_MS * 0.54 / (4 * ::c.Pi * 1e-7);
				else
					Fp = Fp * MU0_MS * MU0_MS * 0.14 / (4 * ::c.Pi * 1e-7);
				Fp = Fp * 4./3. * ::c.Pi * R0_NUMBER * R0_NUMBER * R0_NUMBER;

				for(int j = 0; j < ::c.ls.nN; j++)
				{
					if(::c.ls.O(j) == i)
						e(j) = Fp / N0_NUMBER * FORCE_CONVERSION;
				}
#endif // BEADS_MIXING

				//if (vls.avr.CG(i).z < 1.003 * (R0_NUMBER - 0.5))
				//{
				//	for(int j = 0; j < ::c.ls.nN; j++)
				//		if(::c.ls.O(j) == i)
				//			e(j).z = 0;
				//}
				//e(i) = Fp / N0_NUMBER * MAGNETIC_FORCE;
			}

		}

		clVector3Dd transform(clVector3Dd Ccap)
		{
			clVector3Dd Ccap1;
			Ccap1.x = Ccap.x * ::c.t.oscillator1cos() - Ccap.y * ::c.t.oscillator1sin();
			Ccap1.y = Ccap.x * ::c.t.oscillator1sin() + Ccap.y * ::c.t.oscillator1cos();
			Ccap1.z = Ccap.z;
			return Ccap1;
		}

		clVector3Dd transformR(clVector3Dd B)	//reverse transformation
		{
			clVector3Dd B1;
			B1.x = B.x * ::c.t.oscillator1cos() + B.y * ::c.t.oscillator1sin();
			B1.y = -B.x * ::c.t.oscillator1sin() + B.y * ::c.t.oscillator1cos();
			B1.z = B.z;
			return B1;
		}

		clVector3Dd calculateB(clVector3Dd Ccap1, clArray1Dv3d &Cdisc)	//dimensionless B, B_real= B * mu0 * Ms
		{			
			clVector3Dd B;
			clVector3Dd Ccap;
			for(int i = 0; i < Cdisc.getFullSize(); i++)
			{
				//double dist = (Ccap1 - Cdisc(i)).abs();
				//if( dist > 80 )
				//	continue;	//magnetic force can be neglected when distant is too large
				Ccap = transformR(clVector3Dd(Ccap1.x - Cdisc(i).x, Ccap1.y - Cdisc(i).y, Ccap1.z));	//transform to polar coords wrt. disc
				//Ccap = (Ccap1 - Cdisc(i));
				double Ccap_r = sqrt(Ccap.x * Ccap.x + Ccap.y * Ccap.y);
				double Ccap_theta = atan2(Ccap.y, Ccap.x);

				int n_theta = 20;
				double sum1, sum2, sum3, sum4, sum5, sum6;
				sum1=0;
				sum2=0;
				sum3=0;
				sum4=0;
				sum5=0;
				sum6=0;
				double S_m;
				for(int m=0; m<=n_theta; m++)
				{
					if(m == 0 || m == n_theta)
						S_m = 1./3.;
					else if(m % 2 == 1)
						S_m = 4./3.;
					else
						S_m = 2./3.;

					double g_m1, g_m2;
					double theta_m = m*2*PI_NUMBER/n_theta;
					double a = sqrt( Ccap_r*Ccap_r + PILLAR_RADIUS*PILLAR_RADIUS - 2.*Ccap_r*PILLAR_RADIUS*cos(Ccap_theta - theta_m) + (Ccap.z - PILLAR_H1)*(Ccap.z - PILLAR_H1) );
					if(a != 0)
						g_m1 = 1./a;
					else
						g_m1 = 0;
					a = sqrt( Ccap_r*Ccap_r + PILLAR_RADIUS*PILLAR_RADIUS - 2.*Ccap_r*PILLAR_RADIUS*cos(Ccap_theta - theta_m) + (Ccap.z - PILLAR_H2)*(Ccap.z - PILLAR_H2) );
					if(a != 0)
						g_m2 = 1./a;
					else
						g_m2 = 0;
					sum1 = sum1 + S_m * cos(theta_m) * g_m1;
					sum2 = sum2 + S_m * cos(theta_m) * g_m2;

					double temp = (Ccap_r*Ccap_r + PILLAR_RADIUS*PILLAR_RADIUS - 2.*Ccap_r*PILLAR_RADIUS*cos(Ccap_theta - theta_m));
					double I_m1, I_m2;
					if( temp !=0)
					{
						I_m1 = (Ccap.z - PILLAR_H1)*g_m1 / temp;
						I_m2 = (Ccap.z - PILLAR_H2)*g_m2 / temp;
					}
					else if(Ccap.z!= PILLAR_H1 && Ccap.z!= PILLAR_H2)
					{
						I_m1 = -1./(2.*(Ccap.z - PILLAR_H1)*(Ccap.z - PILLAR_H1));
						I_m2 = -1./(2.*(Ccap.z - PILLAR_H2)*(Ccap.z - PILLAR_H2));
					}
					else
					{
						I_m1 = 0;
						I_m2 = 0;
					}
					sum3 = sum3 + S_m * cos(theta_m) * (Ccap_r - cos(Ccap_theta - theta_m)*PILLAR_RADIUS) * I_m1;
					sum4 = sum4 + S_m * cos(theta_m) * (Ccap_r - cos(Ccap_theta - theta_m)*PILLAR_RADIUS) * I_m2;
					sum5 = sum5 + S_m * cos(theta_m) * sin(Ccap_theta - theta_m) * I_m1;
					sum6 = sum6 + S_m * cos(theta_m) * sin(Ccap_theta - theta_m) * I_m2;
				}

				B.x = B.x + (sum3 - sum4) * cos(Ccap_theta) * PILLAR_RADIUS/(2.*n_theta) - (sum5 - sum6) * sin(Ccap_theta) * PILLAR_RADIUS * PILLAR_RADIUS/(2.*n_theta);
				B.y = B.y + (sum3 - sum4) * sin(Ccap_theta) * PILLAR_RADIUS/(2.*n_theta) + (sum5 - sum6) * cos(Ccap_theta) * PILLAR_RADIUS * PILLAR_RADIUS/(2.*n_theta);
				B.z = B.z + (sum2 - sum1) * PILLAR_RADIUS/(2.*n_theta);
			}
			return B;
		}

		void updateW()
		{
			for(int i = 0; i < ::c.ls.nN; i++)
			{
				//Wb(i) += b_m1(i) * vls.C.d(i);
				Wb(i) += b(i) * vls.C.d(i);
				Wc(i) += c(i) * vls.C.d(i);
				We(i) += e(i) * vls.C.d(i);
			}
//			b_m1.copy(b);
		}
		void updateA(int lsStep)
		{
			if(::c.t.lsSteps == 1 || lsStep == 0)
			{
				ba.copy(b);
				ca.copy(c);
				ea.copy(e);
				sa.copy(s);
				return;
			}
			for(int i = 0; i < ::c.ls.nN; i++)
			{
				ba(i) = (ba(i) * (double)lsStep + b(i)) / ((double)lsStep + 1.);
				ca(i) = (ca(i) * (double)lsStep + c(i)) / ((double)lsStep + 1.);
				ea(i) = (ea(i) * (double)lsStep + e(i)) / ((double)lsStep + 1.);
				sa(i) = (sa(i) * (double)lsStep + s(i)) / ((double)lsStep + 1.);
			}
		}
		void zerosB()
		{
//			b_old.copy(b);
			b.zeros();
		}
		void zerosP()
		{
			p.zeros();
		}
		void zerosS()
		{
			s.zeros();
		}
		void zerosW()
		{
			Wb.zeros();
			Wc.zeros();
			We.zeros();
		}

		void zerosC()
		{
			c.zeros();
		}

		void save2file()
		{
			save("lsf");
			b.savetxt("lsfb");
			c.savetxt("lsfc");
			e.savetxt("lsfe");
			p.savetxt("lsfp");
			s.savetxt("lsfs");
			Wb.savetxt("lsfwb");
			Wc.savetxt("lsfwc");
			We.savetxt("lsfwe");
#ifdef LS_DEBUG
			Pb.savetxt("lsfpb");
			Pc.savetxt("lsfpc");
#endif //LS_DEBUG
		};
	} F;

	class clR
	{
	public:
		clArray1Dv3d CG, V, P, O, On, L, F, T;
		clArray1Di NR, NRc;
		clArray1Dd MR;
		clArray1Dts3d II1;

		void ini()
		{
			CG.zeros(c.ls.nO);
			V.zeros(c.ls.nO);
			P.zeros(c.ls.nO);
			O.zeros(c.ls.nO);
			On.zeros(c.ls.nO);
			L.zeros(c.ls.nO);
			F.zeros(c.ls.nO);
			T.zeros(c.ls.nO);

			NR.zeros(c.ls.nO);
			NRc.zeros(c.ls.nO);
			MR.zeros(c.ls.nO);
			II1.zeros(c.ls.nO);

			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				NR(o)++;
				MR(o)+= c.ls.M(n0);
			}

			update();
		};

		void update()
		{
			calcCG();
			calcV();
			calcO();
			calcII();

			calcF();
		}

		void calcO()
		{
			O.zeros();
			On.zeros();
			L.zeros();
			NRc.copy(NR);
			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);	
				clVector3Dd r = vls.C(n0) - CG(o);
				clVector3Dd v = vls.V(n0) - V(o);
				clVector3Dd om = r ^ v;
				On(o) += om;
			}
			for(int i = 0; i < c.ls.nO; i++)
				On(i) = On(i).n();

			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				clVector3Dd r = vls.C(n0) - CG(o);
				clVector3Dd rt = On(o) * (r * On(o));
				clVector3Dd rn = r - rt;
				double rnabs2 = rn.r2();
				if(rnabs2 < 0.1)
				{
					NRc(o)--;
					continue;
				}
				clVector3Dd v = vls.V(n0) - V(o);
				clVector3Dd om = (rn ^ v) / rnabs2;
				O(o) += om;
				L(o) += om * c.ls.M(n0);
			}

			for(int i = 0; i < c.ls.nO; i++)
				if(NRc(i) > 0)
					O(i) /= (double)NRc(i);
		};

		void calcV()
		{
			V.zeros();
			P.zeros();
			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				V(o) += vls.V(n0);
				P(o) += vls.V(n0) * c.ls.M(n0);
			}

			for(int i = 0; i < c.ls.nO; i++)
				if(NR(i) > 0)
					V(i) /= (double)NR(i);
		};

		void calcCG()
		{
			CG.zeros();
			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				CG(o) += vls.C(n0) * c.ls.M(n0);
			}
			for(int i = 0; i < c.ls.nO; i++)
				if(MR(i) > 0.)
					CG(i) /= MR(i);
		};

		void calcII()
		{
			II1.zeros();
			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				clVector3Dd r = vls.C(n0) - CG(o);
				double m = c.ls.M(n0);
				II1(o).xx += (r.y * r.y + r.z * r.z) * m;
				II1(o).yy += (r.x * r.x + r.z * r.z) * m;
				II1(o).zz += (r.x * r.x + r.y * r.y) * m;
				II1(o).xy -= r.x * r.y * m;
				II1(o).xz -= r.x * r.z * m;
				II1(o).yz -= r.y * r.z * m;
			}
			for(int i = 0; i < c.ls.nO; i++)
				II1(i) = II1(i).inv();
		};
		void calcF()
		{
			F.zeros();
			T.zeros();
			for(int i = 0; i < c.ls.nTr; i++)
			{
				int n0 = c.ls.Tr(i);
				int o = c.ls.O(n0);
				double m = c.ls.M(n0);

				clVector3Dd r = vls.C(n0) - CG(o);

				F(o) += vls.F(n0);
				T(o) += r ^ vls.F(i);
			}
			//for(int i = 0; i < c.ls.nO; i++)
			//	if(NR(i) > 0.)
			//		F(i) /= NR(i);
		};

		void save2file()
		{
		};
	} R;

	//	class clTa
	//	{
	//	public:
	//		clArray2Di Tp; 
	//		clArray1Di nTp;
	//		int maxNode;
	//		double Lcut;
	//		double Sx, Sy, Smax;
	//		void ini()
	//		{
	//			maxNode = 15;
	//			Lcut = 6.;
	//			Sx = 0.;
	//			Sy = 0.;
	//			Smax = 1.;
	//			Tp.zeros(c.ls.nTa, maxNode);
	//			nTp.zeros(c.ls.nTa);
	//
	//			updateFull();
	//
	//#ifdef LS_DEBUG
	//			save2file();
	//#endif //LS_DEBUG
	//		}
	//
	//		void update()
	//		{
	//			Sx += vls.V.Xmax * c.t.dTls;
	//			Sy += vls.V.Ymax * c.t.dTls;
	//			if(Sx < Smax && Sy < Smax)
	//				return;
	//
	//			Sx = 0.;
	//			Sy = 0.;
	//
	//			updateFull();
	//		}
	//
	//		void updateFull()
	//		{
	//			for(int i = 0; i < c.ls.nTa; i++)
	//			{
	//				int n0 = c.ls.Ta(i);
	//				double x0 = vls.C.x(n0), y0 = vls.C.y(n0);
	//				int np = 0;
	//				for(int j = 0; j < c.ls.nTp; j++)
	//				{
	//					int n1 = c.ls.Tp(j);
	//					double x1 = vls.C.x(n1), y1 = vls.C.y(n1);
	//					double dX = x1 - x0, dY = y1 - y0;
	//					double r = sqrt(dX * dX + dY * dY);
	//					if(r > Lcut)
	//						continue;
	//
	//					if(np == maxNode)
	//					{
	//						//p.error("Too many nodes");
	//						maxNode += 5;
	//						Tp.zeros(c.ls.nTa, maxNode);
	//						updateFull();
	//						return;
	//					}
	//
	//					Tp(i,np++) = n1;
	//				}
	//				nTp(i) = np;
	//			}
	//		};
	//		void save2file()
	//		{
	//			Tp.savetxt("lstatp");
	//			nTp.savetxt("lstantp");
	//		}
	//	} Ta;
	class clP
	{
	public:
		//		class clN: public clArray1Di
		//		{
		//		public:
		//			int nN, nA;
		//			clN()
		//			{
		//				ini();
		//			}
		//
		//			void ini()
		//			{
		//				nN = 0;
		//				nA = 0;
		//			}
		//			void increase(int nIncr)
		//			{
		//				clArray1Di::reallocate(nA + nIncr);
		//				nA += nIncr;
		//			}
		//			void add(int n)
		//			{
		//				if(nN >= nA)
		//					increase(LS_ARRAY_INCREASE_STEP);
		//				el(nN++) = n;
		//			}
		//			void reset()
		//			{
		//				nN = 0;
		//			}
		//			int find(int n)
		//			{
		//				for(int i = 0; i < nN; i++)
		//					if(el(i) == n)
		//						return i;
		//				return -1;
		//			}
		////#ifdef LS_DEBUG
		////			int& clN::operator()(const int n1)
		////			{
		////				if(n1 >= nN || n1 >= nA || n1 < 0)
		////					printf("Error in clN, %d from %d tot %d", n1, nN, nA);
		////				return el(n1);
		////			}
		////#endif
		//		};
		clArray3D <clArray1Dsti> M;
		clArray1D <clArray1Dsti> N;
		clArray3D <clArray1Dsti> Ma;
		clArray1D <clArray1Dsti> Na;
		clArray1D <clArray1Dsti> BL;
		clArray1D <clArray1Dsti> BC;
		clArray1Dd F;
		int nBC;
		int step_counter;
		int stepR, stepA;
		int updR, updB;
		clVector3Di nnR, nnA;
		BOOL flRepulsionB;
		void ini()
		{
			updR = 0;
			updB = updR - 1;

			force_update();
			flRepulsionB = FALSE;

			stepR = LS_SEARCH_REPULSIVE_STEP;
			nnR = c.lb.nn / stepR;
			M.allocate(nnR);
			N.allocate(c.ls.nN);
			BL.allocate(c.ls.nN);

			nBC = c.ls.nN;
			BC.allocate(nBC);

			stepA = LS_SEARCH_ACTIVE_STEP;
			nnA = c.lb.nn / stepA;
			Ma.allocate(nnA);
			Na.allocate(c.ls.nN);
		};
		void force_update()
		{
			step_counter = c.ls.nPupd;
		}
		void update()
		{
#ifdef LS_POTENCIAL
			if(++step_counter < c.ls.nPupd)
				return;
			step_counter = 0;

			updateM();

			updateA();
#endif //LS_POTENCIAL
		};
		void updateF()
		{
#ifdef LS_POTENCIAL
			force_repulsion();
			force_repulsionB();
			force_repulsion_wall();
			force_attraction();
#endif //LS_POTENCIAL
#ifdef LS_BEADS_FRICTION
			clArray1Dv3d Fe;
			Fe.zeros(c.ls.nO);

			for(int i = 0; i < c.ls.nN; i++)	// Apply forces
			{
				int o = c.ls.O(i);
				Fe(o) += vls.F.e(i);
			}
			for(int i = 0; i < c.ls.nN; i++)
			{
				int o = c.ls.O(i);
				if (vls.C(i).z < -0.4 && o < N_BEAD)	//rolling, in contact with wall
				{
					vls.F.c(i) = - vls.V(i).n() * (-Fe(o).z) * 0.4;
					vls.F.c(i).z = 0.;
					//break;
				}
				
				if (o < N_BEAD)		// Test if node is part of a bead
				{
				vls.F.c(i) += - vls.V(i).n() * TANG_FORCE / N0_NUMBER;	// Add uniformly distributed force to bead in opposite direction of velocity
				vls.F.c(i).z = 0.;
				}
			}
#endif // LS_BEADS_FRICTION
//#ifdef SET_ORBIT
//			for(int i = 0; i < c.ls.nN; i++)	// Apply force to cause circular orbit
//			{
//			int o = c.ls.O(i);	// Bead object number
//			int od = ceil((o + 1) / 2.) + N_BEAD;	// Determine object number of corresponding disk
//			clVector3Dd R_bd = vls.avr.CG(od) - vls.avr.CG(od);	// Vector from bead center to disk center
//			vls.F.c(i) = c.ls.M(i) * BEAD_OMEGA_SQ * R_bd - (vls.F.s(i) + vls.F.b(i) + vls.F.e(i));
//			}
//#endif // SET_ORBIT
		};
		void updateM()
		{	
			updR++;

			for(int i = 0; i < nnR.x; i++)
				for(int j = 0; j < nnR.y; j++)
					for(int k = 0; k < nnR.z; k++)
						M(i,j,k).reset();

			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n = c.ls.BN(i);
				if(c.ls.N.NL(n) > 0)
					M(indexR(vls.C.s(n))).add(n);
			}

			for(int i = 0; i < c.ls.nN; i++)
				N(i).reset();

			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n = c.ls.BN(i);
				if(c.ls.N.NL(n) > 0)
					findM(n);
			}
		};
		void updateA()
		{
			if(c.ls.nTa == 0)
				return;

			for(int i = 0; i < nnA.x; i++)
				for(int j = 0; j < nnA.y; j++)
					for(int k = 0; k < nnA.z; k++)
						Ma(i,j,k).reset();

			for(int i = 0; i < c.ls.nTa; i++)
			{
				int n = c.ls.Ta(i);
				if(c.ls.N.NL(n) > 0)
					Ma(indexA(vls.C.s(n))).add(n);
			}

			for(int i = 0; i < c.ls.nN; i++)
				Na(i).reset();

			for(int i = 0; i < c.ls.nTa; i++)
			{
				int n = c.ls.Ta(i);
				if(c.ls.N.NL(n) > 0)
					findA(n);
			}
		};
		void updateB()
		{
			updB = updR;

			for(int i = 0; i < c.ls.nN; i++)
				BL(i).reset();

			for(int i = 0; i < c.ls.nBL; i++)
			{
				BL(c.ls.BL(i,0)).add(i);
				BL(c.ls.BL(i,1)).add(i);
				BL(c.ls.BL(i,2)).add(i);
			}
			if(c.ls.nBL > nBC)
			{
				nBC = c.ls.nBL + 100;
				BC.allocate(nBC);
			}
			for(int i = 0; i < c.ls.nBL; i++)
				BC(i).reset();

			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n0 = c.ls.BN(i);

				clArray1Dsti &N0 = N(n0);
				clArray1Dsti &BL0 = BL(n0);

				for(int j = 0; j < N0.nN; j++)
				{
					int n1 = N0(j);
					clArray1Dsti &BL1 = BL(n1);
					for(int k = 0; k < BL0.nN; k++)
						for(int m = 0; m < BL1.nN; m++)
							if(BC(BL0(k)).find(BL1(m)) < 0)
								BC(BL0(k)).add(BL1(m));
				}
			}
			//testBL();
		};
		void testBL()
		{
			for(int i = 0; i < c.ls.nBL; i++)
			{
				clArray1Dsti &BC0 = BC(i);
				for(int j = 0; j < BC0.nN; j++)
					if(BC0(j) >= c.ls.nBL)
						printf("\nError in BC %d, %d of %d", i, BC0(j), c.ls.nBL);
			}
		};
		clVector3Di indexR(clVector3Dd cc0)
		{
			return index(cc0, stepR, nnR);
		};
		clVector3Di indexA(clVector3Dd cc0)
		{
			return index(cc0, stepA, nnA);
		};
		clVector3Di index(clVector3Dd cc0, double step, clVector3Di size)
		{
			return (clVector3Di)MOD((cc0 + clVector3Dd(0.5, 0.5, 0.5)) / step, size);
		};
		void findM(int n0)
		{
			const int size = LS_SEARCH_REPULSIVE_SIZE;
			clVector3Di ii0 = indexR(vls.C.s(n0));

			clArray1Dsti &Nd = N(n0);
			//			int O0 = c.ls.O(n0);
			int Op0 = c.ls.Op(n0);

			for(int i = -size; i <= size; i++)
				for(int j = -size; j <= size; j++)
					for(int k = -size; k <= size; k++)
					{
						clArray1Dsti &MNs = M(MOD(ii0 + clVector3Di(i,j,k),nnR));
						for(int m = 0; m < MNs.nN; m++)
						{
							int n1 = MNs(m);
							//							int O1 = c.ls.O(n1);
							int Op1 = c.ls.Op(n1);
							//							if(n1 > n0 && c.ls.N.find_next(n0, n1) == FALSE)
							if(n1 > n0 && (Op0 != Op1 || vls.C.p(n0) != vls.C.p(n1) || c.ls.N.find_next(n0, n1) == FALSE || c.ls.NPN(n0) != n1))
								Nd.add(n1);
						}
					}
		};
		void findA(int n0)
		{
			const int size = LS_SEARCH_ACTIVE_SIZE;
			clVector3Di ii0 = indexA(vls.C.s(n0));

			clArray1Dsti &Nd = Na(n0);
			int O0 = c.ls.O(n0);
			//			int Op0 = c.ls.Op(n0);

			for(int i = -size; i <= size; i++)
				for(int j = -size; j <= size; j++)
					for(int k = -size; k <= size; k++)
					{
						clArray1Dsti &MNs = Ma(MOD(ii0 + clVector3Di(i,j,k),nnA));
						for(int m = 0; m < MNs.nN; m++)
						{
							int n1 = MNs(m);
							int O1 = c.ls.O(n1);
							//							int Op1 = c.ls.Op(n1);
							//							if(n1 > n0 && c.ls.N.find_next(n0, n1) == FALSE)
							//							if(n1 > n0 && (O0 != O1 || Op0 != Op1))
							if(n1 > n0 && O0 != O1)
								Nd.add(n1);
						}
					}
		};
		void force_repulsion()
		{
			flRepulsionB = FALSE;
			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n0 = c.ls.BN(i);
				if(c.ls.N.NL(n0) == 0)
					continue;

				clVector3Dd c0 = vls.C.s(n0);
				int o0 = c.ls.O(n0);
				clArray1Dsti &N0 = N(n0);

				for(int j = 0; j < N0.nN; j++)
				{
					int n1 = N0(j);
					clVector3Dd c1 = vls.C.s(n1);
					int o1 = c.ls.O(n1);
					clForceV Fv;

					Fv = force_repulsionC(c0, c1, n0, n1);

					if(o0 != o1 || (vls.C.p(n0) != vls.C.p(n1) && c.ls.NPN(n0) != n1))
						Fv = Fv + force_repulsionR(c0, c1, n0, n1);

					set_force(Fv, n0, n1);
				}
			}
		};
		void force_repulsionB()
		{
			if(flRepulsionB == FALSE)
				return;

			//testBL();
			if(updR > updB)
				updateB();


			const int add_nodes = 7;
			clVector3Dd C0[add_nodes], C1[add_nodes];
			clVector3Dd C0s[NUMBER_DIMENSIONS], C1s[NUMBER_DIMENSIONS];

			int N0[NUMBER_DIMENSIONS],  N1[NUMBER_DIMENSIONS];
			double D[add_nodes][NUMBER_DIMENSIONS] = {		{0.5, 0.25, 0.25},
			{0.25, 0.5, 0.25},
			{0.25, 0.25, 0.5},
			{1./3., 1./3., 1./3.},
			{0.5, 0.5, 0.},
			{0.5, 0., 0.5},
			{0., 0.5, 0.5}
			};
			for(int i = 0; i < c.ls.nBL; i++)
			{
				clArray1Dsti &BC0 = BC(i);
				if(BC0.nN == 0)
					continue;

				for(int m = 0; m < NUMBER_DIMENSIONS; m++)
				{
					N0[m] = c.ls.BL(i,m); 
					C0s[m] = vls.C.s(N0[m]);
				}
				for(int m = 0; m < add_nodes; m++)
				{
					C0[m] = 0.;
					for(int q = 0; q < NUMBER_DIMENSIONS; q++)
						C0[m] += C0s[q] * D[m][q];
				}

				for(int j = 0; j < BC0.nN; j++)
				{
					int n1 = BC0(j);
					for(int m = 0; m < NUMBER_DIMENSIONS; m++)
					{
						N1[m] = c.ls.BL(n1,m); 
						C1s[m] = vls.C.s(N1[m]);
					}
					for(int m = 0; m < add_nodes; m++)
					{
						C1[m] = 0.;
						for(int q = 0; q < NUMBER_DIMENSIONS; q++)
							C1[m] += C1s[q] * D[m][q];
					}
					for(int m0 = 0; m0 < add_nodes; m0++)
						for(int m1 = 0; m1 < add_nodes; m1++)
						{
							clForceV Fv = force_repulsionC(C0[m0], C1[m1], m0, m1);

							if(Fv.p == 0.)
								continue;

							set_forceP(Fv * D[m0][0], N0[0]);
							set_forceP(Fv * D[m0][1], N0[1]);
							set_forceP(Fv * D[m0][2], N0[2]);

							set_forceM(Fv * D[m1][0], N1[0]);
							set_forceM(Fv * D[m1][1], N1[1]);
							set_forceM(Fv * D[m1][2], N1[2]);
						}
						for(int m0 = 0; m0 < add_nodes; m0++)
							for(int q = 0; q < NUMBER_DIMENSIONS; q++)
							{
								clForceV Fv = force_repulsionC(C0[m0], C1s[q], m0, q);

								if(Fv.p == 0.)
									continue;

								set_forceP(Fv * D[m0][0], N0[0]);
								set_forceP(Fv * D[m0][1], N0[1]);
								set_forceP(Fv * D[m0][2], N0[2]);

								set_forceM(Fv, N1[q]);
							}
							for(int m1 = 0; m1 < add_nodes; m1++)
								for(int q = 0; q < NUMBER_DIMENSIONS; q++)
								{
									clForceV Fv = force_repulsionC(C0s[q], C1[m1], q, m1);

									if(Fv.p == 0.)
										continue;

									set_forceP(Fv, N0[q]);

									set_forceM(Fv * D[m1][0], N1[0]);
									set_forceM(Fv * D[m1][1], N1[1]);
									set_forceM(Fv * D[m1][2], N1[2]);
								}
				}
			}
		};
		void force_repulsion_wall()
		{
			if(c.lb.bcWallX0 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallX0(i);

			if(c.lb.bcWallX1 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallX1(i);

			if(c.lb.bcWallY0 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallY0(i);

			if(c.lb.bcWallY1 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallY1(i);

			if(c.lb.bcWallZ0 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallZ0(i);

			if(c.lb.bcWallZ1 == TRUE)
				for(int i = 0; i < c.ls.nBN; i++)
					FwallZ1(i);
		};
		void FwallX0(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double xw = -0.5;
			double r = vls.C.s(na1).x - xw;
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeX0 + c.ls.TaDe(ia));

			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRX0);

			vls.F.c(na1).x -= Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;
		};
		void FwallX1(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double xw = c.lb.nX - 0.5;
			double r = -(vls.C.s(na1).x - xw);
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeX1 + c.ls.TaDe(ia));

			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRX1);

			vls.F.c(na1).x += Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;
		};
		void FwallY0(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double yw = -0.5;
			double r = vls.C.s(na1).y - yw;
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeY0 + c.ls.TaDe(ia));


			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRY0);

			vls.F.c(na1).y -= Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;
		};
		void FwallY1(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double yw = c.lb.nY - 0.5;
			double r = -(vls.C.s(na1).y - yw);
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeY1 + c.ls.TaDe(ia));

			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRY1);


			vls.F.c(na1).y += Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;
		};
		void FwallZ0(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double zw = -0.5;
			double r = vls.C.s(na1).z - zw;
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeZ0 + c.ls.TaDe(ia));

			double v = sqrt(vls.V(na1).x2() + vls.V(na1).y2());
			clVector3Dd v1 = vls.V(na1).n() * v;
			double alf = 0;
			if( r < R_EQ_NUMBER)
				alf = RATIO / r;
			clVector3Dd Ff = v1 * alf;

			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRZ0);


			vls.F.c(na1).z -= Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;

			vls.F.c(na1).x -= Ff.x;
			vls.F.c(na1).y -= Ff.y;
		};
		void FwallZ1(int ia)
		{
			int na1 = c.ls.BN(ia);
			if(!(c.ls.T(na1) & LS_WALL_INTER_NODE))
				return;

			double zw = c.lb.nZ - 0.5;
			double r = -(vls.C.s(na1).z - zw);
			if(r <= 0. || r > LS_SEARCH_ACTIVE_SIZE * LS_SEARCH_ACTIVE_STEP)
				return;

			clForce Fa, Fr;
			if(c.ls.T(na1) & LS_ACTIVE_NODE)
				Fa = Fattract(r, c.ls.FwDeZ1 + c.ls.TaDe(ia));

			Fr = FrepulsC(r) + Frepuls(r, c.ls.FwDeRZ1);


			vls.F.c(na1).z += Fa.f + Fr.f;
			vls.F.p(na1) += Fa.p + Fr.p;
		};
		void force_attraction()
		{
			for(int i = 0; i < c.ls.nTa; i++)
			{
				int n0 = c.ls.Ta(i);
				if(c.ls.N.NL(n0) == 0)
					continue;

				clVector3Dd c0 = vls.C.s(n0);
				clArray1Dsti &N0 = Na(n0);

				for(int j = 0; j < N0.nN; j++)
				{
					int n1 = N0(j);
					clVector3Dd c1 = vls.C.s(n1);
					clForceV Fv = force_attraction(c0, c1, n0, n1);
					set_force(Fv, n0, n1);
				}
			}
		};
		clForceV force_repulsionR(clVector3Dd c0, clVector3Dd c1, int n0, int n1)
		{
			clVector3Dd d = c1 - c0;
			double r = d.r();

			clForce Fr = Frepuls(r);

			return clForceV(d * Fr.f / r, Fr.p);
		};
		clForceV force_repulsionC(clVector3Dd c0, clVector3Dd c1, int n0, int n1)
		{
			clVector3Dd d = c1 - c0;
			double r = d.r();

			if(r < 0.1 * c.ls.FwReC)
			{
				if(c.ls.NPN(n0) == n1)
					return clForceV();

				printf("\nRepulsion: two nodes are too close: r=%f\tn0=%d\tn1=%d", r, n0, n1);
				if(r <= 0.)
					return clForceV();
			}

			if(r < 1.25 * c.ls.FwReC)
				flRepulsionB = TRUE;
			else
				return clForceV();

			clForce Fc = FrepulsC(r);

			return clForceV(d * Fc.f / r, Fc.p);
		};
		clForceV force_attraction(clVector3Dd c0, clVector3Dd c1, int n0, int n1)
		{
			clVector3Dd d = c1 - c0;
			double r = d.r();

			if(r < 0.05 * c.ls.FwRe)
			{
				if(c.ls.NPN(n0) == n1)
					return clForceV();

				printf("\nAttraciton: nodes are too close: r=%f\tn0=%d\tn1=%d", r, n0, n1);
				if(r <= 0.)
					return clForceV();
			}

			clForce Fc = Fattract(r, c.ls.TaDe(n0) + c.ls.TaDe(n1));

			return clForceV(d * Fc.f / r, Fc.p);
		};
		void set_force(clForceV Fc, int n0, int n1)
		{
			set_forceP(Fc, n0);
			set_forceM(Fc, n1);
		};
		void set_forceP(clForceV Fc, int n)
		{
			vls.F.p(n) += Fc.p;
			vls.F.c(n) += Fc.f;
		};
		void set_forceM(clForceV Fc, int n)
		{
			vls.F.p(n) += Fc.p;
			vls.F.c(n) -= Fc.f;
		};
		clForce Frepuls(double r)
		{
			return Frepuls(r, c.ls.FwDeR);
		};
		clForce Frepuls(double r, double De)
		{
			if(r > c.ls.FwRe)
				return clForce();
			return Fmorse(r, De, c.ls.FwRe, c.ls.FwKe, 0.);
		};
		clForce FrepulsC(double r)
		{
			if(r > c.ls.FwReC)
				return clForce();
			return Fmorse(r, c.ls.FwDeC, c.ls.FwReC, c.ls.FwKe, 0.);
		};
		clForce Fattract(double r, double De)
		{
			return Fmorse(r, De, c.ls.FwRe, c.ls.FwKe, De);
		};
		clForce Fmorse(double r, double De, double re, double ke, double De_inf)
		{
			clForce ret;
			if(De == 0.)
				return ret;

			double e = exp(-(r - re)/ke);
			double e1 = 1. - e;
			ret = clForce(2. * De/ke * e1 * e, De * e1 * e1 - De_inf);
			return ret;
		};
		void save2file()
		{
		};

	} P;


	class clAvr
	{
	public:
		clArray1Dv3d CG, Xmax, Xmin, V, P, Fb, Fc, Fe, O, On, L;
		clArray1Dd Ravr, Rmax, Rmin, Vmax;
		clArray1Dd Ek, Ec, Ev, Smax, Smin;
		clArray1Dd Wb, Wc, We;
		clArray1Di NOc;

		void ini()
		{
			CG.zeros(c.ls.nO);
			Xmax.zeros(c.ls.nO);
			Xmin.zeros(c.ls.nO);
			Ravr.zeros(c.ls.nO);
			V.zeros(c.ls.nO);
			Vmax.zeros(c.ls.nO);
			P.zeros(c.ls.nO);
			O.zeros(c.ls.nO);
			On.zeros(c.ls.nO);
			L.zeros(c.ls.nO);
			Fb.zeros(c.ls.nO);
			Fc.zeros(c.ls.nO);
			Fe.zeros(c.ls.nO);
			Rmax.zeros(c.ls.nO);
			Rmin.zeros(c.ls.nO);
			Ek.zeros(c.ls.nO);
			Ec.zeros(c.ls.nO);
			Ev.zeros(c.ls.nO);
			Smax.zeros(c.ls.nO);
			Smin.zeros(c.ls.nO);
			Wb.zeros(c.ls.nO);
			Wc.zeros(c.ls.nO);
			We.zeros(c.ls.nO);

			NOc.zeros(c.ls.nO);

			update();
		};

		void update()
		{
			calcNO();
			calcCG();
			calcRavr();
			calcV();
			calcOmega();
			calcF();
			calcE();
			calcS();
			calcW();
		}

		void calcNO()
		{
			c.ls.NO.zeros(c.ls.nO);
			for(int i = 0; i < c.ls.nN; i++)
				if(c.ls.N.NL(i) > 0)
					c.ls.NO(c.ls.O(i))++;
		}

		void calcS()
		{
			Smax.zeros();
			Smin.zeros();

			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				for(int j = 0; j < c.ls.nL; j++)
				{
					int n = c.ls.N(i,j);
					if(n >= 0 && i < n)
					{
						double S = vls.C.S(i,j);
						if(S > Smax(o))
							Smax(o) = S;
						else if(S < Smin(o))
							Smin(o) = S;
					}
				}
			}
		}

		void calcE()
		{
			Ek.zeros();
			Ev.zeros();

			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				double Ekcur = 0.;
				for(int j = 0; j < c.ls.nL; j++)
				{
					int n = c.ls.N(i,j);
					if(n >= 0 && i < n)
						Ekcur += c.ls.K.E(i,j);
				}
				clVector3Dd V = vls.V(i);
				double Evcur = 0.5 * c.ls.M(i) * V.r2();

				int o = c.ls.O(i);
				Ek(o) += Ekcur;
				Ev(o) += Evcur;
			}
			for(int i = 0; i < c.ls.nAN; i++)
			{
				int o = c.ls.O(c.ls.AN(i,1));
				Ek(o) += c.ls.AK.E(i);
			}
			calcEc();
		}

		void calcEc()
		{
			Ec.zeros();

			for(int i = 0; i < c.ls.nBN; i++)
			{
				int n = c.ls.BN(i);
				Ec(c.ls.O(n)) += vls.F.p(n);
			}
		}


		void calcF()
		{
			Fb.zeros();
			Fc.zeros();
			Fe.zeros();
			for(int i = 0; i < c.ls.nN; i++)
			{
				if((c.ls.T(i) & LS_PERIODIC_NODE) || c.ls.N.NL(i) == 0)
					continue;

				int o = c.ls.O(i);
				Fb(o) += vls.F.ba(i);
				Fc(o) += vls.F.ca(i);
				Fe(o) += vls.F.ea(i);
			}
		};

		void calcOmega()
		{
			O.zeros();
			On.zeros();
			L.zeros();
			NOc.copy(c.ls.NO);
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				clVector3Dd r = vls.C(i) - CG(o);
				clVector3Dd v = vls.V(i) - V(o);
				clVector3Dd om = r ^ v;
				On(o) += om;
			}
			for(int i = 0; i < c.ls.nO; i++)
				On(i) = On(i).n();

			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				clVector3Dd r = vls.C(i) - CG(o);
				clVector3Dd rt = On(o) * (r * On(o));
				clVector3Dd rn = r - rt;
				//clVector3Dd rn = r;
				double rnabs2 = rn.r2();
				if(rnabs2 < c.eps)
				{
					NOc(o)--;
					continue;
				}
				clVector3Dd v = vls.V(i) - V(o);
				clVector3Dd om = (rn ^ v) / rnabs2;
				O(o) += om;
				L(o) += om * c.ls.M(i);
			}


			for(int i = 0; i < c.ls.nO; i++)
				if(NOc(i) > 0)
					O(i) /= (double)NOc(i);
		};

		void calcV()
		{
			Vmax.zeros();
			V.zeros();
			P.zeros();
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				double Vmax_tmp = vls.V(i).abs();
				if(Vmax_tmp > Vmax(o))
					Vmax(o) = Vmax_tmp;
				V(o) += vls.V(i);
				P(o) += vls.V(i) * c.ls.M(i);
			}

			for(int i = 0; i < c.ls.nO; i++)
				if(c.ls.NO(i) > 0)
					V(i) /= (double)c.ls.NO(i);
		};

		void calcCG()
		{
			CG.zeros();
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				CG(o) += vls.C(i) * c.ls.M(i);
			}
			for(int i = 0; i < c.ls.nO; i++)
				if(c.ls.M.O(i) > 0)
					CG(i) /= c.ls.M.O(i);
		};

		void calcRavr()
		{
			const double max_val = 1e10, min_val = -1e10;
			Ravr.zeros();
			Rmax.zeros();
			Rmin.fill(c.nX + c.nY);

			Xmax.fill(clVector3Dd(min_val, min_val, min_val));
			Xmin.fill(clVector3Dd(max_val, max_val, max_val));

			double R, Ra = 0;
			int o;
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				o = c.ls.O(i);
				clVector3Dd Ccur = vls.C(i);
				R = clVector3Dd(Ccur - CG(o)).r();
				if(R > Rmax(o))
					Rmax(o) = R;
				if(R < Rmin(o))
					Rmin(o) = R;

				Ravr(o) += R;
				if(Ccur.x > Xmax(o).x)
					Xmax(o).x = Ccur.x;
				if(Ccur.y > Xmax(o).y)
					Xmax(o).y = Ccur.y;
				if(Ccur.z > Xmax(o).z)
					Xmax(o).z = Ccur.z;

				if(Ccur.x < Xmin(o).x)
					Xmin(o).x = Ccur.x;
				if(Ccur.y < Xmin(o).y)
					Xmin(o).y = Ccur.y;
				if(Ccur.z < Xmin(o).z)
					Xmin(o).z = Ccur.z;
			}
			for(int i = 0; i < c.ls.nO; i++)
			{
				if(c.ls.NO(i) > 0)
				{
					Ravr(i) /= (double)c.ls.NO(i);
				}
				else
				{
					Rmax(o) = 0.;
					Rmin(o) = 0.;
				}
			}
		};
		void calcW()
		{
			Wc.zeros();
			Wb.zeros();
			We.zeros();
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.N.NL(i) == 0)
					continue;
				int o = c.ls.O(i);
				Wc(o) += vls.F.Wc(i);
				Wb(o) += vls.F.Wb(i);
				We(o) += vls.F.We(i);
			}
			vls.F.zerosW();
		};
	} avr;

};
////////////////////////////////////////////////////////////////////////////
class clVariablesLB
{
public:
	clVariablesLB(){};
	virtual ~clVariablesLB(){};

	void ini()
	{
		M.ini();

		F.ini();
		Ro.ini();
		Ro.update();
		J.ini();
		PP.ini();
		B.ini();
#ifdef BINARY_FLUID
		G.ini();
		N.ini();
#endif		
		
		M.update();


#ifdef INI_SHEAR
		F.inishear();
#endif //INI_SHEAR
#ifdef INI_POISEUILLE
		F.iniPoiseuille();
#endif //INI_POISEUILLE
		//F.iniflow();

//		Ro.update();
		J.update();
		PP.update();
#ifdef BINARY_FLUID
		N.update();
#endif

		updateoutput();

#ifdef LB_COARSE
		C.ini();
//		C.update();
#endif //LB_COARSE
	};

	void update()
	{
	};
	void updateoutput()
	{
		Ro.updateoutput();
		J.updateoutput();
	};

	void updateoutputmix()
	{

	};

	class clF3D: public clArray3Dd
	{
	protected:
		//#ifdef LB_THIN_WALL
		//		clArray3Dd na;
		//#endif //LB_THIN_WALL

		int ind;
		int shiftX, shiftY, shiftZ;
		int shiftXn, shiftYn, shiftZn;
		int incX, incY, incZ;
		int nX, nY, nZ;
		int nXg, nYg, nZg;
		char FileNameDistr[80], FileName[80];

		void make_copy(clArray3Dd &copy)
		{
			copy.zeros(nXg, nYg, nZg);
			for(int i = 0; i < nXg; i++)
				for(int j = 0; j < nYg; j++)
					for(int m = 0; m < nZg; m++)
						copy(i,j,m) = ff(i,j,m);

		}

		void make_copy_full(clArray3Dd &copy)
		{
			copy.zeros(nX, nY, nZ);
			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int m = 0; m < nZ; m++)
						copy(i,j,m) = ff(i,j,m);

		}

		int propagation(int x, int inc, int nx)
		{
			x -= inc;
			if(x >= nx)
				x -= nx;
			else if(x < 0)
				x += nx;
			return x;
		}
	public:
		void ini(char* Name, int n, clVector3Di size)
		{
			ind = n;
			incX = c.lb.Cx[ind];
			incY = c.lb.Cy[ind];
			incZ = c.lb.Cz[ind];
			nXg = size.x;
			nYg = size.y;
			nZg = size.z;
			sprintf(FileName, "%s", Name);
			sprintf(FileNameDistr, "%s%d", Name, ind);
			ini();
		}
		virtual void ini()
		{
			nX = nXg + 1;
			nY = nYg + 1;
			nZ = nZg + 1;

			shiftXn = 0;
			shiftYn = 0;
			shiftZn = 0;

			propagation();

			zeros(nX, nY, nZ);
			fill(c.lb.F0[ind]);
			//#ifdef LB_THIN_WALL
			//			na.zeros(nX, nY, nZ);
			//			na.fill(c.lb.F0[ind]);
			//#endif //LB_THIN_WALL

			load();
			//copy();
			savetxt();
		}

		void load()
		{

			char Buf[80], Bufbin[80];
			sprintf(Buf, "load" SLASH "%s", FileNameDistr);
			sprintf(Bufbin, "load" SLASH "%s.bin", FileName);

			clArray3Dd loadArray;
			loadArray.zeros(nXg, nYg, nZg);

			FILE *stream = fopen(Bufbin, "r+b" );
			BOOL retbin;
			for(int i = 0; i <= ind; i++)
				retbin = loadArray.loadbin(stream);

			if(retbin == TRUE)
				fclose(stream);


			if((retbin == FALSE) && (loadArray.loadtxt(Buf) == FALSE))
			{
				if ( strncmp(FileName,"lbg",3) == 0 )
					loadg();
				else
					loadf();

			}
			else
			{
				for(int i = 0; i < nXg; i++)
					for(int j = 0; j < nYg; j++)
						for(int m = 0; m < nZg; m++)
							el(i,j,m) = loadArray(i,j,m);
			}

		}

		void loadf() 
		{
			for(int i = 0; i < nXg; i++)
				for(int j = 0; j < nYg; j++)
					for(int m = 0; m < nZg; m++)
						el(i,j,m) = c.lb.F0[ind];
		}
		void loadg() {
			clArray3Dd loadN;
			loadN.zeros(nXg, nYg, nZg);
			if(loadN.loadtxt("load" SLASH "lbn") == FALSE)
				loadN.fill(0.);
			
			for(int i = 0; i < nXg; i++)
				for(int j = 0; j < nYg; j++)
					for(int k = 0; k < nZg; k++)
					{
						double G;
						clVector3Dd  J = vlb.J(i,j,k); 
						double Ro = vlb.Ro(i,j,k);
						double N = loadN(i,j,k);
						clVector3Dd U = J / Ro;
						clVector3Dd C = c.lb.C[ind];



						if(Ro == 0.)
							continue;
#ifdef BINARY_SCALAR
						G = c.lb.F0[ind] * N;
#else
						if(ind > 0)
						{

							//double Ux = vlb.J(i,j,k).x / Ro, Uy = vlb.J.y(i,j) / Ro;
							//double U = (c.lb.Cx[ind] * Ux + c.lb.Cy[ind] * Uy);
							double Gamma = LBB_GAMMA_NUMBER;
							//double d2N = vlb.N.d2(i,j,k).sum();
							double GMu = Gamma * ((-LBB_A + LBB_B * N * N) * N );
							//								double GMu = c.lb.G * ((-c.lb.a + c.lb.b * N * N) * N);
							G = c.lb.F0[ind] * (GMu + U * C * N) / c.lb.Cs2;

						}
						else
						{
							double Gtot = 0.;
							for(int n = 1; n < c.lb.nL; n++)
								Gtot += vlb.G(n)(i,j,k);
							G = N - Gtot;

						}
#endif // BINARY_SCALAR
						el(i,j,k) = G;
					}
					//char Buftxt[80], Bufbin[80];
					//sprintf(Buftxt, "load" SLASH "%s.txt", FileName);
					//sprintf(Bufbin, "load" SLASH "%s.bin", FileName);

					//clArray3Dd loadArray;
					//loadArray.zeros(nXg, nYg, nZg);

					//// if text file doesn't exist, assume we have no Species-B fluid
					//if (loadArray.loadtxt(Buftxt) == FALSE) {
					//	loadArray.fill(0);
					//	return;
					//}

					//// otherwise populate the N values
					//for(int i = 0; i < nXg; i++)
					//	for(int j = 0; j < nYg; j++)
					//		for(int m = 0; m < nZg; m++) {
					//			el(i,j,m) = loadArray(i,j,m);

					//		}
		}


		void update()
		{
		}

		void copy()
		{
		}

		void change()
		{
		}

		void propagation()
		{
			shiftX = shiftXn;
			shiftY = shiftYn;
			shiftZ = shiftZn;

			shiftXn = propagation(shiftXn, incX, nX);
			shiftYn = propagation(shiftYn, incY, nY);
			shiftZn = propagation(shiftZn, incZ, nZ);
		}

		double& f(int i, int j, int m)
		{
			i = MODFAST(i, nX);
			j = MODFAST(j, nY);
			m = MODFAST(m, nZ);

			return ff(i,j,m);
		}

		double& ff(int i, int j, int m)
		{
			int is, js, ms;
			is = MODFAST(i + shiftX, nX);
			js = MODFAST(j + shiftY, nY);
			ms = MODFAST(m + shiftZ, nZ);

			return el(is,js,ms);
		}

		double& f(clVector3Di ii)
		{
			return f(ii.x, ii.y, ii.z);
		}

		double& ff(clVector3Di ii)
		{
			return ff(ii.x, ii.y, ii.z);
		}
		//#ifdef LB_THIN_WALL
		//		double& nc(clVector3Di ii)
		//		{
		//			return nc(ii.x, ii.y, ii.z);
		//		}
		//		double& nc(int i, int j, int m)
		//		{
		//			int is, js, ms;
		//			is = MOD(i + shiftX, nX);
		//			js = MOD(j + shiftY, nY);
		//			ms = MOD(m + shiftZ, nZ);
		//
		//			return na(is,js,ms);
		//		}
		//#endif //LB_THIN_WALL

		double& n(clVector3Di ii)
		{
			return n(ii.x, ii.y, ii.z);
		}

		double& n(int i, int j, int m)
		{
			i = MODFAST(i, nX);
			j = MODFAST(j, nY);
			m = MODFAST(m, nZ);
			int is, js, ms;
			is = MODFAST(i + shiftXn, nX);
			js = MODFAST(j + shiftYn, nY);
			ms = MODFAST(m + shiftZn, nZ);
			//#ifndef LB_THIN_WALL
			return el(is,js,ms);
			//#else 
			//			return na(is,js,ms);
			//#endif //LB_THIN_WALL
		}

		double& npx(int i, int j, int z)
		{
			j = MODFAST(j + incY, nYg);
			z = MODFAST(z + incZ, nZg);
			return n(i,j,z);
		}

		double& npy(int i, int j, int z)
		{
			i = MODFAST(i + incX, nXg);
			z = MODFAST(z + incZ, nZg);
			return n(i,j,z);
		}

		double& npz(int i, int j, int z)
		{
			i = MODFAST(i + incX, nXg);
			j = MODFAST(j + incY, nYg);
			return n(i,j,z);
		}

		double& operator()(int i, int j, int m)
		{
			return f(i,j,m);
		}

		double& operator()(clVector3Di ii)
		{
			return f(ii);
		}

		void savetxt()
		{
#ifdef LB_DEBUG
			clArray3Dd copy;
			make_copy(copy);
			copy.savetxt(FileNameDistr);
#endif //LB_DEBUG
		};
		void savebin()
		{
			clArray3Dd copy;
			make_copy(copy);

			char FileNameBin[80];
			FILE *stream;
#ifdef SAVE_BIN_IN_LOAD
			sprintf(FileNameBin, "load" SLASH "%s.bin", FileName);
#else //SAVE_BIN_IN_LOAD
			sprintf(FileNameBin, "%s.bin", FileName);
#endif //SAVE_BIN_IN_LOAD

			if(ind == 0)
				stream = fopen(FileNameBin, "w+b");
			else
				stream = fopen(FileNameBin, "a+b");

			BOOL ret = copy.savebin(stream);
			if(ret == TRUE && stream != NULL)
				fclose(stream);
			else
				printf("\nProblem with output to file: %s\n", FileNameBin);


		};
	};

	class clF3Dn: public clF3D
	{
	protected:
		clArray3Dd na;

	public:
		void ini(char* Name, int n, clVector3Di size)
		{
			clF3D::ini(Name, n, size);
			na.zeros(nX, nY, nZ);
			na.fill(c.lb.F0[ind]);
			copy();
		}
		void ini()
		{
			clF3D::ini();
			na.zeros(nX, nY, nZ);
			na.fill(c.lb.F0[ind]);
			copy();
		}
		void copy()
		{
			na.copy(*this);
		}

		void change()
		{
			na.change(*this);
		}

		double& nc(clVector3Di ii)
		{
			return nc(ii.x, ii.y, ii.z);
		}
		double& nc(int i, int j, int m)
		{
			i = MODFAST(i, nX);
			j = MODFAST(j, nY);
			m = MODFAST(m, nZ);
			int is, js, ms;
			is = MODFAST(i + shiftX, nX);
			js = MODFAST(j + shiftY, nY);
			ms = MODFAST(m + shiftZ, nZ);

			return na(is,js,ms);
		}

		double& n(clVector3Di ii)
		{
			return n(ii.x, ii.y, ii.z);
		}

		double& n(int i, int j, int m)
		{
			i = MODFAST(i, nX);
			j = MODFAST(j, nY);
			m = MODFAST(m, nZ);
			int is, js, ms;
			is = MODFAST(i + shiftXn, nX);
			js = MODFAST(j + shiftYn, nY);
			ms = MODFAST(m + shiftZn, nZ);
			return na(is,js,ms);
		}

		double& npx(int i, int j, int z)
		{
			j = MODFAST(j + incY, nYg);
			z = MODFAST(z + incZ, nZg);
			return n(i,j,z);
		}

		double& npy(int i, int j, int z)
		{
			i = MODFAST(i + incX, nXg);
			z = MODFAST(z + incZ, nZg);
			return n(i,j,z);
		}

		double& npz(int i, int j, int z)
		{
			i = MODFAST(i + incX, nXg);
			j = MODFAST(j + incY, nYg);
			return n(i,j,z);
		}
	};
	class clF3Dc: public clF3D
	{
	protected:
		clArray3Dd ca;

	public:
		void ini(char* Name, int n, clVector3Di size)
		{
			clF3D::ini(Name, n, size);
			ca.zeros(nX, nY, nZ);
			ca.fill(c.lb.F0[ind]);
			copy();
		}
		//void ini()
		//{
		//	clF3D::ini();
		//	ca.zeros(nX, nY, nZ);
		//	ca.fill(c.lb.F0[ind]);
		//	copy();
		//}
		void copy()
		{
			ca.copy(*this);
		}

		void change()
		{
			ca.change(*this);
		}

		double& o(clVector3Di ii)
		{
			return o(ii.x, ii.y, ii.z);
		}
		double& o(int i, int j, int m)
		{
			i = MODFAST(i, nX);
			j = MODFAST(j, nY);
			m = MODFAST(m, nZ);
			int is, js, ms;
			is = MODFAST(i + shiftX, nX);
			js = MODFAST(j + shiftY, nY);
			ms = MODFAST(m + shiftZ, nZ);

			return ca(is,js,ms);
		}
	};

	class clF
	{
	public:
#ifndef LB_THIN_WALL
		clF3D F[LB_NUMBER_CONNECTIONS];
#else 
		clF3Dn F[LB_NUMBER_CONNECTIONS];
#endif //LB_THIN_WALL


#ifndef LB_THIN_WALL
		clF3D& operator()(int n)
#else 
		clF3Dn& operator()(int n)
#endif //LB_THIN_WALL
		{
			return F[n];
		};
		virtual void ini()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].ini("lbf", i, c.lb.nn);

		}

		void inishear()
		{
			for(int i0 = 0; i0 < c.lb.nX; i0++)
				for(int j0 = 0; j0 < c.lb.nY; j0++)
					for(int m0 = 0; m0 < c.lb.nZ; m0++)
					{
						double Ux = c.lb.UY0.x + (c.lb.UY1.x - c.lb.UY0.x) * ((double)j0 + 0.5) / (double)c.lb.nY;
						double Uy = 0.;
						for(int i = 0; i < c.lb.nL; i++)
						{
							double U = c.lb.F0[i] * (c.lb.Cx[i] * Ux + c.lb.Cy[i] * Uy) / c.lb.Cs2;
							vlb.F(i)(i0,j0,m0) = c.lb.F0[i] + U;
						}
					}

		}
		void iniflow()
		{
			for(int i0 = 0; i0 < c.lb.nX; i0++)
				for(int j0 = 0; j0 < c.lb.nY; j0++)
					for(int m0 = 0; m0 < c.lb.nZ; m0++)
					{
						if(vlb.M(i0,j0,m0) == LB_SOLID)
							continue;
						clVector3Dd UU = clVector3Dd(0., 0.01, 0.);
						//clVector3Dd UU = -clVector3Dd(i0 - c.lb.nX/2., j0 - c.lb.nY/2., 0.) * 0.0001;
						for(int i = 0; i < c.lb.nL; i++)
						{
							double U = c.lb.F0[i] * (c.lb.C[i] * UU) / c.lb.Cs2;
							vlb.F(i)(i0,j0,m0) = c.lb.F0[i] + U;
						}
					}

		}
		void iniPoiseuille()
		{
			for(int i0 = 0; i0 < c.lb.nX; i0++)
				for(int j0 = 0; j0 < c.lb.nY; j0++)
					for(int m0 = 0; m0 < c.lb.nZ; m0++)
					{
						double a = 0.5 * (double)c.lb.nY;
						double y = ((double)j0 + 0.5) - a;
						double Fx = c.lb.F[vlb.M(i0,j0,m0)].x * c.lb.FM(i0,j0,m0);
						double Ux = c.lb.F[vlb.M(i0,j0,m0)].x * c.lb.FM(i0,j0,m0) * (a * a - y * y) / 2. * 6.;
						double Uy = 0.;
						for(int i = 0; i < c.lb.nL; i++)
						{
							double U = c.lb.F0[i] * (c.lb.Cx[i] * Ux + c.lb.Cy[i] * Uy) / c.lb.Cs2;
							vlb.F(i)(i0,j0,m0) = c.lb.F0[i] + U;
						}
					}

		}
		void update()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].update();
		}

		void propagation()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].propagation();
		}

		void copy()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].copy();
		}

		void change()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].change();
		}

		void save2file()
		{
#ifdef CREATE_BIN_FILES
			savetxt();
			savebin();
#endif //CREATE_BIN_FILES
		};
		void savetxt()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].savetxt();
		};
		void savebin()
		{
			for(int i = 0; i < c.lb.nL; i++)
				F[i].savebin();
		};
	} F;

	class clG: public clF {

	public:
		void ini()
		{


			for(int i = c.lb.nL-1; i >= 0; i--)
				F[i].ini("lbg", i, c.lb.nn);


			//for(int i0 = 0; i0 < c.lb.nX; i0 = i0+4)
			//	for(int j0 = 0; j0 < c.lb.nY; j0 = j0+4)
			//		for(int m0 = 0; m0 < c.lb.nZ; m0 = m0+4)
			//			for(int i = 0; i < c.lb.nL; i++)
			//			{
			//				vlb.G(0)(i0,j0,m0) = -0;
			//			}
//			save2file();
		}

	} G;

	class clN {
	public: 
		clArray3Dv3d d, d2;
//		clArray3Dd x, y, z, x2, y2, z2;
		clArray3Dd a;

		 double operator()(int i, int j, int k) {
			return a(i,j,k);
		 }
		 double operator()(clVector3Di ii) {
			return a(ii);
		 }
		 //double operator()(int i, int j, int k) {
			//return sqrt( x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k));
		 //}

		void ini() {
			a.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			d.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			d2.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//x.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//y.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//z.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//x2.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//y2.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			//z2.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);

			update();
			save2file();
		}

		void update() {

			for (int i=0; i < c.lb.nX; i++) {
				for (int j=0; j < c.lb.nY; j++ ) {
					for (int k=0; k < c.lb.nZ; k++) {

						if ( vlb.M(i,j,k) == LB_SOLID ) continue;

						double F = 0;
						double t = 0;

						for (int n=0; n < c.lb.nL; n++) {
							double tF = vlb.G(n)(i,j,k);
							if ( n > 0 ) t+= tF;
							F += vlb.G(n)(i,j,k);
						}
#ifdef LBN_CORRECTION
						if(F > 1.0)
							F = 1.0;

						if(F < -1.0)
							F = -1.0;
#endif //LBN_CORRECTION
						a(i,j,k) = F;
					}
				}
			}

			updatexyz();
#ifdef LB_DEBUG
			save2file();
#endif //LB_DEBUG
		}

		void updatexyz() {
			for (int i=0; i < c.lb.nX; i++) 
				for (int j=0; j < c.lb.nY; j++) 
					for (int k=0; k < c.lb.nZ; k++) 
						updatexyz(i,j,k);

		if(c.lb.bcWallX0 == TRUE) boundaryX0();
		if(c.lb.bcWallX1 == TRUE) boundaryX1();
		if(c.lb.bcWallY0 == TRUE) boundaryY0();
		if(c.lb.bcWallY1 == TRUE) boundaryY1();
		if(c.lb.bcWallZ0 == TRUE) boundaryZ0();
		if(c.lb.bcWallZ1 == TRUE) boundaryZ1();
		}

		void updatexyz(clVector3Di ii) {
			updatexyz(ii.x,ii.y,ii.z);
		}


		void updatexyz(int i, int j, int k) {

			char M = vlb.M(i,j,k);

			if ( M == LB_SOLID ) return;

			int im = MODFAST(i-1,c.lb.nX);
			int ip = MODFAST(i+1,c.lb.nX);
			int jm = MODFAST(j-1,c.lb.nY);
			int jp = MODFAST(j+1,c.lb.nY);
			int km = MODFAST(k-1,c.lb.nZ);
			int kp = MODFAST(k+1,c.lb.nZ);

			double dxm = 0,
				   dxp = 0,
				   dym = 0,
				   dyp = 0,
				   dzm = 0,
				   dzp = 0;

			if ( vlb.M(im,j,k) == M )	 dxm = a(i,j,k) - a(im,j,k);
			if ( vlb.M(ip,j,k) == M )	 dxp = a(ip,j,k) - a(i,j,k);
			if ( vlb.M(i,jm,k) == M )	 dym = a(i,j,k) - a(i,jm,k);
			if ( vlb.M(i,jp,k) == M )	 dyp = a(i,jp,k) - a(i,j,k);
			if ( vlb.M(i,j,km) == M )	 dzm = a(i,j,k) - a(i,j,km);
			if ( vlb.M(i,j,kp) == M )	 dzp = a(i,j,kp) - a(i,j,k);

			// compute dN/dx
			d(i,j,k).x = (dxm + dxp)/2.;
			d(i,j,k).y = (dym + dyp)/2.;
			d(i,j,k).z = (dzm + dzp)/2.;

			// compute d^2/dN2
			d2(i,j,k).x = dxp-dxm;
			d2(i,j,k).y = dyp-dym;
			d2(i,j,k).z = dzp-dzm;
		}

		double d2abs (int i, int j, int k) {
			return d2(i,j,k).abs();
		}

		void save2file()
		{
			a.savetxt("lbn");
			d.savetxt("lbnd");
			d2.savetxt("lbnd2");
			//y.savetxt("lbny");
			//z.savetxt("lbnz");
		};

		void boundaryX0() 
		{
			int i, j, k;
			for(j = 0; j < c.lb.nY; j++)
				for (k = 0; k < c.lb.nZ; k++) {
					i = 0;
					d(i,j,k).x = 0.;
					d2(i,j,k).x = 0.;
			}
		}
		void boundaryX1() 
		{
			int i, j, k;
			for(j = 0; j < c.lb.nY; j++)
				for (k = 0; k < c.lb.nZ; k++) {
					i = c.nX - 1;
					d(i,j,k).x = 0.;
					d2(i,j,k).x = 0.;
			}
		}
		void boundaryY0() 
		{
			int i, j, k;
			for(i = 0; i < c.lb.nX; i++)
				for (k = 0; k < c.lb.nZ; k++) {
					j = 0;
					d(i,j,k).y = 0.;
					d2(i,j,k).y = 0.;
			}
		}
		void boundaryY1() 
		{
			int i, j, k;
			for(i = 0; i < c.lb.nX; i++)
				for (k = 0; k < c.lb.nZ; k++) {
					j = c.nY - 1;
					d(i,j,k).y = 0.;
					d2(i,j,k).y = 0.;
			}
		}
		void boundaryZ0() 
		{
			int i, j, k;
			for(i = 0; i < c.lb.nX; i++)
				for(j = 0; j < c.lb.nY; j++) {
					k = 0;
					d(i,j,k).z = 0.;
					d2(i,j,k).z = 0.;
			}
		}
		void boundaryZ1() 
		{
			int i, j, k;
			for(i = 0; i < c.lb.nX; i++)
				for(j = 0; j < c.lb.nY; j++) {
					k = c.nZ - 1;
					d(i,j,k).z = 0.;
					d2(i,j,k).z = 0.;
			}
		}


	} N;
	

	class clNM {
	public:
		double S, sum, avg;
		int count;

		void updateoutput()
		{
			sum = 0;
			avg = 0;
			count = 0;
			S = 0;
			
#ifdef BINARY_FLUID
			for (int i=0; i < c.lb.nX; i++) {
				for (int j=0; j < c.lb.nY; j++ ) {
					for (int k=0; k < c.lb.nZ; k++) {
						if ( vlb.M(i,j,k) != LB_SOLID )
						{
							sum = sum + vlb.N.a(i,j,k);
							count = count + 1;
						}
					}
				}
			}
			avg = sum/count;
			sum = 0;
			for (int i=0; i < c.lb.nX; i++) {
				for (int j=0; j < c.lb.nY; j++ ) {
					for (int k=0; k < c.lb.nZ; k++) {
						if ( vlb.M(i,j,k) != LB_SOLID )
						{
							sum = sum + (vlb.N.a(i,j,k) - avg) * (vlb.N.a(i,j,k) - avg);
						}
					}
				}
			}
			S = sqrt(sum/count);
#endif // BINARY_FLUID
		}

	} NM;

	class clRo: public clArray3Dd
	{
	public:
		clArray1Di counter;
		clArray1Dd avr, mass;
		clArray3Dd cor;

		//Modified. Assume 1 type of fluid
		double max, avgabs;
		//Modified
		void ini()
		{
			cor.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			avr.zeros(c.lb.nF);
			mass.zeros(c.lb.nF);
			counter.zeros(c.lb.nF);
			clArray3Dd::zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			//update();
			//for(int i = 0; i < c.lb.nX; i++)
			//	for(int j = 0; j < c.lb.nY; j++)
			//		for(int m = 0; m < c.lb.nZ; m++)
			//			if(vlb.M(i,j,m) != LB_SOLID)
			//				update(clVector3Di(i, j, m));

			//Modified
			max = 0.0;
			avgabs = 0.0;
			//Modified
		}

		void update()
		{
//#ifndef LB_NO_MASS_CORRECTION
//			updateoutput();
//#endif LB_NO_MASS_CORRECTION
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						if(vlb.M(i,j,m) != LB_SOLID)
							update(clVector3Di(i, j, m));
		}

		void update(clVector3Di ii)
		{
			//if(vlb.M(ii) == LB_SOLID)
			//	return;
			double dF = 0.;
#ifndef LB_NO_MASS_CORRECTION
			if(vlb.M.OS(ii) == 0){
				dF = (cor(ii) - avr(vlb.M(ii)) / (double)(c.t.CurrentSaveTimeStep)) / (double)c.lb.nL;
			}
//			dF = (cor(ii) - avr(vlb.M(ii))) / (double)c.lb.nL;
#endif //LB_NO_MASS_CORRECTION
			cor(ii) = 0.;

			double F = 0.;
			for(int k = 0; k < c.lb.nL; k++)
				F += (vlb.F(k)(ii) += dF);

			el(ii) = F * c.lb.Ro[vlb.M(ii)];

		}

		double d(int i, int j, int k)
		{
			return el(i,j,k) - c.lb.Ro[vlb.M(i,j,k)];
		}

		void save2file()
		{
			clArray3Dd copy(c.lb.nX, c.lb.nY, c.lb.nZ);
			for(int i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						copy(i,j,m) = d(i,j,m);
			copy.savetxt("lbro");
			//a.savetxt("lbro");
		}

		void updateoutput()
		{
			mass.zeros();
			counter.zeros();
			//Modified
			max = 0.0;
			avgabs = 0.0;
			//Modified
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int k = 0; k < c.lb.nZ; k++)
					{
						char map = vlb.M(i,j,k);
						if(map != LB_SOLID)
						{
#pragma omp atomic
							mass(map) += d(i,j,k);
#pragma omp atomic

							counter(map)++;
						}
						//Modified
						if(max < fabs(d(i,j,k)) )
							max = fabs(d(i,j,k));
						
						avgabs += fabs(d(i,j,k));
						//Modified
					}

			for(int i = 1; i < c.lb.nF; i++)
				if(counter(i) > 0)
				{
					avr(i) = mass(i) / (double)counter(i);
					mass(i) += c.lb.Ro[i] * (double)counter(i);
				}

				//Modified
				avgabs /= counter(1);
				//Modified
		}
	} Ro;

	class clJ: public clArray3Dv3d
	{
	public:
		clArray1Dv3d avr, tot;
		clArray1Di counter;
		clArray2Dv3d add, rem, cor;
		//clArray2Di Nadd, Nrem;
		clArray2Di Ncor;
		double invel;
		//clArray3Dv3d x;

		virtual void ini()
		{
			add.zeros(c.ls.nO, c.lb.nF);
			rem.zeros(c.ls.nO, c.lb.nF);
			cor.zeros(c.ls.nO, c.lb.nF);
			//Nadd.zeros(c.ls.nO, c.lb.nF);
			//Nrem.zeros(c.ls.nO, c.lb.nF);
			Ncor.zeros(c.ls.nO, c.lb.nF);
			avr.zeros(c.lb.nF);
			tot.zeros(c.lb.nF);
			counter.zeros(c.lb.nF);
			clArray3Dv3d::zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			//update();					
		}

		void update()
		{
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						if(vlb.M(i,j,m) != LB_SOLID)
							update(clVector3Di(i, j, m));
			updatecorrectoin();
		}

		void update(clVector3Di ii)
		{
			clVector3Dd J;
			for(int k = 1; k < c.lb.nL; k++)
			{
				double F = vlb.F(k)(ii);
				J += c.lb.C[k] * F;
			}
			el(ii) = J * c.lb.Ro[vlb.M(ii)];
		}
		void updateinvel()
		{
#ifdef OSC_FLOW
			if(MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD) <= INHALE)
				invel = ::c.t.oscillator1();
			else
				invel = ::c.t.oscillator2();
#else
#ifdef INVEL_BC
			invel = INLET_VEL_MAG;
#else
			invel = 1.;
#endif // INVEL_BC
#endif // OSC_FLOW
		}

		void updatecorrectoin()
		{
#ifdef LB_NO_FLUX_CORRECTION
			return;
#endif //LB_NO_FLUX_CORRECTION
			for(int i = 0; i < c.ls.nO; i++)
				for(int j = 0; j < c.lb.nF; j++)
				{
					//double nc = fabs((double)(Nrem(i,j) - Nadd(i,j)));
					double nc = fabs((double)Ncor(i,j));
					if(nc == 0.) 
						nc = 1.;

					cor(i,j) = (add(i,j) - rem(i,j)) / nc;
				}
		}

		void cor_add(clVector3Dd v, int no, char map)
		{
#ifdef LB_NO_FLUX_CORRECTION
			return;
#endif //LB_NO_FLUX_CORRECTION
			add(no, map) += v - vls.avr.V(no) * c.lb.Ro[map];
			//Nadd(no, map)++;
			Ncor(no, map)++;
		}

		void cor_rem(clVector3Dd v, int no, char map)
		{
#ifdef LB_NO_FLUX_CORRECTION
			return;
#endif //LB_NO_FLUX_CORRECTION
			rem(no, map) += v - vls.avr.V(no) * c.lb.Ro[map];
			//Nrem(no, map)++;
			Ncor(no, map)--;
		}
		void updateoutput()
		{
			avr.zeros();
			counter.zeros();
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int k = 0; k < c.lb.nZ; k++)
					{
						char map = vlb.M(i,j,k);
						if(map != LB_SOLID)
						{
#pragma omp atomic
							avr(map).x += el(i,j,k).x;
#pragma omp atomic
							avr(map).y += el(i,j,k).y;
#pragma omp atomic
							avr(map).z += el(i,j,k).z;
#pragma omp atomic
							counter(map)++;
						}
					}
					for(int i = 1; i < c.lb.nF; i++)
					{
						if(counter(i) > 0)
						{
							tot(i) = avr(i);
							avr(i) = avr(i) / (double)counter(i);
						}
					}
		}

		void save2file()
		{
			clArray3Dv3d copy;
			copy.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						copy(i,j,m) = el(i,j,m) + c.lb.F[vlb.M(i,j,m)] * 0.5 * c.lb.FM(i,j,m);
			copy.savetxt("lbj");
			//copy.savetxt_xyz("lbj");
		};
	} J;

	class clPP: public clArray3Dts3d
	{
	public:
		virtual void ini()
		{
			clArray3Dts3d::zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			//update();
		}

		void update()
		{
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						if(vlb.M(i,j,m) != LB_SOLID)
							update(clVector3Di(i, j, m));
		}

		void update(clVector3Di ii)
		{
			//if(vlb.M(ii) == LB_SOLID)
			//	return;
			clTensorS3Dd PP;
			for(int k = 1; k < c.lb.nL; k++)
			{
				clTensorS3Dd C = clTensorS3Dd(c.lb.C[k]);
				double F = vlb.F(k)(ii);
				PP += C * F;
			}
			el(ii) = PP * c.lb.Ro[vlb.M(ii)];
		}

		void save2file()
		{
#ifdef LB_DEBUG
			savetxt("lbpp");
#endif //LB_DEBUG
		};
	} PP;

	class clM: public clArray3Dc
	{
	public:
		clArray3Dc s, sf;
		clArray2Di S, Sf;
		clArray1Dc Sfdel;
		clArray3Di Bflag, PBC, O, OS;
		clArray3Dv3d D, Dold;

		int nS, nSf;

		void ini()
		{
			clArray3Dc::zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			clArray3Dc::fill(-1);
			s.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			sf.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			Bflag.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			PBC.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			D.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			Dold.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			O.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			OS.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);

			O.fill(-1);

			if(Sf.loadnew("load" SLASH "lbmsf", NUMBER_DIMENSIONS + 1) == FALSE)
				inifill_Sf();

			nSf = Sf.getSize(1);

			Sfdel.zeros(nSf);
			Sfdel.load("load" SLASH "lbmsfd");

			inifill_S();
			//initest_S();

			c.ls.iniBL();

#ifdef LS_RUPTURE
			c.ls.N.iniS();
#endif //LS_RUPTURE

			Sf.save("lbmsf");

			if(load("load" SLASH "lbm") == FALSE)
				//				clArray3Dc::fill(1);
				inifill_map();

			//update();

			//save2file();
		}

		BOOL initest_S()
		{
			clArray3Dc scopy;
			scopy.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			scopy.copy(s);
			fill_s();
			for(int i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						if(sf(i,j,m) == 0 && s(i,j,m) != scopy(i,j,m))
						{
							printf("s-map test failed!!!\n");
							return FALSE;
						}
						return TRUE;
		}

		void inivol()
		{
			int nScount = 0;
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.T(i) & LS_VOLUME_NODE)
					for(int j = 0; j < c.ls.nL; j++)
						if(c.ls.N(i,j) >= 0 && c.ls.T(c.ls.N(i,j)) & LS_VOLUME_NODE)
							for(int m = j + 1; m < c.ls.nL; m++)
								if(c.ls.N(i,m) >= 0 && c.ls.T(c.ls.N(i,m)) & LS_VOLUME_NODE)
									for(int k = m + 1; k < c.ls.nL; k++)
									{
										if(c.ls.N(i,k) >= 0 && c.ls.T(c.ls.N(i,k)) & LS_VOLUME_NODE)
										{
											int n0 = i, n1 = c.ls.N(i,j), n2 = c.ls.N(i,m), n3 = c.ls.N(i,k);
											if(n0 < n1 && n0 < n2 && n0 < n3 && test_connection(n1, n2, n3) && test_volume(n0, n1, n2, n3))
												//if(test_connection(n1, n2, n3) && test_volume(n0, n1, n2, n3))
												//if(test_connection(n1, n2, n3))
											{
												if(nScount >= nS)
												{
													nS += c.ls.nN;
													S.reallocate(nS, NUMBER_DIMENSIONS + 1);
												}

												S(nScount, 0) = n0;
												S(nScount, 1) = n1;
												S(nScount, 2) = n2;
												S(nScount, 3) = n3;
												nScount++;
											}
										}
									}

			}
			nS = nScount;
			S.reallocate(nS, NUMBER_DIMENSIONS + 1);
			//S.savetxt("lbmsf");
		}
		BOOL test_connection(int n0, int n1, int n2)
		{
			return test_connection(n0, n1) && test_connection(n0, n2) && test_connection(n1, n2);
		}
		BOOL test_connection(int n0, int n1)
		{
			for(int i = 0; i < c.ls.nL; i++)
			{
				if(c.ls.N(n0,i) == n1)
					return TRUE;
			}
			return FALSE;
		}

		BOOL test_volume(int n0, int n1, int n2, int n3)
		{
			clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
			if(volume(v1-v0, v2-v0, v3-v0) < c.eps)
				return FALSE;
			return TRUE;
		}

		clVector3Di findCmax()
		{
			clVector3Dd Cmax;
			for(int n = 0; n < c.ls.nN; n++)
				Cmax = Cmax.max(vls.C(n));

			return (clVector3Di)Cmax.ceil() * 4 + 1.;
		}
		clVector3Di findCmin()
		{
			clVector3Dd Cmin = c.lb.nn;
			for(int n = 0; n < c.ls.nN; n++)
				Cmin = Cmin.min(vls.C(n));

			return (clVector3Di)Cmin.floor() * 4;
		}
		void inifill_Sf()
		{
			nS = 0;
			nSf = 0;

			inivol();

			clVector3Di Cmax = findCmax(), Cmin = findCmin();

			clVector3Di nCs = Cmax - Cmin;
			if(nCs.x < 0 || nCs.y < 0 || nCs.z < 0)
				return;

			clArray3Dc s4, s4_test;
			s4.zeros(nCs);
			s4_test.zeros(nCs);

			clArray1Dc s_remove;
			s_remove.zeros(nS);
			int s_rem_count = 0;

			for(int n = 0; n < nS; n++)
			{
				int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
				if(test_boundary_node3(n0, n1, n2, n3) != TRUE)
					continue;

				clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
				clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3) * 4. - Cmin).ceil()).max(clVector3Di(0,0,0));
				clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3) * 4. - Cmin).floor() + 1).min(nCs);

				clVector3Di ii;
				for(ii.x = imin.x; ii.x < imax.x; ii.x++)
					for(ii.y = imin.y; ii.y < imax.y; ii.y++)
						for(ii.z = imin.z; ii.z < imax.z; ii.z++)
						{
							////if(test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
							////	s4_test(ii) = c.ls.O(n0) + 1;

							if(s4(ii) != 0 && s_remove(n) == FALSE 
								&& test_inside_noboundary((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
							{
								s_remove(n) = TRUE;
								s_rem_count++;
								break;
							}
						}

				if(s_remove(n) == TRUE)
					continue;

				for(ii.x = imin.x; ii.x < imax.x; ii.x++)
					for(ii.y = imin.y; ii.y < imax.y; ii.y++)
						for(ii.z = imin.z; ii.z < imax.z; ii.z++)
							if(test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
								s4(ii) = c.ls.O(n0) + 1;
			}

			for(int n = 0; n < nS; n++)
			{
				int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
				if(test_boundary_node3(n0, n1, n2, n3) == TRUE)
					continue;

				clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
				clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3) * 4. - Cmin).ceil()).max(clVector3Di(0,0,0));
				clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3) * 4. - Cmin).floor() + 1).min(nCs);

				clVector3Di ii;
				for(ii.x = imin.x; ii.x < imax.x; ii.x++)
					for(ii.y = imin.y; ii.y < imax.y; ii.y++)
						for(ii.z = imin.z; ii.z < imax.z; ii.z++)
						{
							////if(test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
							////	s4_test(ii) = c.ls.O(n0) + 1;

							if(s4(ii) != 0 && s_remove(n) == FALSE 
								&& test_inside_noboundary((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
							{
								s_remove(n) = TRUE;
								s_rem_count++;
								break;
							}
						}

				if(s_remove(n) == TRUE)
					continue;

				for(ii.x = imin.x; ii.x < imax.x; ii.x++)
					for(ii.y = imin.y; ii.y < imax.y; ii.y++)
						for(ii.z = imin.z; ii.z < imax.z; ii.z++)
							if(test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
								s4(ii) = c.ls.O(n0) + 1;
			}

			//for(int n = 0; n < nS; n++)
			//{
			//	int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
			//	clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
			//	clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3) * 4. - Cmin).ceil()).max(clVector3Di(0,0,0));
			//	clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3) * 4. - Cmin).floor() + 1).min(nCs);
			//	clVector3Di ii;
			//	for(ii.x = imin.x; ii.x < imax.x; ii.x++)
			//		for(ii.y = imin.y; ii.y < imax.y; ii.y++)
			//			for(ii.z = imin.z; ii.z < imax.z; ii.z++)
			//				if(test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
			//					s4_test(ii) = c.ls.O(n0) + 1;
			//}

			////for(int n = 0; n < nS; n++)
			////{
			////	int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
			////	if(test_boundary_node3(n0, n1, n2, n3) != TRUE)
			////		continue;

			////	clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
			////	clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3) * 4. - Cmin).ceil()).max(clVector3Di(0,0,0));
			////	clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3) * 4. - Cmin).floor() + 1).min(nCs);
			////	clVector3Di ii;
			////	for(ii.x = imin.x; ii.x < imax.x; ii.x++)
			////		for(ii.y = imin.y; ii.y < imax.y; ii.y++)
			////			for(ii.z = imin.z; ii.z < imax.z; ii.z++)
			////				if((s4(ii) != s4_test(ii)) 
			////					&& test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
			////				{
			////					s4(ii) = c.ls.O(n0) + 1;
			////					if(s_remove(n) == TRUE)
			////					{
			////						s_remove(n) = FALSE;
			////						s_rem_count--;
			////					}
			////				}
			////}
			////for(int n = 0; n < nS; n++)
			////{
			////	int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
			////	if(test_boundary_node3(n0, n1, n2, n3) == TRUE)
			////		continue;

			////	clVector3Dd v0 = vls.C(n0),  v1 = vls.C(n1), v2 = vls.C(n2), v3 = vls.C(n3);
			////	clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3) * 4. - Cmin).ceil()).max(clVector3Di(0,0,0));
			////	clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3) * 4. - Cmin).floor() + 1).min(nCs);
			////	clVector3Di ii;
			////	for(ii.x = imin.x; ii.x < imax.x; ii.x++)
			////		for(ii.y = imin.y; ii.y < imax.y; ii.y++)
			////			for(ii.z = imin.z; ii.z < imax.z; ii.z++)
			////				if((s4(ii) != s4_test(ii)) 
			////					&& test_inside((clVector3Dd)(ii + Cmin)/4., v0, v1, v2, v3) == TRUE)
			////				{
			////					s4(ii) = c.ls.O(n0) + 1;
			////					if(s_remove(n) == TRUE)
			////					{
			////						s_remove(n) = FALSE;
			////						s_rem_count--;
			////					}
			////				}
			////}

			nSf = nS - s_rem_count;
			Sf.zeros(nSf, NUMBER_DIMENSIONS + 1);
			int j = 0;
			for(int i = 0; i < nS; i++)
			{
				if(s_remove(i) == FALSE)
				{
					Sf(j,0) = S(i,0);
					Sf(j,1) = S(i,1);
					Sf(j,2) = S(i,2);
					Sf(j,3) = S(i,3);
					j++;
				}
			}

			nS = nSf;
			S.reallocate(nS, NUMBER_DIMENSIONS + 1);
			S.copy(Sf);

			fill_s();
			fill_sf();

			//Sf.savetxt("lbmsfull");
			//s4.savetxt("lbms4");
			//s4_test.savetxt("lbmstest4");
			//s.savetxt("lbmsini");
		}
		void inifill_S()
		{
			clArray1Dc s_remove;
			s_remove.zeros(nSf);
			int s_rem_count = 0;

			for(int i = 0; i < nSf; i++)
			{
				//if(s_remove(i) == FALSE && (test_stationary(i) == TRUE || test_boundary(i) == FALSE))
				//{
				//	s_remove(i) = TRUE;
				//	s_rem_count++;
				//}
				if(Sfdel(i) == TRUE || (s_remove(i) == FALSE && test_boundary(i) == FALSE))
				{
					s_remove(i) = TRUE;
					s_rem_count++;
				}
			}

			nS = nSf - s_rem_count;
			S.reallocate(nS, NUMBER_DIMENSIONS + 1);
			int j = 0;
			for(int i = 0; i < nSf; i++)
			{
				if(s_remove(i) == FALSE)
				{
					S(j,0) = Sf(i,0);
					S(j,1) = Sf(i,1);
					S(j,2) = Sf(i,2);
					S(j,3) = Sf(i,3);
					j++;
				}
			}

			//S.savetxt("lbmsshort");
		}
		void fill_s()
		{
			s.zeros();
			int n;
#ifdef LS_FILL_USE_OPENMP	// not LS_USE_OPENMP: fill_s accumulates flux corrections into shared LB cells (race)
#pragma omp parallel for default(shared) private(n)
#endif //LS_FILL_USE_OPENMP
			for(n = 0; n < nS; n++)
			{
				int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
				clVector3Di p0 = vls.C.p(n0), p1 = vls.C.p(n1), p2 = vls.C.p(n2), p3 = vls.C.p(n3);
				if(p0 == p1 && p1 == p2 && p2 == p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					continue;
				}
				if(p0 != p1 && p1 == p2 && p2 == p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n1);
					continue;
				}
				if(p0 == p1 && p1 != p2 && p2 == p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n2);
					continue;
				}
				if(p0 == p1 && p1 == p2 && p2 != p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 != p1 && p1 != p2 && p2 == p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n1);
					fill_s(n0, n1, n2, n3, n2);
					continue;
				}
				if(p0 != p1 && p1 == p2 && p2 != p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n1);
					fill_s(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 == p1 && p1 != p2 && p2 != p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n2);
					fill_s(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 != p1 && p1 != p2 && p2 != p3)
				{
					fill_s(n0, n1, n2, n3, n0);
					fill_s(n0, n1, n2, n3, n1);
					fill_s(n0, n1, n2, n3, n2);
					fill_s(n0, n1, n2, n3, n3);
					continue;
				}

			}
		}
		void fill_s(int n0, int n1, int n2, int n3, int p)
		{
			clVector3Dd c0 = vls.C.s(n0,p),  c1 = vls.C.s(n1,p), c2 = vls.C.s(n2,p), c3 = vls.C.s(n3,p);
			clVector3Di imin = clVector3Di(clVector3Dd(min(c0, c1, c2, c3)).ceil()).max(clVector3Di(0,0,0));
			clVector3Di imax = clVector3Di(clVector3Dd(max(c0, c1, c2, c3)).floor() + 1).min(c.lb.nn);
			clVector3Di ii;
			for(ii.x = imin.x; ii.x < imax.x; ii.x++)
				for(ii.y = imin.y; ii.y < imax.y; ii.y++)
					for(ii.z = imin.z; ii.z < imax.z; ii.z++)
						if(test_inside((clVector3Dd)ii, c0, c1, c2, c3) == TRUE)
						{
#ifndef LB_NO_FLUX_CORRECTION
							if(el(ii) > LB_SOLID)
								vlb.J.cor_rem(vlb.J(ii), c.ls.O(n0), el(ii));
#endif //LB_NO_FLUX_CORRECTION

							s(ii) = c.ls.O(n0) + 1;
							el(ii) = LB_SOLID;
						}

		}
		void fill_sf()
		{
			sf.zeros();
			int n;
#ifdef LS_FILL_USE_OPENMP	// not LS_USE_OPENMP: fill_s accumulates flux corrections into shared LB cells (race)
#pragma omp parallel for default(shared) private(n)
#endif //LS_FILL_USE_OPENMP
			for(n = 0; n < nS; n++)
			{
				if(test_stationary(n) == FALSE)
					continue;

				int n0 = Sf(n,0), n1 = Sf(n,1), n2 = Sf(n,2), n3 = Sf(n,3); 
				clVector3Dd c0 = vls.C(n0),  c1 = vls.C(n1), c2 = vls.C(n2), c3 = vls.C(n3);
				clVector3Di imin = clVector3Di(clVector3Dd(min(c0, c1, c2, c3)).ceil()).max(clVector3Di(0,0,0));
				clVector3Di imax = clVector3Di(clVector3Dd(max(c0, c1, c2, c3)).floor() + 1).min(c.lb.nn);
				clVector3Di ii;
				for(ii.x = imin.x; ii.x < imax.x; ii.x++)
					for(ii.y = imin.y; ii.y < imax.y; ii.y++)
						for(ii.z = imin.z; ii.z < imax.z; ii.z++)
							if(test_inside(clVector3Dd(ii), c0, c1, c2, c3) == TRUE)
								sf(ii) = c.ls.O(n0) + 1;
			}
			//s_stat.savetxt("lbmsstat");
		}

		void inifill_map()
		{
			if(c.ls.nN == 0)
				fill(1);

			for(int n = 0; n < c.ls.nBL; n++)
			{
				int n0 = c.ls.BL(n,0), n1 = c.ls.BL(n,1), n2 = c.ls.BL(n,2);
				clVector3Dd c0 = vls.C.s(n0, n0), c1 = vls.C.s(n1, n0), c2 = vls.C.s(n2, n0);
				clVector3Dd cn = clVector3Dd((c1 - c0) ^ (c2 - c0)).n();
				//clVector3Dd NL = c.ls.directionBL(n0) + c.ls.directionBL(n1) + c.ls.directionBL(n2); //fix for capsules
				//if(cn * NL > 0.)
				//	fill_m_bn(cn + (c0 + c1 + c2) / 3., c.ls.BF(c.ls.BNind(n0)) + 1);
				//else
					fill_m_bn(cn + (c0 + c1 + c2) / 3., c.ls.BF(c.ls.BNind(n0)));
			}
			//save2file();
			BOOL flRepeate, flChange;
			do
			{
				flRepeate = FALSE;
				flChange = FALSE;
				for(int i = 0; i < c.lb.nX; i++)
					for(int j = 0; j < c.lb.nY; j++)
						for(int z = 0; z < c.lb.nZ; z++)
						{
							if(el(i,j,z) == -1)
							{
								flRepeate = TRUE;
								if(fill_m_around(i,j,z) == TRUE)
									flChange = TRUE;
							}
						}
						for(int i = c.lb.nX-1; i >= 0; i--)
							for(int j = c.lb.nY-1; j >= 0; j--)
								for(int z = c.lb.nZ-1; z >= 0; z--)
								{
									if(el(i,j,z) == -1)
									{
										flRepeate = TRUE;
										if(fill_m_around(i,j,z) == TRUE)
											flChange = TRUE;
									}
								}
			}
			while(flRepeate == TRUE && flChange == TRUE);

			//save2file();

			for(int i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int z = 0; z < c.lb.nZ; z++)
						if(el(i,j,z) == -1)
							el(i,j,z) = LB_SOLID;
		}

		void fill_m_bn(clVector3Dd v, char m_fluid)
		{
			clVector3Di imin = clVector3Di(v.floor()).max(clVector3Di(0,0,0));
			clVector3Di imax = clVector3Di(v.ceil()).min(c.lb.nn-1);
			clVector3Di ii;
			for(ii.x = imin.x; ii.x <= imax.x; ii.x++)
				for(ii.y = imin.y; ii.y <= imax.y; ii.y++)
					for(ii.z = imin.z; ii.z <= imax.z; ii.z++)
						if(s(ii) == FALSE && el(ii) == -1)
							el(ii) = m_fluid;
		}

		BOOL fill_m_around(int x, int y, int z)
		{
			char m_fluid = find_m_fluid(x, y, z);

			if(m_fluid <= 0)
				return FALSE;

			el(x, y, z) = m_fluid;
			return TRUE;

			//BOOL flChange = FALSE;
			//for(int i = MAX(x-1, 0); i < MIN(x+2, c.lb.nX); i++)
			//	for(int j = MAX(y-1, 0); j < MIN(y+2, c.lb.nY); j++)
			//		for(int k = MAX(z-1, 0); k < MIN(z+2, c.lb.nZ); k++)
			//			if(s(i,j,k) == FALSE && el(i,j,k) == -1)
			//			{
			//				el(i,j,k) = m_fluid;
			//				flChange = TRUE;
			//			}
			//return flChange;
		}

		char find_m_fluid(int x, int y, int z)
		{
			int xm = MAX(x-1, 0), xp = MIN(x+1, c.lb.nX-1);
			if(el(xm,y,z) > LB_SOLID)
				return el(xm,y,z);
			if(el(xp,y,z) > LB_SOLID)
				return el(xp,y,z);
			int ym = MAX(y-1, 0), yp = MIN(y+1, c.lb.nY-1);
			if(el(x,ym,z) > LB_SOLID)
				return el(x,ym,z);
			if(el(x,yp,z) > LB_SOLID)
				return el(x,yp,z);
			int zm = MAX(z-1, 0), zp = MIN(z+1, c.lb.nZ-1);
			if(el(x,y,zm) > LB_SOLID)
				return el(x,y,zm);
			if(el(x,y,zp) > LB_SOLID)
				return el(x,y,zp);
			//for(int i = MAX(x-1, 0); i < MIN(x+2, c.lb.nX); i++)
			//	for(int j = MAX(y-1, 0); j < MIN(y+2, c.lb.nY); j++)
			//		for(int k = MAX(z-1, 0); k < MIN(z+2, c.lb.nZ); k++)
			//			if(s(i,j,k) == FALSE && el(i,j,k) != -1)
			//				return el(i,j,k);
			return -1;
		}

		BOOL test_stationary(int n)
		{
			if(c.ls.SO(c.ls.O(Sf(n,0))) == TRUE && c.ls.SO(c.ls.O(Sf(n,1))) == TRUE && 
				c.ls.SO(c.ls.O(Sf(n,2))) == TRUE && c.ls.SO(c.ls.O(Sf(n,3))) == TRUE)
				return TRUE;
			if(c.ls.Ts(Sf(n,0)) == TRUE && c.ls.Ts(Sf(n,1)) == TRUE && 
				c.ls.Ts(Sf(n,2)) == TRUE && c.ls.Ts(Sf(n,3)) == TRUE)
				return TRUE;
			return FALSE;
		}

		BOOL test_boundary(int n)
		{
			if(c.ls.BNind(Sf(n,0)) >= 0 || c.ls.BNind(Sf(n,1)) >= 0 || c.ls.BNind(Sf(n,2)) >= 0 || c.ls.BNind(Sf(n,3)) >= 0)
				return TRUE;
			return FALSE;
		}
		BOOL test_boundary_node3(int n0, int n1, int n2, int n3)
		{
			int count = 0;
			if(c.ls.T(n0) & LS_OUTER_NODE)
				count++;
			if(c.ls.T(n1) & LS_OUTER_NODE)
				count++;
			if(c.ls.T(n2) & LS_OUTER_NODE)
				count++;
			if(c.ls.T(n3) & LS_OUTER_NODE)
				count++;

			if(count >= 3)
				return TRUE;
			return FALSE;
		}
		BOOL test_inside(clVector3Dd vt, clVector3Dd v0, clVector3Dd v1, clVector3Dd v2, clVector3Dd v3)
		{
			clVector3Dd v10 = v1 - v0,  v20 = v2 - v0,  v30 = v3 - v0;
			clVector3Dd vt0 = vt - v0,  vt1 = vt - v1,  vt2 = vt - v2,  vt3 = vt - v3;
			double Stot = volume(v10, v20, v30);
			double S012 = volume(vt0, vt1, vt2), S023 = volume(vt0, vt2, vt3), S013 = volume(vt0, vt1, vt3), S123 = volume(vt1, vt2, vt3);
			double Stest = S012 + S023 + S013 + S123;
			if(fabs(Stot - Stest) < c.eps)
				return TRUE;
			return FALSE;
		}
		BOOL test_inside_noboundary(clVector3Dd vt, clVector3Dd v0, clVector3Dd v1, clVector3Dd v2, clVector3Dd v3)
		{
			clVector3Dd v10 = v1 - v0,  v20 = v2 - v0,  v30 = v3 - v0;
			clVector3Dd vt0 = vt - v0,  vt1 = vt - v1,  vt2 = vt - v2,  vt3 = vt - v3;
			double S012 = volume(vt0, vt1, vt2), S023 = volume(vt0, vt2, vt3), S013 = volume(vt0, vt1, vt3), S123 = volume(vt1, vt2, vt3);
			if(fabs(S012) < c.eps || fabs(S023) < c.eps || fabs(S013) < c.eps || fabs(S123) < c.eps)
				return FALSE;

			double Stot = volume(v10, v20, v30);
			double Stest = S012 + S023 + S013 + S123;
			if(fabs(Stot - Stest) < c.eps)
				return TRUE;
			return FALSE;
		}
		double volume(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
		{
			return fabs(v0 * (v1 ^ v2) / 6.);
		}
		clVector3Dd max(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2, clVector3Dd v3)
		{
			return clVector3Dd(max(v0.x,v1.x,v2.x,v3.x), max(v0.y,v1.y,v2.y,v3.y), max(v0.z,v1.z,v2.z,v3.z));
		}
		clVector3Dd min(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2, clVector3Dd v3)
		{
			return clVector3Dd(min(v0.x,v1.x,v2.x,v3.x), min(v0.y,v1.y,v2.y,v3.y), min(v0.z,v1.z,v2.z,v3.z));
		}
		double max(double n0, double n1, double n2, double n3)
		{
			if(n0 > n1) n1 = n0;
			if(n1 > n2) n2 = n1;
			if(n2 > n3) return n2;
			return n3;
		}
		double min(double n0, double n1, double n2, double n3)
		{
			if(n0 < n1) n1 = n0;
			if(n1 < n2) n2 = n1;
			if(n2 < n3) return n2;
			return n3;
		}
		void update()
		{
			fill_s();

			for(int z = 0; z < c.lb.nZ; z++)
				for(int j = 0; j < LB_WALL_SHIFT; j++)
					for(int i = 0; i < c.lb.nX/2; i++)
						if(vlb.M.s(i,j,z) == 0)
						{
							vlb.M.s(i,j,z) = 1;
							vlb.M(i,j,z) = LB_SOLID;
						}
						else
							break;

			for(int z = 0; z < c.lb.nZ; z++)
				for(int j = 0; j < LB_WALL_SHIFT; j++)
					for(int i = c.lb.nX-1; i >= c.lb.nX/2; i--)
						if(vlb.M.s(i,j,z) == 0)
						{
							vlb.M.s(i,j,z) = 1;
							vlb.M(i,j,z) = LB_SOLID;
						}
						else
							break;

			for(int z = 0; z < c.lb.nZ; z++)
				for(int j = c.lb.nY - 1; j >= c.lb.nY - LB_WALL_SHIFT; j--)
					for(int i = 0; i < c.lb.nX/2; i++)
						if(vlb.M.s(i,j,z) == 0)
						{
							vlb.M.s(i,j,z) = 1;
							vlb.M(i,j,z) = LB_SOLID;
						}
						else
							break;

			for(int z = 0; z < c.lb.nZ; z++)
				for(int j = c.lb.nY - 1; j >= c.lb.nY - LB_WALL_SHIFT; j--)
					for(int i = c.lb.nX-1; i >= c.lb.nX/2; i--)
						if(vlb.M.s(i,j,z) == 0)
						{
							vlb.M.s(i,j,z) = 1;
							vlb.M(i,j,z) = LB_SOLID;
						}
						else
							break;


			Bflag.zeros();
			PBC.zeros();
			Dold.change(D);
			D.zeros();
			O.fill(-1);
			OS.zeros();
#ifdef LB_DEBUG
			save2file();
#endif //LB_DEBUG
		}

		void test()
		{
			int i;
			//for(i = 0; i < c.lb.nX; i++)
			//	for(int j = 0; j < c.lb.nY; j++)
			//		for(int z = 0; z < c.lb.nZ; z++)
			//		{
			//			clVector3Di ii = clVector3Di(i, j, z);
			//			if(el(ii) != LB_SOLID)
			//				for(int dir = 0; dir < c.lb.nL; dir++)
			//				{
			//					//clVector3Di ii1 = MOD(ii+(clVector3Di)c.lb.C[dir], c.lb.nn);
			//					clVector3Di ii1 = ii+(clVector3Di)c.lb.C[dir];
			//					ii1 = MOD(ii1, c.lb.nn);
			//					//if(ii1 != MOD(ii1, c.lb.nn))
			//					//	continue;
			//					int dirp = c.lb.rev[dir];
			//					if(((Bflag(ii) & c.lb.LF[dir]) && !(Bflag(ii1) & c.lb.LF[dirp])) || (!(Bflag(ii) & c.lb.LF[dir]) && (Bflag(ii1) & c.lb.LF[dirp])))
			//					{
			//						printf("Error in BC!!! n0 %d %d %d n1 %d %d %d dir0 %d dir1 %d\n", ii.x, ii.y, ii.z, ii1.x, ii1.y, ii1.z, dir, dirp);
			//					}
			//					//if(((Bflag(ii) & c.lb.LF[dir]) && (Bflag(ii1) & c.lb.LF[dirp])))
			//					//{
			//					//	printf("BC is OK!!! n0 %d %d %d n1 %d %d %d dir0 %d dir1 %d\n", ii.x, ii.y, ii.z, ii1.x, ii1.y, ii1.z, dir, dirp);
			//					//}
			//				}
			//		}



#ifdef LB_TEST_BC
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int z = 0; z < c.lb.nZ; z++)
					{
						clVector3Di ii = clVector3Di(i, j, z);
						if(el(ii) != LB_SOLID)
							for(int dir = 0; dir < c.lb.nL; dir++)
							{
								//clVector3Di ii1 = MOD(ii+(clVector3Di)c.lb.C[dir], c.lb.nn);
								clVector3Di ii1 = ii+(clVector3Di)c.lb.C[dir];
								/*if(ii1 != MOD(ii1, c.lb.nn))
									continue;*/
								if((i == 0 && c.lb.C[dir].x == -1) || (i == c.lb.nX-1 && c.lb.C[dir].x == 1))
									continue;
								if((j == 0 && c.lb.C[dir].y == -1) || (j == c.lb.nY-1 && c.lb.C[dir].y == 1))
									continue;
								if((z == 0 && c.lb.C[dir].z == -1) || (z == c.lb.nZ-1 && c.lb.C[dir].z == 1))
									continue;
								if(el(ii1) == LB_SOLID && !(Bflag(ii) & c.lb.LF[dir]))
								{
									double F = vlb.F(dir)(ii);
									double Ro = vlb.Ro(ii); 
									double UC = vlb.J(ii) * c.lb.C[dir] / Ro;
									int dirp = c.lb.rev[dir];
									//double Fn = F - 2. * c.lb.F0[dir] * UC / c.lb.Cs2;
									//if(c.lb.Cy[dir] == -1)
									//	printf("\nError Fn in test\n");
									//if(Fn < c.eps)
									//	printf("\nError Fn in test\n");
									vlb.F(dirp).n(ii) = F - 2. * c.lb.F0[dir] * UC / c.lb.Cs2;
									//if(ii.x == 0 && ii.y == 4 && ii.z == 0 && dir == 4)
									//{
									//	double Fn0 = F - 2. * c.lb.F0[dir] * UC / c.lb.Cs2;
									//	double Fn1 = vlb.F(dirp).n(ii);
									//	double Fn2 = vlb.F(dirp)(ii.x, ii.y-1, ii.z);
									//	Fn1 = Fn1;
									//}

#ifdef LB_DEBUG
									vlb.B.D[dir](ii) = 0.5;
									vlb.B.Ux[dir](ii) = vlb.J(ii).x / Ro;
									vlb.B.Uy[dir](ii) = vlb.J(ii).y / Ro;
									vlb.B.Uz[dir](ii) = vlb.J(ii).z / Ro;
#endif //LB_DEBUG
								}
							}
					}
#endif //LB_TEST_BC
		}

		void remove_s(int d0, int d1)
		{
			if(d0 == d1)
				return;

			int bf = 1, bs = 0;
			if(c.ls.BNind(d0) >= 0)
			{
				bf = c.ls.BF(c.ls.BNind(d0));
				bs = c.ls.BS(c.ls.BNind(d0));
			}
			else if(c.ls.BNind(d1) >= 0)
			{
				bf = c.ls.BF(c.ls.BNind(d1));
				bs = c.ls.BS(c.ls.BNind(d1));
			}

			for(int n = 0; n < nSf; n++)
			{
				if(Sfdel(n) == TRUE)
					continue;

				int n0 = Sf(n,0), n1 = Sf(n,1), n2 = Sf(n,2), n3 = Sf(n,3);
				int s0 = -1, s1;
				if((d0 == n0 || d0 == n1 || d0 == n2 || d0 == n3) 
					&& (d1 == n0 || d1 == n1 || d1 == n2 || d1 == n3))
				{
					if(d0 != n0 && d1 != n0)
						s0 = n0;
					if(d0 != n1 && d1 != n1)
					{
						if(s0 == -1)
							s0 = n1;
						else
							s1 = n1;
					}
					if(d0 != n2 && d1 != n2)
					{
						if(s0 == -1)
							s0 = n2;
						else
							s1 = n2;
					}
					if(d0 != n3 && d1 != n3)
						s1 = n3;

					c.ls.addB(s0, bf, bs);
					c.ls.addB(s1, bf, bs);

					add_BL(s0, s1, d0, d1);
					add_BL(s0, s1, d1, d0);

					Sfdel(n) = TRUE;

					//remove_Sfull(n--);
				}
			}
		}
		void add_BL(int n0, int n1, int n2, int d)
		{
			clVector3Dd C0 = vls.C(n0),  C1 = vls.C(n1), C2 = vls.C(n2), Cd = vls.C(d);
			clVector3Dd Cc = Cd - (C0 + C1 + C2) / 3.;
			clVector3Dd Cn = (C1 - C0) ^ (C2 - C0);
			if(Cc * Cn > 0.)
				c.ls.addBL(n0, n1, n2);
			else
				c.ls.addBL(n0, n2, n1);
		}
		//void remove_Sfull(int d)
		//{
		//	for(int n = d + 1; n < nSf; n++)
		//	{
		//		Sf(n-1, 0) = Sf(n, 0);
		//		Sf(n-1, 1) = Sf(n, 1);
		//		Sf(n-1, 2) = Sf(n, 2);
		//		Sf(n-1, 3) = Sf(n, 3);
		//	}
		//	nSf--;
		//	Sf.reallocate(nSf, NUMBER_DIMENSIONS + 1);
		//}
		void fill_f(int n, int f)
		{
			int n0 = S(n,0), n1 = S(n,1), n2 = S(n,2), n3 = S(n,3); 
			clVector3Di p0 = vls.C.p(n0), p1 = vls.C.p(n1), p2 = vls.C.p(n2), p3 = vls.C.p(n3);
			if(p0 == p1 && p1 == p2 && p2 == p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				return;
			}
			if(p0 != p1 && p1 == p2 && p2 == p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n1, f);
				return;
			}
			if(p0 == p1 && p1 != p2 && p2 == p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n2, f);
				return;
			}
			if(p0 == p1 && p1 == p2 && p2 != p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n3, f);
				return;
			}
			if(p0 != p1 && p1 != p2 && p2 == p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n1, f);
				fill_f(n0, n1, n2, n3, n2, f);
				return;
			}
			if(p0 != p1 && p1 == p2 && p2 != p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n1, f);
				fill_f(n0, n1, n2, n3, n3, f);
				return;
			}
			if(p0 == p1 && p1 != p2 && p2 != p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n2, f);
				fill_f(n0, n1, n2, n3, n3, f);
				return;
			}
			if(p0 != p1 && p1 != p2 && p2 != p3)
			{
				fill_f(n0, n1, n2, n3, n0, f);
				fill_f(n0, n1, n2, n3, n1, f);
				fill_f(n0, n1, n2, n3, n2, f);
				fill_f(n0, n1, n2, n3, n3, f);
				return;
			}

		}
		void fill_f(int n0, int n1, int n2, int n3, int p, int map)
		{
			clVector3Dd v0 = vls.C.s(n0,p),  v1 = vls.C.s(n1,p), v2 = vls.C.s(n2,p), v3 = vls.C.s(n3,p);
			clVector3Di imin = clVector3Di(clVector3Dd(min(v0, v1, v2, v3)).ceil()).max(clVector3Di(0,0,0));
			clVector3Di imax = clVector3Di(clVector3Dd(max(v0, v1, v2, v3)).floor() + 1).min(c.lb.nn);
			clVector3Di ii;
			for(ii.x = imin.x; ii.x < imax.x; ii.x++)
				for(ii.y = imin.y; ii.y < imax.y; ii.y++)
					for(ii.z = imin.z; ii.z < imax.z; ii.z++)
						if(test_inside((clVector3Dd)ii, v0, v1, v2, v3) == TRUE)
							fill_f(ii, (vls.V(n0) + vls.V(n1) + vls.V(n2) + vls.V(n3)) / 4., c.ls.O(n0), map);
		}
		void fill_f(clVector3Di ii, clVector3Dd UU, int no, int map)
		{
			double Ro = c.lb.Ro[map];
			for(int i = 0; i < c.lb.nL; i++)
			{
				clVector3Di ii1 = MOD(ii + c.lb.C[i], c.lb.nn);
				clVector3Dd UUc = UU - vlb.J.cor(no, map) / Ro;
				double U = (c.lb.C[i] * UUc) * c.lb.F0[i] / c.lb.Cs2;
				double F = Ro * (c.lb.F0[i] + U) / c.lb.Ro[map];
				vlb.F(i)(ii) = F;
#ifdef LB_THIN_WALL
				vlb.F(i).nc(ii) = F;
#endif //LB_THIN_WALL
			}

			s(ii) = 0;
			el(ii) = map;

			vlb.Ro.update(ii);
			vlb.J.update(ii);
			vlb.J.cor_add(vlb.J(ii), no, map);

		}
		void save2file()
		{
			clArray3Dc::save("lbm");
			s.savetxt("lbms");
			Sf.save("lbmsf");
			Sfdel.save("lbmsfd");
		};
	} M;
	class clB
	{
	public:
		clArray3Dd D[LB_NUMBER_CONNECTIONS], Ux[LB_NUMBER_CONNECTIONS], Uy[LB_NUMBER_CONNECTIONS], Uz[LB_NUMBER_CONNECTIONS];
		void ini()
		{
#ifdef LB_DEBUG
			for(int i = 1; i < c.lb.nL; i++)
			{
				D[i].zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
				Ux[i].zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
				Uy[i].zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
				Uz[i].zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			}
#endif //LB_DEBUG
		}
		void zeros()
		{
#ifdef LB_DEBUG
			for(int i = 1; i < c.lb.nL; i++)
			{
				D[i].zeros();
				Ux[i].zeros();
				Uy[i].zeros();
				Uz[i].zeros();
			}
#endif //LB_DEBUG
		}

		void update()
		{
		}

		void save2file()
		{
			char Buf[80];
			for(int i = 1; i < c.lb.nL; i++)
			{
				sprintf(Buf, "lbbd%d", i);
				D[i].savetxt(Buf);
				sprintf(Buf, "lbbux%d", i);
				Ux[i].savetxt(Buf);
				sprintf(Buf, "lbbuy%d", i);
				Uy[i].savetxt(Buf);
				sprintf(Buf, "lbbuz%d", i);
				Uz[i].savetxt(Buf);
			} 
		};
	} B;
	class clMeq
	{
	public:
		void ini()
		{
		}
		void update()
		{
		}
		double operator()(const int ind, clVector3Di ii)
		{
			return operator()(ind, vlb.Ro(ii), vlb.J(ii));
		}
		double operator()(const int ind, const double Ro, clVector3Dd J)
		{
			switch (ind)
			{
			case 0:
				return Ro;
			case 3:
				return J.x;
			case 5:
				return J.y;
			case 7:
				return J.z;
			case 1:
				return J.r2() / Ro;
			case 9:
				return (3. * J.x2() - J.r2()) / Ro;
			case 11:
				return (J.y2() - J.z2()) / Ro;
			case 13:
				return (J.x * J.y) / Ro;
			case 14:
				return (J.y * J.z) / Ro;
			case 15:
				return (J.x * J.z) / Ro;
			}
			return 0.;
		}
		void update(double *a, const double Ro, clVector3Dd J)
		{
			double Jx2 = J.x2(),  Jy2 = J.y2(),  Jz2 = J.z2();
			a[0] = Ro;
			a[1] = (Jx2 + Jy2 + Jz2) / Ro;
			a[2] = 0;
			a[3] = J.x;
			a[4] = 0;
			a[5] = J.y;
			a[6] = 0;
			a[7] = J.z;
			a[8] = 0;
			a[9] = (2. * Jx2 - Jy2 - Jz2) / Ro;
			a[10] = 0;
			a[11] = (Jy2 - Jz2) / Ro;
			a[12] = 0;
			a[13] = (J.x * J.y) / Ro;
			a[14] = (J.y * J.z) / Ro;
			a[15] = (J.x * J.z) / Ro;
			a[16] = 0;
			a[17] = 0;
			a[18] = 0;
		}

		void save2file()
		{
		};
	} Meq;
	class clS
	{
	public:
		void ini()
		{
		}
		void update()
		{
		}
		double operator()(const int ind, clVector3Di ii)
		{
			return operator()(ind, c.lb.lamda[vlb.M(ii)]);
		}
		double operator()(const int ind, const double lamda)
		{
			switch (ind)
			{
			case 1:
			case 2:
			case 4:
			case 6:
			case 8:
			case 10:
			case 12:
			case 16:
			case 17:
			case 18:
				return -1.;
			case 9:
			case 11:
			case 13:
			case 14:
			case 15:
				return lamda;
			}
			return 0.;
		}

		void save2file()
		{
		};
	} S;


	class clCoarse
	{
	public:
		int nX, nY, nZ, step;
		int X0, Y0, Z0;
		double lamda, lamdab;
		clVector3Di nn, c0;
		int dt;

		void ini()
		{
			dt = 0;
			step = c.lb.CoarseStep;
			double tauf = -1./c.lb.lamda[LB_FLUID_COARSE];
			double tauc = (tauf - 0.5) * step + 0.5;
			double taubf = -1./c.lb.lamdab[LB_FLUID_COARSE];
			double taubc = (taubf - 0.5) * step + 0.5;
//			lamda = -1./tauf;
			lamda = -1./tauc;
//			lamdab = c.lb.lamdab[LB_FLUID_COARSE];
			lamdab = -1./taubc;
			nX = c.lb.nXc;
			nY = c.lb.nYc;
			nZ = c.lb.nZc;

			X0 = c.lb.Xc0;
			Y0 = c.lb.Yc0;
			Z0 = c.lb.Zc0;

			nn = clVector3Di(nX, nY, nZ);
			c0 = clVector3Di(X0, Y0, Z0);

			M.ini();

			F.ini();
			Ro.ini();
			J.ini();
			PP.ini();
		};

		//int index_f(clVector3Di ii)
		//{
		//	return MOD(ii.x,step) + MOD(ii.y,step) * 2 + MOD(ii.z,step) * 4;
		//}
//		clVector3Di c2f(clVector3Di ii)
//		{
//			return (ii + c0) * step;
//		}
//		clVector3Di f2c(clVector3Di ii)
//		{
//			//clVector3Dd ret = (ii - 1);
//			//ret = ret / step - c0;
//			//ret = ret.floor();
//			//clVector3Di reti = (ii - 1 + step) / step - c0 - 1;
//			return (ii - 1 + step) / step - c0 - 1;
//		}
//		double Fc2f(int ind, clVector3Di iif)
//		{
//			clVector3Di iic = f2c(iif);
//			double wt_new = 0.25 + (double)dt * 0.5;
//			double wt_old = 1. - wt_new;
//			int ixf = MOD(iif.x,step), iyf = MOD(iif.y,step), izf = MOD(iif.z,step); 
//			double Ff = 0.;
//			for(int ix = 0; ix < 2; ix++)
//				for(int iy = 0; iy < 2; iy++)
//					for(int iz = 0; iz < 2; iz++)
//					{
//						double w = 1. / 64.;
//						w *= 1. + abs(ixf - ix) * 2.;
//						w *= 1. + abs(iyf - iy) * 2.;
//						w *= 1. + abs(izf - iz) * 2.;
//						clVector3Di iict = iic + clVector3Di(ix, iy, iz);
//
//						Ff += w * (wt_new * F.F[ind](iict) + wt_old * F.F[ind].o(iict));
////						Ff += w * F.F[ind](iict);
//					}
//
//			return Ff;
//		}
		clVector3Di c2f(clVector3Di ii)
		{
			return (ii + c0) * step - 1;
		}
		clVector3Di f2c(clVector3Di ii)
		{
			return ii / step - c0;
		}
		double Fc2f(int ind, clVector3Di iif)
		{
			clVector3Di iic = f2c(iif);
			double wt_new = 0. + (double)dt * 0.5;
			double wt_old = 1. - wt_new;
			int ixf = MOD(iif.x,step), iyf = MOD(iif.y,step), izf = MOD(iif.z,step); 
			double Ff = 0.;
			for(int ix = 0; ix < 2; ix++)
				for(int iy = 0; iy < 2; iy++)
					for(int iz = 0; iz < 2; iz++)
					{
						double w = 1. / 64.;
						w *= 1. + (1. - abs(ixf - ix)) * 2.;
						w *= 1. + (1. - abs(iyf - iy)) * 2.;
						w *= 1. + (1. - abs(izf - iz)) * 2.;
						clVector3Di iict = iic + clVector3Di(ix, iy, iz);

						Ff += w * (wt_new * F.F[ind](iict) + wt_old * F.F[ind].o(iict));
//						Ff += w * F.F[ind](iict);
					}

			return Ff;
		}
		double Ff2c(int ind, clVector3Di iic)
		{
			clVector3Di iif = c2f(iic);
			double Fc = 0.;
			for(int ix = 0; ix < 2; ix++)
				for(int iy = 0; iy < 2; iy++)
					for(int iz = 0; iz < 2; iz++)
						Fc += (vlb.F.F[ind](iif + clVector3Di(ix, iy, iz) - c.lb.C[ind])) / 8.;
//						Fc += (vlb.F.F[ind](iif + clVector3Di(ix, iy, iz)) + vlb.F.F[ind](iif + clVector3Di(ix, iy, iz))) / 16.;
			return Fc;
		}
		BOOL test_cinf(clVector3Di iic)
		{
			clVector3Di iif = c2f(iic);
			if(iif.x < 0 || iif.x+1 >= c.lb.nX)
				return FALSE;
			if(iif.y < 0 || iif.y+1 >= c.lb.nY)
				return FALSE;
			if(iif.z < 0 || iif.z+1 >= c.lb.nZ)
				return FALSE;
			return TRUE;
		}
		BOOL test_finc(clVector3Di iif)
		{
			clVector3Di iic = f2c(iif);
			if(iic.x < 0 || iic.x+1 >= nX)
				return FALSE;
			if(iic.y < 0 || iic.y+1 >= nY)
				return FALSE;
			if(iic.z < 0 || iic.z+1 >= nZ)
				return FALSE;
			return TRUE;
		}
		void update()
		{
			dt = MOD((int)c.t.Time,step);
			if(dt != 0)
				return;

			bcLB();
			propagationLB();
			collisionLBGK();
		};
		void propagationLB()
		{
			F.propagation();
		}
		void debug()
		{
			Ro.update();
			J.update();
			PP.update();
			save2file();
		}

		void collisionLBGK()
		{
			Ro.update();
			J.update();
			PP.update();
			F.change();
//			return;
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel default(shared) private(i)
#endif //LB_USE_OPENMP
			{
				clVector3Dd Fbody; 
				clTensorS3Dd Tunit;
				Tunit.unit();

#ifdef LB_USE_OPENMP
#pragma omp for
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
				{
					for(int j = 0; j < nY; j++)
					{
						for(int z = 0; z < nZ; z++)
						{
							if(M(i,j,z) != LB_SOLID)
							{
								Fbody = c.lb.F[LB_FLUID_COARSE];// * c.lb.FM(i,j,z);
								double Roc = Ro(i,j,z);
								clVector3Dd  Jc = J(i,j,z); 
								clTensorS3Dd PPc = PP(i,j,z);
								clTensorS3Dd RoUU = clTensorS3Dd(Jc) / Roc;

								clTensorS3Dd PPneq = PPc - (RoUU + Tunit * 1./3. * Roc);

								clTensorS3Dd PPneqtr = Tunit * 1./3. * PPneq.tr_sum();
								clTensorS3Dd PPneq_s = (PPneq - PPneqtr) * (1.+lamda) + PPneqtr * (1.+lamdab);

								clTensorS3Dd FFbody2 = clTensorS3Dd(Jc/Roc, Fbody);

								clVector3Dd JF = Jc + Fbody;
								clTensorS3Dd RoUUPF = RoUU + PPneq_s + FFbody2;

								for(int k = 0; k < c.lb.nL; k++)
								{
									const clVector3Dd C = c.lb.C[k];
									clTensorS3Dd CC = clTensorS3Dd(C) - Tunit * 1./3.;
									double Fnew = c.lb.F0[k] * (Roc + 3. * (JF * C) + RoUUPF * CC * 4.5);
									F(k)(i,j,z) = Fnew / c.lb.Ro[LB_FLUID_COARSE];
								}
							}
						}
					}
				}
			}
		}
		void bcLB()
		{
			bcFineLB();

			if(c.lb.bcPeriodicX == TRUE) bcPeriodicLBx();
			if(c.lb.bcPeriodicY == TRUE) bcPeriodicLBy();
			if(c.lb.bcPeriodicZ == TRUE) bcPeriodicLBz();

			if(c.lb.bcFreeX0 == TRUE) bcFreeFlowLBx0();
			if(c.lb.bcFreeX1 == TRUE) bcFreeFlowLBx1();
			if(c.lb.bcFreeY0 == TRUE) bcFreeFlowLBy0();
			if(c.lb.bcFreeY1 == TRUE) bcFreeFlowLBy1();
			if(c.lb.bcFreeZ0 == TRUE) bcFreeFlowLBz0();
			if(c.lb.bcFreeZ1 == TRUE) bcFreeFlowLBz1();

#ifdef LB_USE_OPENMP
#pragma omp parallel sections
			{
#pragma omp section
			{if(c.lb.bcWallX0 == TRUE) bcWallLBx0();}
#pragma omp section
			{if(c.lb.bcWallX1 == TRUE) bcWallLBx1();}
#pragma omp section
			{if(c.lb.bcWallY0 == TRUE) bcWallLBy0();}
#pragma omp section
			{if(c.lb.bcWallY1 == TRUE) bcWallLBy1();}
#pragma omp section
			{if(c.lb.bcWallZ0 == TRUE) bcWallLBz0();}
#pragma omp section
			{if(c.lb.bcWallZ1 == TRUE) bcWallLBz1();}

#pragma omp section
			{if(c.lb.bcSymmetryX0 == TRUE) bcSymmetryLBx0();}
#pragma omp section
			{if(c.lb.bcSymmetryX1 == TRUE) bcSymmetryLBx1();}
#pragma omp section
			{if(c.lb.bcSymmetryY0 == TRUE) bcSymmetryLBy0();}
#pragma omp section
			{if(c.lb.bcSymmetryY1 == TRUE) bcSymmetryLBy1();}
#pragma omp section
			{if(c.lb.bcSymmetryZ0 == TRUE) bcSymmetryLBz0();}
#pragma omp section
			{if(c.lb.bcSymmetryZ1 == TRUE) bcSymmetryLBz1();}
			}

#else //LB_USE_OPENMP
			if(c.lb.bcWallX0 == TRUE) bcWallLBx0();
			if(c.lb.bcWallX1 == TRUE) bcWallLBx1();
			if(c.lb.bcWallY0 == TRUE) bcWallLBy0();
			if(c.lb.bcWallY1 == TRUE) bcWallLBy1();
			if(c.lb.bcWallZ0 == TRUE) bcWallLBz0();
			if(c.lb.bcWallZ1 == TRUE) bcWallLBz1();

			if(c.lb.bcSymmetryX0 == TRUE) bcSymmetryLBx0();
			if(c.lb.bcSymmetryX1 == TRUE) bcSymmetryLBx1();
			if(c.lb.bcSymmetryY0 == TRUE) bcSymmetryLBy0();
			if(c.lb.bcSymmetryY1 == TRUE) bcSymmetryLBy1();
			if(c.lb.bcSymmetryZ0 == TRUE) bcSymmetryLBz0();
			if(c.lb.bcSymmetryZ1 == TRUE) bcSymmetryLBz1();
#endif //LB_USE_OPENMP

		}
		void bcFineLB()
		{
			int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
			for(i = 0; i < M.nB; i++)
				for(int k = 0; k < c.lb.nL; k++)
					F(k).n(M.B(i)) = Ff2c(k, M.B(i));
		}
		////////////////////////////////solid wall BC/////////////////////////////////////////////////////////////////////////
		void bcWallLB(int i, int j, int z, int dir, clVector3Dd UU2Cs2)
		{
			int dirp = c.lb.rev[dir];
			F(dirp).n(i,j,z) = F(dir)(i,j,z) - c.lb.F0[dir] * (c.lb.C[dir] * UU2Cs2);
		}
		void bcWallLBx0()
		{
			clVector3Dd U0 = c.lb.UX0 * 2. / c.lb.Cs2;

			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == -1)
							bcWallLB(0, j, z, k, U0);
				}
		}
		void bcWallLBx1()
		{
			int nXm = nX-1;
			clVector3Dd U1 = c.lb.UX1 * 2. / c.lb.Cs2;

			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == 1)
							bcWallLB(nXm, j, z, k, U1);
				}
		}

		void bcWallLBy0()
		{
			clVector3Dd U0 = c.lb.UY0 * 2. / c.lb.Cs2;

			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == -1)
							bcWallLB(i, 0, z, k, U0);
				}
		}
		void bcWallLBy1()
		{
			int nYm = nY-1;
			clVector3Dd U1 = c.lb.UY1 * 2. / c.lb.Cs2;

			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == 1)
							bcWallLB(i, nYm, z, k, U1);
				}
		}

		void bcWallLBz0()
		{
			clVector3Dd U0 = c.lb.UZ0 * 2. / c.lb.Cs2;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == -1)
							bcWallLB(i, j, 0, k, U0);
				}
		}
		void bcWallLBz1()
		{
			int nZm = nZ-1;
			clVector3Dd U1 = c.lb.UZ1 * 2. / c.lb.Cs2;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
				{
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == 1)
							bcWallLB(i, j, nZm, k, U1);
				}
		}
		////////////////////////////////symmetry BC/////////////////////////////////////////////////////////////////////////
		void bcSymmetryLB(int i, int j, int z, int ip, int jp, int zp, int dir, int ind)
		{
			int dirp = c.lb.rev[dir];
			if(dir > 2 * NUMBER_DIMENSIONS)
			{
				int dirn = c.lb.per[dir];
				if(c.lb.C[dir].e(ind) == -c.lb.C[dirn].e(ind))
					dirn = c.lb.rev[dirn];
				F(dirp).n(i,j,z) = F(dirn)(ip,jp,zp);
			}
			else
				F(dirp).n(i,j,z) = F(dir)(i,j,z);
		}
		void bcSymmetryLBx0()
		{
			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == -1)
							bcSymmetryLB(0, j, z, 1, j+c.lb.Cy[k], z+c.lb.Cz[k], k, 0);
		}
		void bcSymmetryLBx1()
		{
			int nXm = nX-1;

			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == 1)
							bcSymmetryLB(nXm, j, z, nXm-1, j+c.lb.Cy[k], z+c.lb.Cz[k], k, 0);
		}

		void bcSymmetryLBy0()
		{
			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == -1)
							bcSymmetryLB(i, 0, z, i+c.lb.Cx[k], 1, z+c.lb.Cz[k], k, 1);
		}
		void bcSymmetryLBy1()
		{
			int nYm = nY-1;

			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == 1)
							bcSymmetryLB(i, nYm, z, i+c.lb.Cx[k], nYm-1, z+c.lb.Cz[k], k, 1);
		}

		void bcSymmetryLBz0()
		{
			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == -1)
							bcSymmetryLB(i, j, 0, i+c.lb.Cx[k], j+c.lb.Cy[k], 1, k, 2);
		}
		void bcSymmetryLBz1()
		{
			int nZm = nZ-1;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == 1)
							bcSymmetryLB(i, j, nZm, i+c.lb.Cx[k], j+c.lb.Cy[k], nZm-1, k, 2);
		}
		/////////////////////////////////////free flow BC//////////////////////////////////////////	void bcFreeFlowLBx0()
		void bcFreeFlowLB(int i, int j, int z, int dir, clVector3Dd UU2Cs2)
		{
			int dirp = c.lb.rev[dir];
			F(dirp).n(i,j,z) = F(dir)(i,j,z) - c.lb.F0[dir] * (c.lb.C[dir] * UU2Cs2);
		}
		void bcFreeFlowLBx0()
		{
			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
				{
					clVector3Dd U0 = J(0, j, z) / Ro(0, j, z) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == -1)
							bcFreeFlowLB(0, j, z, k, U0);
				}
		}
		void bcFreeFlowLBx1()
		{
			int nXm = nX-1;
			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
				{
					clVector3Dd U1 = J(nXm, j, z) / Ro(nXm, j, z) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == 1)
							bcFreeFlowLB(nXm, j, z, k, U1);
				}
		}

		void bcFreeFlowLBy0()
		{
			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
				{
					clVector3Dd U0 = J(i, 0, z) / Ro(i, 0, z) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == -1)
							bcFreeFlowLB(i, 0, z, k, U0);
				}
		}
		void bcFreeFlowLBy1()
		{
			int nYm = nY-1;
			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
				{
					clVector3Dd U1 = J(i, nYm, z) / Ro(i, nYm, z) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == 1)
							bcFreeFlowLB(i, nYm, z, k, U1);
				}
		}
/*
		void bcFreeFlowLBz0()
		{
			clVector3Dd U0 = c.lb.UZ0 * 2. / c.lb.Cs2;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
				{
					clVector3Dd U0 = J(i, j, 0) / Ro(i, j, 0) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == -1)
							bcFreeFlowLB(i, j, 0, k, U0);
				}
		}
		void bcFreeFlowLBz1()
		{
			int nZm = nZ-1;
			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
				{
					clVector3Dd U1 = J(i, j, nZm) / Ro(i, j, nZm) * 2. / c.lb.Cs2;
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == 1)
							bcFreeFlowLB(i, j, nZm, k, U1);
				}
		}
*/
		//void bcFreeFlowLBx0()
		//{
		//	for(int j = 0; j < nY; j++)
		//		for(int z = 0; z < nZ; z++)
		//			for(int k = 1; k < c.lb.nL; k++)
		//				if(c.lb.Cx[k] == 1)
		//					F.F[k].npx(0,j,z) = F.F[k](1,j,z);
		//}
		//void bcFreeFlowLBx1()
		//{
		//	int nXm = nX-1;
		//	for(int j = 0; j < nY; j++)
		//		for(int z = 0; z < nZ; z++)
		//			for(int k = 1; k < c.lb.nL; k++)
		//				if(c.lb.Cx[k] == -1)
		//					F.F[k].npx(nXm,j,z) = F.F[k](nXm-1,j,z);
		//}

		//void bcFreeFlowLBy0()
		//{
		//	for(int i = 0; i < nX; i++)
		//		for(int z = 0; z < nZ; z++)
		//			for(int k = 1; k < c.lb.nL; k++)
		//				if(c.lb.Cy[k] == 1)
		//					F.F[k].npy(i,0,z) = F.F[k](i,1,z);
		//}
		//void bcFreeFlowLBy1()
		//{
		//	int nYm = nY-1;
		//	for(int i = 0; i < nX; i++)
		//		for(int z = 0; z < nZ; z++)
		//			for(int k = 1; k < c.lb.nL; k++)
		//				if(c.lb.Cy[k] == -1)
		//					F.F[k].npy(i,nYm,z) = F.F[k](i,nYm-1,z);
		//}

		void bcFreeFlowLBz0()
		{
			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == 1)
							F.F[k].npz(i,j,0) = F.F[k](i,j,1);
		}

		void bcFreeFlowLBz1()
		{
			int nZm = nZ-1;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == -1)
							F.F[k].npz(i,j,nZm) = F.F[k](i,j,nZm-1);
		}
		////////////////////////////////////////////periodic BC//////////////////////////////////////////
		void bcPeriodicLBx()
		{
			int nXm = nX-1;

			for(int j = 0; j < nY; j++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
					{
						if(c.lb.Cx[k] == -1)
							F(k).npx(nXm,j,z) = F(k)(0,j,z);
						if(c.lb.Cx[k] == 1)
							F(k).npx(0,j,z) = F(k)(nXm,j,z);
					}
		}

		void bcPeriodicLBy()
		{
			int nYm = nY-1;

			for(int i = 0; i < nX; i++)
				for(int z = 0; z < nZ; z++)
					for(int k = 1; k < c.lb.nL; k++)
					{
						if(c.lb.Cy[k] == -1)
							F(k).npy(i,nYm,z) = F(k)(i,0,z);
						if(c.lb.Cy[k] == 1)
							F(k).npy(i,0,z) = F(k)(i,nYm,z);
					}
		}

		void bcPeriodicLBz()
		{
			int nZm = nZ-1;

			for(int i = 0; i < nX; i++)
				for(int j = 0; j < nY; j++)
					for(int k = 1; k < c.lb.nL; k++)
					{
						if(c.lb.Cz[k] == -1)
							F(k).npz(i,j,nZm) = F(k)(i,j,0);
						if(c.lb.Cz[k] == 1)
							F(k).npz(i,j,0) = F(k)(i,j,nZm);
					}
		}
//////////////////////////////////////////////////////////////////////////////////////////////
		class clFc
		{
		public:
			clF3Dc F[LB_NUMBER_CONNECTIONS];


			clF3Dc& operator()(int n)
			{
				return F[n];
			};
			void ini()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].ini("lbfc", i, vlb.C.nn);

			}

			void update()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].update();
			}

			void propagation()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].propagation();
			}

			void copy()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].copy();
			}

			void change()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].change();
			}

			void save2file()
			{
#ifdef CREATE_BIN_FILES
				savetxt();
				savebin();
#endif //CREATE_BIN_FILES
			};
			void savetxt()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].savetxt();
			};
			void savebin()
			{
				for(int i = 0; i < c.lb.nL; i++)
					F[i].savebin();
			};
		} F;
		class clRo: public clArray3Dd
		{
		public:
			int nX, nY, nZ;
			clArray1Di counter;
			clArray1Dd avr, mass;
			void ini()
			{
				nX = vlb.C.nX;
				nY = vlb.C.nY;
				nZ = vlb.C.nZ;
				avr.zeros(c.lb.nF);
				mass.zeros(c.lb.nF);
				counter.zeros(c.lb.nF);
				clArray3Dd::zeros(nX, nY, nZ);
				clArray3Dd::fill(1);
				update();
			}

			void update()
			{
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int m = 0; m < nZ; m++)
							if(vlb.C.M(i,j,m) != LB_SOLID)
								update(clVector3Di(i, j, m));
							else
								el(i,j,m) = Rof2c(clVector3Di(i, j, m));
			}
			double Rof2c(clVector3Di iic)
			{
				clVector3Di iif = vlb.C.c2f(iic);
				double Jc = 0.;
				for(int ix = 0; ix < 2; ix++)
					for(int iy = 0; iy < 2; iy++)
						for(int iz = 0; iz < 2; iz++)
							Jc += vlb.Ro(iif + clVector3Di(ix, iy, iz));
				return Jc / 8.;
			}

			void update(clVector3Di ii)
			{
				double dF = 0.;
#ifndef LB_NO_MASS_CORRECTION
//				dF = -vlb.C.Ro.avr(vlb.C.M(ii)) / (double)(c.t.CurrentSaveTimeStep) / (double)c.lb.nL;
#endif //LB_NO_MASS_CORRECTION

				double F = 0.;
				for(int k = 0; k < c.lb.nL; k++)
					F += (vlb.C.F(k)(ii) += dF);

				el(ii) = F * c.lb.Ro[LB_FLUID_COARSE];
				
			}

			double d(int i, int j, int k)
			{
				return el(i,j,k) - c.lb.Ro[LB_FLUID_COARSE];
			}

			void save2file()
			{
				clArray3Dd copy(nX, nY, nZ);
				for(int i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int m = 0; m < nZ; m++)
							copy(i,j,m) = d(i,j,m);
				copy.savetxt("lbroc");
			}

			void updateoutput()
			{
				mass.zeros();
				counter.zeros();
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int k = 0; k < nZ; k++)
						{
							char map = vlb.C.M(i,j,k);
							if(map != LB_SOLID)
							{
#pragma omp atomic
								mass(map) += d(i,j,k);
#pragma omp atomic
								counter(map)++;
							}
						}
						for(int i = 1; i < c.lb.nF; i++)
							if(counter(i) > 0)
							{
								avr(i) = mass(i) / (double)counter(i);
								mass(i) += c.lb.Ro[i] * (double)counter(i);
							}
			}
		} Ro;

		class clJ: public clArray3Dv3d
		{
		public:
			int nX, nY, nZ;
			clArray1Dv3d avr, tot;
			clArray1Di counter;
			virtual void ini()
			{
				nX = vlb.C.nX;
				nY = vlb.C.nY;
				nZ = vlb.C.nZ;
				avr.zeros(c.lb.nF);
				tot.zeros(c.lb.nF);
				counter.zeros(c.lb.nF);
				clArray3Dv3d::zeros(nX, nY, nZ);
				update();
			}

			void update()
			{
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int m = 0; m < nZ; m++)
							if(vlb.C.M(i,j,m) != LB_SOLID)
								update(clVector3Di(i, j, m));
							else
								el(i,j,m) = Jf2c(clVector3Di(i, j, m));
			}
			clVector3Dd Jf2c(clVector3Di iic)
			{
				clVector3Di iif = vlb.C.c2f(iic);
				clVector3Dd Jc;
				for(int ix = 0; ix < 2; ix++)
					for(int iy = 0; iy < 2; iy++)
						for(int iz = 0; iz < 2; iz++)
							Jc += vlb.J(iif + clVector3Di(ix, iy, iz));
				return Jc / 8.;
			}

			void update(clVector3Di ii)
			{
				clVector3Dd J;
				for(int k = 1; k < c.lb.nL; k++)
				{
					double F = vlb.C.F(k)(ii);
					J += c.lb.C[k] * F;
				}
				el(ii) = J * c.lb.Ro[LB_FLUID_COARSE];
			}

			void updateoutput()
			{
				avr.zeros();
				counter.zeros();
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int k = 0; k < nZ; k++)
						{
							char map = vlb.C.M(i,j,k);
							if(map != LB_SOLID)
							{
#pragma omp atomic
								avr(map).x += el(i,j,k).x;
#pragma omp atomic
								avr(map).y += el(i,j,k).y;
#pragma omp atomic
								avr(map).z += el(i,j,k).z;
#pragma omp atomic
								counter(map)++;
							}
						}
						for(int i = 1; i < c.lb.nF; i++)
						{
							if(counter(i) > 0)
							{
								tot(i) = avr(i);
								avr(i) = avr(i) / (double)counter(i);
							}
						}
			}

			void save2file()
			{
				clArray3Dv3d copy;
				copy.zeros(nX, nY, nZ);
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int m = 0; m < nZ; m++)
							copy(i,j,m) = el(i,j,m) + c.lb.F[LB_FLUID_COARSE] * 0.5;// * c.lb.FM(i,j,m);
				copy.savetxt("lbjc");
			};
		} J;

		class clPP: public clArray3Dts3d
		{
		public:
			int nX, nY, nZ;
			virtual void ini()
			{
				nX = vlb.C.nX;
				nY = vlb.C.nY;
				nZ = vlb.C.nZ;
				clArray3Dts3d::zeros(nX, nY, nZ);
			}

			void update()
			{
				int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif //LB_USE_OPENMP
				for(i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int m = 0; m < nZ; m++)
							if(vlb.C.M(i,j,m) != LB_SOLID)
								update(clVector3Di(i, j, m));
			}

			void update(clVector3Di ii)
			{
				clTensorS3Dd PP;
				for(int k = 1; k < c.lb.nL; k++)
				{
					clTensorS3Dd C = clTensorS3Dd(c.lb.C[k]);
					double F = vlb.C.F(k)(ii);
					PP += C * F;
				}
				el(ii) = PP * c.lb.Ro[vlb.C.M(ii)];
			}

			void save2file()
			{
#ifdef LB_DEBUG
				savetxt("lbppc");
#endif //LB_DEBUG
			};
		} PP;

		class clM: public clArray3Dc
		{
		public:
			clArray1Dv3i B;
			int nB;
			int nX, nY, nZ;
			void ini()
			{
				nX = vlb.C.nX;
				nY = vlb.C.nY;
				nZ = vlb.C.nZ;
				nB = 0;
				clArray3Dc::zeros(nX, nY, nZ);
				fill(LB_FLUID_COARSE);
				for(int i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int k = 0; k < nZ; k++)
							if(vlb.C.test_cinf(clVector3Di(i,j,k)) == TRUE)
								el(i,j,k) = 0;

				for(int i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int k = 0; k < nZ; k++)
							if(el(i,j,k) == 0 && testfineboudnary(clVector3Di(i,j,k)) == TRUE)
								nB++;

				B.zeros(nB);
				int ind = 0;
				for(int i = 0; i < nX; i++)
					for(int j = 0; j < nY; j++)
						for(int k = 0; k < nZ; k++)
							if(el(i,j,k) == 0 && testfineboudnary(clVector3Di(i,j,k)) == TRUE)
								B(ind++) = clVector3Di(i,j,k);

				for(int i = 0; i < nB; i++)
					for(int k = 0; k < c.lb.nL; k++)
						el(B(i)) = LB_FLUID_COARSE;
				B.savetxt("lbmcb");
			}

			BOOL testfineboudnary(clVector3Di ii)
			{
				int xm = ii.x - 1, xp = ii.x + 1;
				int ym = ii.y - 1, yp = ii.y + 1;
				int zm = ii.z - 1, zp = ii.z + 1;
				if(xm < 0) xm = 0;
				if(xp >= nX) xp = nX - 1;
				if(ym < 0) ym = 0;
				if(yp >= nY) yp = nY - 1;
				if(zm < 0) zm = 0;
				if(zp >= nZ) zp = nZ - 1;

				for(int ix = xm; ix <= xp; ix++)
					for(int iy = ym; iy <= yp; iy++)
						for(int iz = zm; iz <= zp; iz++)
							if(el(ix, iy, iz) != 0)
								return TRUE;
				return FALSE;
			}
			void update()
			{
			};

			void save2file()
			{
				savetxt("lbmc");
#ifdef LB_DEBUG
#endif //LB_DEBUG
			};
		} M;

		void save2file()
		{
#ifdef LB_SOLVER
			F.save2file();
			Ro.save2file();
			J.save2file();
			PP.save2file();
			M.save2file();
#endif //LB_SOLVER
		}

	} C;

	void save2file()
	{
#ifdef LB_SOLVER
		F.save2file();
		Ro.save2file();
		J.save2file();
		PP.save2file();
		M.save2file();
		G.save2file();
		N.save2file();
#ifdef LB_COARSE
		C.save2file();
#endif //LB_COARSE
#endif //LB_SOLVER
	}
};
/////////////////////////////////////////////////////////////////////////

class clVariablesTr
{
	int updatetime;

public:
	void ini()
	{
		C.ini();
		V.ini();
		U.ini();
		updatetime = BEADS_UPDATE_TIME;


#ifdef TR_USE_DIFFUSION
		W.ini();
#endif

	};

	void update1()
	{		
#ifdef TR_BEADS
		if(MOD((int)c.t.Time,updatetime) == 0)
		{
			vtr.A.update_tables();
		}
#endif

		U.update();
#ifdef TR_BEADS
		A.Test_Inside_Beads();
#endif
		C.update();

	}

	void update2()
	{
		V.update();
#ifdef TR_USE_DIFFUSION
		W.update();
#endif

	};

	void updateoutput()
	{
		vtr.C.Calc_Concentration();
	}

	void save2file()
	{
		C.save2file();
		//V.save2file();
		A.save2file();

#ifdef TR_USE_DIFFUSION
		//W.save2file();
#endif

	}

	class clC: public clArray1Dv3d
	{
	public:
		clArray1Dv3i p;
		clArray1Dv3d d;
		clArray3Dd Conc;
		clVector3Dd OutputVal;

		void ini()
		{
			srand(1);
			c.tr.nN = TRC_NUM_TRACERS;
			clArray1Dv3d::zeros(c.tr.nN);
			p.zeros(c.tr.nN);
			c.tr.remain = c.tr.nN;
			d.zeros(c.tr.nN);
			vtr.A.ini();
			iniVolume_half();
			OutputVal = 0.;
			Conc.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			Calc_Concentration();
			shift();
			save2file();
		}

		void iniVolume()
		{
			double maxim = (double)RAND_MAX;
			double nXm = (double)c.lb.nX, nYm = (double)c.lb.nY, nZm = (double)c.lb.nZ;
			clVector3Dd pos;
			int flag_inside;

			for(int i = 0; i < c.tr.nN; i++)
			{
				flag_inside = 1;
				while(flag_inside == 1)
				{
					flag_inside = 0;
					pos.x = nXm*(double)rand()/maxim - 0.5;
					pos.y = nYm*(double)rand()/maxim - 0.5;
					pos.z = nZm*(double)rand()/maxim - 0.5;
				
#ifdef TR_RIDGES
					if(pos.z >= vtr.A.Max_vec.z)
					{
						if(test_inside_ridges(pos) == TRUE)
							flag_inside = 1;
					}
#endif

#ifdef TR_BEADS			
					if(pos.z <= vtr.A.Min_vec.z)
					{
						if(test_inside_beads(pos) == TRUE)
							flag_inside = 1;
					}
#endif

#ifndef TR_HIDDEN_DISKS
					if(pos.z <= PILLAR_H2)
					{
						if(test_inside_disks(pos) == TRUE)
							flag_inside = 1;
					}
#endif

				}
				el(i) = pos;
			}
		}

		void Calc_Concentration()
		{
			Conc.zeros();
			int i;
			clArray3Di CountA;
			CountA.zeros(c.nn);

#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{
				clVector3Dd Pos = vtr.C.s(i) + 0.5;
				clVector3Dd Posi = Pos.floor();
				double inc;
				if(vlb.M(Posi) == LB_SOLID)
					continue;
				if(vtr.A(i) == 1)
					inc = 1.;
				else
					inc = -1.;	
#pragma omp atomic
				CountA(Posi)++;
#pragma omp atomic
				Conc(Posi) += inc;
			}

			double C_tot = 0.;
			int count = 0;

#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
			{
				for(int j = 0; j < c.lb.nY; j++)
				{
					for(int k = 0; k < c.lb.nZ; k++)
					{
						if(vlb.M(i,j,k) != LB_SOLID)
						{
							if(CountA(i,j,k) != 0)
								Conc(i,j,k) /= (double)CountA(i,j,k);
#pragma omp atomic
							C_tot += Conc(i,j,k);
#pragma omp atomic
							count++;
						}
					}
				}
			}

			double Nnode_f = (double)count;
			double C_avr = C_tot/(double)Nnode_f;

			double C_diff = 0.;

#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
			{
				for(int j = 0; j < c.lb.nY; j++)
				{
					for(int k = 0; k < c.lb.nZ; k++)
					{
						if(vlb.M(i,j,k) != LB_SOLID)
						{
							double C_diff1 = (Conc(i,j,k) - C_avr) * (Conc(i,j,k) - C_avr);
#pragma omp atomic
							C_diff += C_diff1;
						}

					}
				}
			}

			double sigma = sqrt(C_diff/(double)Nnode_f);

			OutputVal = clVector3Dd(sigma,C_avr,C_diff);

		}

		void iniVolume_half()
		{
			double maxim = (double)RAND_MAX+1.;
			double nXm = (double)c.lb.nX, nYm = (double)c.lb.nY, nZm = (double)c.lb.nZ;
			clVector3Dd pos;
			int flag_inside;

			for(int i = 0; i < c.tr.nN; i++)
			{
				flag_inside = 1;
				while(flag_inside == 1)
				{
					flag_inside = 0;
					pos.x = nXm*(double)rand()/maxim - 0.5;
					pos.y = nYm*(double)rand()/maxim - 0.5;
					pos.z = nZm*(double)rand()/maxim - 0.5;
				
#ifdef TR_RIDGES
					if(pos.z >= vtr.A.Max_vec.z)
					{
						if(test_inside_ridges(pos) == TRUE)
							flag_inside = 1;
					}
#endif

#ifdef TR_BEADS			
					if(pos.z <= vtr.A.Min_vec.z)
					{
						if(test_inside_beads(pos) == TRUE)
							flag_inside = 1;
					}
#endif

#ifndef TR_HIDDEN_DISKS
					if(pos.z <= PILLAR_H2)
					{
						if(test_inside_disks(pos) == TRUE)
							flag_inside = 1;
					}
#endif

				}
				el(i) = pos;
				if(pos.y > (double)c.lb.nY/2. - 0.5)
					vtr.A(i) = 1;
				else
					vtr.A(i) = 0;

			}
		}

		BOOL test_inside_beads(clVector3Dd C2_C)
		{
			clVector3Dd vLd;
			int bl;
			clNode<int> *tempt;

			for(int j = 0; j < N_BEAD; j++)
			{
				clVector3Dd Dist = C2_C - vtr.A.Bead_center(j);
				if(Dist.r() <= R0_NUMBER*0.98)
					return TRUE;
				if(Dist.r() <= R0_NUMBER*1.01)
				{
					clVector3Dd C2i = C2_C + 0.5;
					clVector3Di ind = C2i.floor();
					
					tempt = vtr.A.B_Table(ind.x,ind.y,ind.z).firstlocal();
					
					while(tempt != NULL)
					{
						bl = tempt->getData();
						if(c.ls.O(c.ls.BL(bl,0)) < N_BEAD)
						{
							if(TestParInBead(C2_C,vls.C(c.ls.BL(bl,0)),vls.C(c.ls.BL(bl,1)),vls.C(c.ls.BL(bl,2)),vLd) == TRUE)
							{
								return TRUE;
							}
							else
							{
								tempt = tempt->getLink();
							}
						}
						else
						{
							tempt = tempt->getLink();
						}
					}
				}
			}
			return FALSE;
		}

		BOOL TestParInBead(clVector3Dd vL0, clVector3Dd vP0, clVector3Dd vP1, 
			clVector3Dd vP2, clVector3Dd &vLd)
		{
			clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
			clVector3Dd vNm = vP10 ^ vP20;
			clVector3Dd vN = vNm.n();
			clVector3Dd vPL = vP0 - vL0;
			double dist = vPL*vN;
			vLd = vN*dist;
			if(dist >= 0.)
			{
				clVector3Dd vC = vL0 + vLd;
				if(vtr.A.testInside(vC,vP0,vP1,vP2) == TRUE)
					return TRUE;
			}

			return FALSE;

		}

		BOOL test_inside_disks(clVector3Dd C2_C)
		{
			//clVector3Dd Disk_center = clVector3Dd(CX0_DISK,CY0_DISK,0.);
			if(C2_C.z > PILLAR_H2)
				return FALSE;

			for(int j = 0; j < N_DISK; j++)
			{
				clVector3Dd Dist = C2_C - vtr.A.Disk_center(j);
				Dist.z = 0.;
				if(Dist.r() <= PILLAR_RADIUS)
				{
					return TRUE;
				}
			}
			return FALSE;
		}

		BOOL test_inside_ridges(clVector3Dd vt)
		{
			clVector3Dd v0, v1, v2, v3;
			for(int i = 0; i < vlb.M.nS; i++)
			{
				v0 = vls.C(vlb.M.S(i,0));
				v1 = vls.C(vlb.M.S(i,1));
				v2 = vls.C(vlb.M.S(i,2));
				v3 = vls.C(vlb.M.S(i,3));
				clVector3Dd v10 = v1 - v0,  v20 = v2 - v0,  v30 = v3 - v0;
				clVector3Dd vt0 = vt - v0,  vt1 = vt - v1,  vt2 = vt - v2,  vt3 = vt - v3;
				double Stot = volume(v10, v20, v30);
				double S012 = volume(vt0, vt1, vt2), S023 = volume(vt0, vt2, vt3), S013 = volume(vt0, vt1, vt3), S123 = volume(vt1, vt2, vt3);
				double Stest = S012 + S023 + S013 + S123;
				if(fabs(Stot - Stest) < c.eps)
					return TRUE;
			}
			vt += 0.5;
			clVector3Di ind = vt.floor();
			int numt = 0;
			if(ind.x == 0)
				numt +=  1;
			if(ind.x == c.lb.nX-1)
				numt += 2;
			if(ind.y == 0)
				numt +=  4;
			if(ind.y == c.lb.nY-1)
				numt += 8;
			if(ind.z == c.lb.nZ-1)
				numt += 16;

			switch(numt)
			{
			case 0:
				{
					if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
						vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
						vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
						vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 1:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 2:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 4:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 8:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 16:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 5:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 9:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 17:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 21:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID) 
							return TRUE;
				}
				break;
			case 25:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(c.lb.nX-1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 6:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID &&
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 10:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID &&
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z+1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 18:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID &&
							vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 22:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 26:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(0,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 20:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y+1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			case 24:
				{
						if(vlb.M(ind.x,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x+1,ind.y,ind.z) == LB_SOLID && vlb.M(ind.x-1,ind.y,ind.z) == LB_SOLID && 
							vlb.M(ind.x,ind.y-1,ind.z) == LB_SOLID && vlb.M(ind.x,ind.y,ind.z-1) == LB_SOLID)
							return TRUE;
				}
				break;
			}

			return FALSE;
		}

		double volume(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
		{
			return fabs(v0 * (v1 ^ v2) / 6.);
		}

		void dzeros()
		{
			d.zeros();
		}
						
		void shift()
		{
			int i;
#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{
				shiftx(i);
				shifty(i);
				shiftz(i);
			}
		}

		void shift_yz()
		{
			int i;
#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{
				shifty(i);
				shiftz(i);
			}
		}

		void shiftx(int i)
		{
			double xi = el(i).x + 0.5;
			int cxi = (int)xi;
			int pxi = (cxi - cxi % c.nX) / c.nX;
			if(xi < 0.) pxi--;
			p(i).x = pxi;
		}

		void shifty(int i)
		{
			double yi = el(i).y + 0.5;
			int cyi = (int)yi;
			int pyi = (cyi - cyi % c.nY) / c.nY;
			if(yi < 0) pyi--;
			p(i).y = pyi;
		}

		void shiftz(int i)
		{
			double zi = el(i).z + 0.5;
			int czi = (int)zi;
			int pzi = (czi - czi % c.nZ) / c.nZ;
			if(zi < 0) pzi--;
			p(i).z = pzi;
		}

		clVector3Dd s(int n)
		{
			return el(n) - (clVector3Dd)(p(n) | c.nn);
			//return clVector3Dd(sx(n), sy(n), sz(n));
		}

		clVector3Dd s(int n, int np)
		{
			return el(n) - (clVector3Dd)(p(np) | c.nn);
			//return clVector3Dd(sx(n,np), sy(n,np), sz(n,np));
		}

		double sx(int i)
		{
			return sx(i, i);
		}

		double sy(int i)
		{
			return sy(i, i);
		}

		double sz(int i)
		{
			return sz(i, i);
		}

		double sx(int ix, int ip)
		{
			return el(ix).x - (double)(p(ip).x * c.nX);
		}

		double sy(int iy, int ip)
		{
			return el(iy).y - (double)(p(ip).y * c.nY);
		}

		double sz(int iz, int ip)
		{
			return el(iz).z - (double)(p(ip).z * c.nZ);
		}

		void update()
		{
			shift();
			//shift_yz();
		}

		void save2file()
		{
#ifdef TR_DEBUG
			p.savetxt("trcp");
#endif //TR_DEBUG
			save("trc");
			Conc.savetxt("trConc");
		};
	} C;


	class clD: public clArray2Dv3d
	{
	public:
		int ind;
		//clArray2Dv3d Part;
		int Part_ind;
		void ini()
		{
			ind = 0;
			zeros(1000,8);
			//Part.zeros(10,2);
			//Part_ind = 0;
		}

		void clear()
		{
			clArray2Dv3d::zeros();
			ind = 0;
		}

		void save2file()
		{
			savetxt("trdebug");
		};
	} D;

	class clV: public clArray1Dv3d
	{
	public:
		void ini()
		{
			
			zeros(c.tr.nN);
			if(load("load" SLASH "trv"))
			{
			}
		}

		void update()
		{
			zeros();
			int i;
#ifdef TR_USE_OPENMP
#pragma omp parallel default(shared) private(i)
#endif //TR_USE_OPENMP
			{
#ifdef TR_USE_OPENMP
#pragma omp for
#endif //TR_USE_OPENMP
				for(i = 0; i < c.tr.nN; i++)
				{
					update(i);
				}
			}
		}

		void update(int n)
		{
			clVector3Dd cs = vtr.C.s(n) + 1.;
			clVector3Di ii = clVector3Di(cs.floor());
			double dX1 = cs.x - (double)ii.x, dX0 = 1. - dX1;
			double dY1 = cs.y - (double)ii.y, dY0 = 1. - dY1;
			double dZ1 = cs.z - (double)ii.z, dZ0 = 1. - dZ1;


			clVector3Dd U00z = vtr.U(ii.x, ii.y, ii.z) * dZ0 + vtr.U(ii.x, ii.y, ii.z+1) * dZ1;
			clVector3Dd U10z = vtr.U(ii.x+1, ii.y, ii.z) * dZ0 + vtr.U(ii.x+1, ii.y, ii.z+1) * dZ1;
			clVector3Dd U01z = vtr.U(ii.x, ii.y+1, ii.z) * dZ0 + vtr.U(ii.x, ii.y+1, ii.z+1) * dZ1;
			clVector3Dd U11z = vtr.U(ii.x+1, ii.y+1, ii.z) * dZ0 + vtr.U(ii.x+1, ii.y+1, ii.z+1) * dZ1;

			clVector3Dd U0yz = U00z * dY0 + U01z * dY1;
			clVector3Dd U1yz = U10z * dY0 + U11z * dY1;

			clVector3Dd Uxyz = U0yz * dX0 + U1yz * dX1;

			el(n) = Uxyz;
		}

		void save2file()
		{
			savetxt("trv");
		};
	} V;

	class clW: public clArray1Dv3d
	{

#pragma warning(disable: 4018 4244 4615 4305 4390)
	public:
		unsigned long kn[128];
		double wn[128], fn[128];
		int thread_num;
		unsigned long int *seed;
		clArray1Di Start, Stop;

		void r4_nor_setup (void)
		{
			double dn = 3.442619855899;
			int i;
			const double m1 = 2147483648.0;
			double q;
			double tn = 3.442619855899;
			const double vn = 9.91256303526217E-03;

			q = vn / exp ( - 0.5 * dn * dn );

			kn[0] = ( int ) ( ( dn / q ) * m1 );
			kn[1] = 0;

			wn[0] = ( double ) ( q / m1 );
			wn[127] = ( double ) ( dn / m1 );

			fn[0] = 1.0;
			fn[127] = ( double ) ( exp ( - 0.5 * dn * dn ) );

			for ( i = 126; 1 <= i; i-- )
			{
				dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
				kn[i+1] = ( int ) ( ( dn / tn ) * m1 );
				tn = dn;
				fn[i] = ( double ) ( exp ( - 0.5 * dn * dn ) );
				wn[i] = ( double ) ( dn / m1 );
			}
			return;
		}
		
		unsigned long int shr3(unsigned long int *jsr)
		{
			unsigned long int value;
			value = *jsr;
			*jsr = (*jsr ^ (*jsr <<   13));
			*jsr = (*jsr ^ (*jsr >>   17));
			*jsr = (*jsr ^ (*jsr <<    5));
			value = value + *jsr;
			return value;
		}

		double r4_nor (unsigned long int *jsr)
		{
			int hz;
			int iz;
			const double r = 3.442620;
			double value;
			double x;
			double y;

			hz = shr3 ( jsr );
			iz = ( hz & 127 );

			if ( abs ( hz ) < kn[iz] )
				value = ( double ) ( hz ) * wn[iz];
			else
			{
				for ( ; ; )
				{
					if ( iz == 0 )
					{
						for ( ; ; )
						{
							x = - 0.2904764 * log ( r4_uni ( jsr ) );
							y = - log ( r4_uni ( jsr ) );
							if ( x * x <= y + y );
								break;
						}
						if ( hz <= 0 )
							value = - r - x;
						else
							value = + r + x;
						
						break;
					}
					x = ( double ) ( hz ) * wn[iz];
					if ( fn[iz] + r4_uni ( jsr ) * ( fn[iz-1] - fn[iz] ) < exp ( - 0.5 * x * x ) )
					{
						value = x;
						break;
					}
					hz = shr3 ( jsr );
					iz = ( hz & 127 );
					if ( abs ( hz ) < kn[iz] )
					{
						value = ( double ) ( hz ) * wn[iz];
						break;
					}
				}
			}
			return value;
		}

		float r4_uni ( unsigned long int *jsr )
		{
		  unsigned long int jsr_input;
		  float value;
		  jsr_input = *jsr;
		  *jsr = ( *jsr ^ ( *jsr <<   13 ) );
		  *jsr = ( *jsr ^ ( *jsr >>   17 ) );
		  *jsr = ( *jsr ^ ( *jsr <<    5 ) );
		  value = fmod ( 0.5 + ( float ) ( jsr_input + *jsr ) / 65536.0 / 65536.0, 1.0 );
		  return value;
		}
		
		void ini()
		{
			zeros(c.tr.nN);
			
			#pragma omp parallel
			{
				#pragma omp master 
				{
					thread_num = omp_get_num_threads();
				}
			}

			r4_nor_setup();
			
#ifdef _DEBUG
			seed = new unsigned long int[1];
			
			seed[0] = (unsigned long int)rand();

			Start.zeros(1);
			Stop.zeros(1);

			Start(0) = 0;
			Stop(0) = c.tr.nN;
			thread_num = 1;
#else
			seed = new unsigned long int[thread_num];
			
			for (int thread = 0; thread < thread_num; thread++)
			{
				seed[thread] = (unsigned long int)rand();
			}

			int count = (int)floor(c.tr.nN/(double)thread_num);

			Start.zeros(thread_num);
			Stop.zeros(thread_num);

			Start(0) = 0;
			Stop(0) = count;

			for(int i = 1; i < thread_num-1; i++)
			{
				Start(i) = i*count;
				Stop(i) = (i+1)*count;
			}

			Start(thread_num-1) = Stop(thread_num-2);
			Stop(thread_num-1) = c.tr.nN;
#endif
		}


		void update()
		{	
			int i;
			double truncate = c.tr.TRUNC;
#pragma omp parallel private (i)
			{
#pragma omp for
				for(i = 0; i < thread_num; i++)
				{
					int thread = omp_get_thread_num();
					unsigned long int jsr = seed[thread];
					for(int j = Start(i); j < Stop(i); j++)
					{
						double X,Y,Z;
						do 
						{
							X = r4_nor(&jsr);
						} 
						while ( X > truncate || X < (0.-truncate));
						
						
						el(j).x = X;

						do 
						{
							Y = r4_nor(&jsr);
						} 
						while ( Y > truncate || Y < (0.-truncate));

						el(j).y = Y;

						do 
						{
							Z = r4_nor(&jsr);
						} 
						while ( Z > truncate || Z < (0.-truncate));
						el(j).z = Z;	
					}
					seed[thread] = jsr;
				}
			}
		}
#pragma warning(restore: 4018 4244 4615 4305 4390)

		void save2file()
		{
			savetxt("trw");
		};
	} W;

	class clA: public clArray1Di
	{
	public:
		clArray1Di BLdep;
		clArray3DLli T_Table, B_Table;
		clArray1Dv3d Bead_center, Disk_center;
		clArray2Di Bead_BL;
		int Bead_nBL;
		clArray1Di Bead_Num;
		clArray2Di Disk_BL;
		int Disk_nBL;
		clArray1Di Disk_Num;
		clArray1Di Surf_Bool; //boolean representing type of surface for each BL
		clArray1Di Bead_BL_list, Disk_BL_list;

		clVector3Dd Min_vec, Max_vec;

		void ini()
		{
			zeros(c.tr.nN);

			FindMinMax();
			ini_tables();

			if(load("load" SLASH "tra") == TRUE)
			{
			}
		}


		void ini_bead_centers()
		{
			Bead_center.zeros(N_BEAD);
			clArray1Dv3d Pos;
			clArray1Di count;
			Pos.zeros(N_BEAD);
			count.zeros(N_BEAD);
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.O(i) < N_BEAD)
				{
					Pos(c.ls.O(i)) += vls.C(i);
					count(c.ls.O(i))++;
				}
			}

			for(int i = 0; i < N_BEAD; i++)
				Bead_center(i) = Pos(i)/((double)count(i));

		}

		void ini_disk_centers()
		{
			Disk_center.zeros(N_DISK);
			clArray1Dv3d Pos;
			clArray1Di count;
			Pos.zeros(N_DISK);
			count.zeros(N_DISK);
			for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.O(i) >= N_BEAD && c.ls.O(i) < N_BEAD + N_DISK)
				{
					Pos(c.ls.O(i)-N_BEAD) += vls.C(i);
					count(c.ls.O(i)-N_BEAD)++;
				}
			}

			for(int i = 0; i < N_DISK; i++)
				Disk_center(i) = Pos(i)/((double)count(i));

		}



		void ini_tables()
		{
			int StartInd = 0;
#ifndef TR_BEADS
			StartInd = 0;
#else
			StartInd = N_BEAD+N_DISK;
#endif
			Surf_Bool.zeros(c.ls.nBL);
#ifdef TR_RIDGES
			T_Table.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);

			for(int i = 0; i < c.ls.nBL; i++)
			{
				int CC = c.ls.BL(i,0);
				int objn = c.ls.O(CC);
				if(objn >= StartInd)
				{
					Surf_Bool(i) = RIDGE_SURF;

					clVector3Dd Maxp = max3(vls.C(c.ls.BL(i,0)),vls.C(c.ls.BL(i,1)),vls.C(c.ls.BL(i,2)));
					clVector3Dd Minp = min3(vls.C(c.ls.BL(i,0)),vls.C(c.ls.BL(i,1)),vls.C(c.ls.BL(i,2)));
					if(Minp.z < Max_vec.z)
						Max_vec.z = Minp.z;

					
					Maxp += 0.5;
					Minp += 0.5;

					clVector3Di Maxi = Maxp.floor();
					clVector3Di Mini = Minp.floor();

					Maxi += 1;
					Mini -= 1;
					
					if(Mini.y >= c.lb.nY)
						continue;

					if(Maxi.y < 0)
						continue;

					if(Mini.y < 0)
						Mini.y = 0;

					if(Maxi.y >= c.lb.nY)
						Maxi.y = c.lb.nY-1;

					if(Mini.z < 0)
						Mini.z = 0;
					
					if(Maxi.z >= c.lb.nZ)
						Maxi.z = c.lb.nZ-1;

					if(Mini.z >= c.lb.nZ)
						Mini.z = c.lb.nZ-1;

					for(int ii = Mini.x; ii <= Maxi.x; ii++)
						for(int jj = Mini.y; jj <= Maxi.y; jj++)
							for(int kk = Mini.z; kk <= Maxi.z; kk++)
								T_Table(MODFAST(ii,c.lb.nX),MODFAST(jj,c.lb.nY),kk).appendItem(i);

				}

			}
			Max_vec.z -= 0.2;
#endif

#ifdef TR_BEADS
			B_Table.zeros(c.lb.nX,c.lb.nY,c.lb.nZ);
			int ind = 0;


#ifndef TR_HIDDEN_DISKS		
			ini_disk_centers();
			clArray2Di Disk_BLt;
			Disk_BLt.zeros(c.ls.nBL,3);
			clArray1Di Disk_Numt, Disk_BL_listt;
			Disk_Numt.zeros(c.ls.nBL);
			Disk_BL_listt.zeros(c.ls.nBL);


			for(int i = 0; i < c.ls.nBL; i++)
			{
				int CC = c.ls.BL(i,0);
				int objn = c.ls.O(CC);
				if(objn >= N_BEAD && objn < N_BEAD+N_DISK)
				{
					int dnum = objn - N_BEAD;
					Surf_Bool(i) = PILLAR_SURF;
					clVector3Dd vP0 = vls.C(c.ls.BL(i,0)), vP1 = vls.C(c.ls.BL(i,1)),  vP2 = vls.C(c.ls.BL(i,2));  
					clVector3Dd vP0C = vP0-Disk_center(dnum);
					clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
					clVector3Dd vNm = vP10 ^ vP20;
					double dist2 = vNm*vP0C;
					clVector3Dd Maxp = max3(vls.C(c.ls.BL(i,0)),vls.C(c.ls.BL(i,1)),vls.C(c.ls.BL(i,2)));
					if(Maxp.z > -0.5 && dist2 > 0.)
					{
						Disk_BLt(ind,0) = c.ls.BL(i,0);
						Disk_BLt(ind,1) = c.ls.BL(i,1);
						Disk_BLt(ind,2) = c.ls.BL(i,2);
						Disk_Numt(ind) = c.ls.O(CC);
						Disk_BL_listt(ind) = i;
						ind++;
					}
				}
			}

			Disk_BL.zeros(ind,3);
			Disk_nBL = ind;
			Disk_Num.zeros(ind);
			Disk_BL_list.zeros(ind);

			for(int i = 0; i < ind; i++)
			{
				Disk_BL(i,0) = Disk_BLt(i,0);
				Disk_BL(i,1) = Disk_BLt(i,1);
				Disk_BL(i,2) = Disk_BLt(i,2);
				Disk_Num(i) = Disk_Numt(i);
				Disk_BL_list(i) = Disk_BL_listt(i);
			}

			//Disk_BL.savetxt("Disk_BL");
			//Disk_BL_list.savetxt("DiskBL_list");

			for(int i = 0; i < Disk_nBL; i++)
			{
				int j = Disk_BL_list(i);

				clVector3Dd Maxp = max3(vls.C(Disk_BL(i,0)),vls.C(Disk_BL(i,1)),vls.C(Disk_BL(i,2)));
				clVector3Dd Minp = min3(vls.C(Disk_BL(i,0)),vls.C(Disk_BL(i,1)),vls.C(Disk_BL(i,2)));
				Maxp += 0.5;
				Minp += 0.5;

				clVector3Di Maxi = Maxp.floor();
				clVector3Di Mini = Minp.floor();
				
				Maxi += 1;
				Mini -= 1;

				if(Mini.z >= c.lb.nZ)
					Mini.z = c.lb.nZ-1;

				if(Mini.z < 0)
					Mini.z = 0;



				for(int ii = Mini.x; ii <= Maxi.x; ii++)
					for(int jj = Mini.y; jj <= Maxi.y; jj++)
						for(int kk = Mini.z; kk <= Maxi.z; kk++)
							B_Table(MODFAST(ii,c.lb.nX),MODFAST(jj,c.lb.nY),MODFAST(kk,c.lb.nZ)).appendItem(j);
			}

#endif

			ini_bead_centers();
			double rad = R0_NUMBER*0.95;
			
			clArray2Di Bead_BLt;
			Bead_BLt.zeros(c.ls.nBL,3);
			ind = 0;
			clArray1Di Bead_Numt, Bead_BL_listt;
			Bead_Numt.zeros(c.ls.nBL);
			Bead_BL_listt.zeros(c.ls.nBL);

			for(int i = 0; i < c.ls.nBL; i++)
			{
				int CC = c.ls.BL(i,0);
				int objn = c.ls.O(CC);
				if(objn < N_BEAD)
				{
					Surf_Bool(i) = BEAD_SURF;
					clVector3Dd distv = Bead_center(objn) - vls.C(CC);
					clVector3Dd vP0 = vls.C(c.ls.BL(i,0)), vP1 = vls.C(c.ls.BL(i,1)),  vP2 = vls.C(c.ls.BL(i,2));  
					clVector3Dd vP0C = vP0-Bead_center(objn);
					clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
					clVector3Dd vNm = vP10 ^ vP20;
					double dist2 = vNm*vP0C;
					double distd = distv.r();
					if(distd >= rad && dist2 > 0.)
					{
						Bead_BLt(ind,0) = c.ls.BL(i,0);
						Bead_BLt(ind,1) = c.ls.BL(i,1);
						Bead_BLt(ind,2) = c.ls.BL(i,2);
						Bead_Numt(ind) = c.ls.O(CC);
						Bead_BL_listt(ind) = i;
						ind++;
					}
				}
			}


			Bead_BL.zeros(ind,3);
			Bead_nBL = ind;
			Bead_Num.zeros(ind);
			Bead_BL_list.zeros(ind);

			for(int i = 0; i < ind; i++)
			{
				Bead_BL(i,0) = Bead_BLt(i,0);
				Bead_BL(i,1) = Bead_BLt(i,1);
				Bead_BL(i,2) = Bead_BLt(i,2);
				Bead_Num(i) = Bead_Numt(i);
				Bead_BL_list(i) = Bead_BL_listt(i);
			}

			for(int i = 0; i < Bead_nBL; i++)
			{
				int j = Bead_BL_list(i);

				clVector3Dd Maxp = max3(vls.C(Bead_BL(i,0)),vls.C(Bead_BL(i,1)),vls.C(Bead_BL(i,2)));
				clVector3Dd Minp = min3(vls.C(Bead_BL(i,0)),vls.C(Bead_BL(i,1)),vls.C(Bead_BL(i,2)));
				if(Maxp.z > Min_vec.z)
						Min_vec.z = Maxp.z;

				Maxp += 0.5;
				Minp += 0.5;

				clVector3Di Maxi = Maxp.floor();
				clVector3Di Mini = Minp.floor();

				Maxi += 1;
				Mini -= 1;

				if(Mini.y >= c.lb.nY)
					continue;

				if(Maxi.y < 0)
					continue;

				if(Mini.y < 0)
					Mini.y = 0;

				if(Maxi.y >= c.lb.nY)
					Maxi.y = c.lb.nY-1;

				if(Mini.z >= c.lb.nZ)
					Mini.z = c.lb.nZ-1;

				if(Mini.z < 0)
					Mini.z = 0;

				for(int ii = Mini.x; ii <= Maxi.x; ii++)
					for(int jj = Mini.y; jj <= Maxi.y; jj++)
						for(int kk = Mini.z; kk <= Maxi.z; kk++)
							B_Table(MODFAST(ii,c.lb.nX),MODFAST(jj,c.lb.nY),kk).appendItem(j);
			}

			Min_vec.z += 0.2;


#endif

		}

		void update_tables()
		{
			B_Table.clear();	

#ifndef TR_HIDDEN_DISKS		
			for(int i = 0; i < Disk_nBL; i++)
			{
				int j = Disk_BL_list(i);

				clVector3Dd Maxp = max3(vls.C(Disk_BL(i,0)),vls.C(Disk_BL(i,1)),vls.C(Disk_BL(i,2)));
				clVector3Dd Minp = min3(vls.C(Disk_BL(i,0)),vls.C(Disk_BL(i,1)),vls.C(Disk_BL(i,2)));
				Maxp += 0.5;
				Minp += 0.5;

				clVector3Di Maxi = Maxp.floor();
				clVector3Di Mini = Minp.floor();
				
				Maxi += 1;
				Mini -= 1;

				if(Mini.z >= c.lb.nZ)
					Mini.z = c.lb.nZ-1;

				if(Mini.z < 0)
					Mini.z = 0;

				for(int ii = Mini.x; ii <= Maxi.x; ii++)
					for(int jj = Mini.y; jj <= Maxi.y; jj++)
						for(int kk = Mini.z; kk <= Maxi.z; kk++)
							B_Table(MODFAST(ii,c.lb.nX),MODFAST(jj,c.lb.nY),MODFAST(kk,c.lb.nZ)).appendItem(i);
			}
#endif

			for(int i = 0; i < Bead_nBL; i++)
			{
				int j = Bead_BL_list(i);
	
				clVector3Dd Maxp = max3(vls.C(Bead_BL(i,0)),vls.C(Bead_BL(i,1)),vls.C(Bead_BL(i,2)));
				clVector3Dd Minp = min3(vls.C(Bead_BL(i,0)),vls.C(Bead_BL(i,1)),vls.C(Bead_BL(i,2)));
				Maxp += 0.5;
				Minp += 0.5;

				clVector3Di Maxi = Maxp.floor();
				clVector3Di Mini = Minp.floor();

				Maxi += 1;
				Mini -= 1;

				if(Mini.y >= c.lb.nY)
					continue;

				if(Maxi.y < 0)
					continue;

				if(Mini.y < 0)
					Mini.y = 0;

				if(Maxi.y >= c.lb.nY)
					Maxi.y = c.lb.nY-1;

				if(Mini.z >= c.lb.nZ)
					Mini.z = c.lb.nZ-1;

				if(Mini.z < 0)
					Mini.z = 0;

				for(int ii = Mini.x; ii <= Maxi.x; ii++)
					for(int jj = Mini.y; jj <= Maxi.y; jj++)
						for(int kk = Mini.z; kk <= Maxi.z; kk++)
							B_Table(MODFAST(ii,c.lb.nX),MODFAST(jj,c.lb.nY),kk).appendItem(j);
			}
		}





		void FindMinMax()
		{
#ifdef TR_RIDGES
	#ifdef TR_BEADS
			Max_vec = clVector3Dd((double)c.lb.nX - 0.5, (double)c.lb.nY - 1., (double)c.lb.nZ);
			Min_vec = clVector3Dd(-0.5, 0., -0.5);
	#else
			Max_vec = clVector3Dd((double)c.lb.nX - 0.5, (double)c.lb.nY - 1.,(double)c.lb.nZ);
			Min_vec = clVector3Dd(-0.5, 0., 0.);
	#endif
#else
			Max_vec = clVector3Dd((double)c.lb.nX - 0.5, (double)c.lb.nY - 1., (double)c.lb.nZ - 1.);
			Min_vec = clVector3Dd(-0.5, 0., -0.5);
#endif
		}

		clVector3Dd min3(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
		{
			double x = v0.x, y = v0.y, z = v0.z;
			if(v1.x < x) x = v1.x;
			if(v2.x < x) x = v2.x;
			if(v1.y < y) y = v1.y;
			if(v2.y < y) y = v2.y;
			if(v1.z < z) z = v1.z;
			if(v2.z < z) z = v2.z;
			return clVector3Dd(x,y,z);
		}

		clVector3Dd max3(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
		{
			double x = v0.x, y = v0.y, z = v0.z;
			if(v1.x > x) x = v1.x;
			if(v2.x > x) x = v2.x;
			if(v1.y > y) y = v1.y;
			if(v2.y > y) y = v2.y;
			if(v1.z > z) z = v1.z;
			if(v2.z > z) z = v2.z;
			return clVector3Dd(x, y, z);
		}

#ifdef TR_RIDGES
	#ifdef TR_BEADS
		void update()
		{
			int i;

#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{								
				clVector3Dd C2_C = vtr.C.s(i);
				
				if(C2_C > Min_vec && C2_C < Max_vec)
					continue;

				if(C2_C.z >= Max_vec.z)
					update_ridges(i,C2_C);
				else
					update_beads(i, C2_C);
			}
			vtr.C.update();
		}
			
	#else			
		void update()
		{
			int i;

#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{								
				clVector3Dd C2_C = vtr.C.s(i);			//position after advection and diffusion		
				if(C2_C > Min_vec && C2_C < Max_vec)
					continue;
				update_ridges(i, C2_C);
			}

			vtr.C.update();

		}
	#endif
#else

		void update()
		{
			int i;
#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{								
				clVector3Dd C2_C = vtr.C.s(i);
				if(C2_C > Min_vec && C2_C < Max_vec)
					continue;
				update_beads(i, C2_C);
			}
			vtr.C.update();
		}

#endif


		void update_ridges(int i, clVector3Dd C2_C)
		{
			double nZm = (double)c.lb.nZ - 0.5;
			double nXm = (double)c.lb.nX - 0.5;
			double nYm = (double)c.lb.nY - 0.5;
			

			clVector3Dd dC = vtr.C.d(i);			//distance traveled since last time step (updated to distance traveled after reflection)
			clVector3Dd C1 = C2_C - dC;				//location at previous timestep (updated to location where each reflection occurs)
			clVector3Dd dC_A = vtr.V(i)*c.t.dTtr;	//advection distplacement
			clVector3Dd dC_D =  dC - dC_A;			//diffusion displacement
				
			double dist;							//||dCc||/||dC|| where dCc is the distance between  
			clVector3Di ind;						//index of particle for use in Tables
			clVector3Dd vCcut, vN, vCvel;			//Intersection pt, Normal vector, Velocity of surf at vCcut

			double dist_use = 100;

			clTensorS3Dd R, U;						//Reflection Tensor and Unit Tensor
			U.unit();							
				
			int bl;									//boundary link index (value stored in each linked list element
			int bl_use;
			clVector3Dd vCcut_use, vN_use;

			clNode<int> *tempt;

			int Eflag = 0;							//Exit flag; when flag = 1 the particle has not crossed a surface
			int Aflag = 0;							//if Aflag = 0 advection velocity needs to be updated
			int Tflag = 0;


			while(Eflag == 0)
			{
				Eflag = 1;
				
				if(C2_C.x > nXm)
				{
					vtr.C.p(i).x++;
					C2_C.x -= (double)c.lb.nX;
					C1.x = C2_C.x - dC.x;
				}
				else if(C2_C.x < -0.5)
				{
					vtr.C.p(i).x--;
					C2_C.x += (double)c.lb.nX;
					C1.x = C2_C.x - dC.x;
				}

				if(C2_C.z < -0.5)
				{
					dist = (-0.5 - C1.z)/dC.z;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 1;
					}
				}

				if(C2_C.y < -0.5)
				{
					dist = (-0.5 - C1.y)/dC.y;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 2;
					}
				}

				if(C2_C.z > nZm)
				{
					dist = (nZm - C1.z)/dC.z;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 3;
					}
				}

				if(C2_C.y > nYm)
				{
					dist = (nYm - C1.y)/dC.y;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 4;
					}
				}
	
				clVector3Dd C2_Cs = C2_C+0.5;
				ind = C2_Cs.floor();
				ind &= c.nn-1;
				tempt = T_Table(ind.x,ind.y,ind.z).firstlocal();

				while(tempt != NULL)
				{
					bl = tempt->getData();
					if(bcFindIntersectionLinePlane(vCcut,dist,C1,dC,vls.C(c.ls.BL(bl,0)),vls.C(c.ls.BL(bl,1)),vls.C(c.ls.BL(bl,2)),vN) == TRUE)
					{
						if(dist < dist_use)
						{
							dist_use = dist;
							Tflag = 5+Aflag;
							bl_use = bl;
							vCcut_use = vCcut;
							vN_use = vN;
						}
					}

					tempt = tempt->getLink();
				}
				switch(Tflag)
				{
				case 1:
					 {
						C1 += dC*dist_use*0.999;
						C2_C.z = -0.5 - (1.-dist_use)*dC.z;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Z_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					 }
				case 2:
					 {
						C1 += dC*dist_use*0.999;
						C2_C.y = -0.5 - (1.-dist_use)*dC.y;
						dC = C2_C-C1;


						Aflag = 1;
						el(i) ^= Y_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					 }
				case 3:
					{
						C1 += dC*dist_use*0.999;
						C2_C.z = nZm - (1.-dist_use)*dC.z;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Z_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					}
				case 4:
					{
						C1 += dC*dist_use*0.999;
						C2_C.y = nYm - (1.-dist_use)*dC.y;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Y_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					}
				case 5:
					 {
						vCvel = vtr.V(i) / 2.;
						dC_A = vtr.V(i)*c.t.dTtr;
						dC = dC_A + dC_D;
						C2_C = C1 + dC;
						

						Aflag = 1;
						Eflag = 0;
						el(i) ^= Surf_Bool(bl_use);
						Tflag = 0;
						break;
					 }
				case 6:
					 {
						C1 = vCcut_use + vN_use*0.01;
						R = U - clTensorS3Dd(vN_use)*2.;
						dC = R * (C2_C - C1);
						C2_C = C1 + dC;

						Eflag = 0;
						el(i) ^= Surf_Bool(bl_use);
						Tflag = 0;
						break;
					}
				}
				dist_use = 100.;
			}

			vtr.C(i).x = (double)(vtr.C.p(i).x*c.lb.nX)+C2_C.x;
			vtr.C(i).y = C2_C.y;
			vtr.C(i).z = C2_C.z;
		}

		void update_beads(int i, clVector3Dd C2_C)
		{
			double nZm = (double)c.lb.nZ - 0.5;
			double nXm = (double)c.lb.nX - 0.5;
			double nYm = (double)c.lb.nY - 0.5;

			clVector3Dd dC = vtr.C.d(i);			//distance traveled since last time step (updated to distance traveled after reflection)
			clVector3Dd C1 = C2_C - dC;				//location at previous timestep (updated to location where each reflection occurs)
			clVector3Dd dC_A = vtr.V(i)*c.t.dTtr;	//advection distplacement
			clVector3Dd dC_D =  dC - dC_A;			//diffusion displacement
				
			double dist;							//||dCc||/||dC|| where dCc is the distance between  
			clVector3Di ind;						//index of particle for use in Tables
			clVector3Dd vCcut, vN, vCvel;			//Intersection pt, Normal vector, Velocity of surf at vCcut

			double dist_use = 100;

			clTensorS3Dd R, U;						//Reflection Tensor and Unit Tensor
			U.unit();							
				
			int bl;									//boundary link index (value stored in each linked list element
			int bl_use;
			clVector3Dd vCcut_use, vN_use;

			clNode<int> *tempt;

			int Eflag = 0;							//Exit flag; when flag = 1 the particle has not crossed a surface
			int Aflag = 0;							//if Aflag = 0 advection velocity needs to be updated
			int Tflag = 0;


			while(Eflag == 0)
			{
				Eflag = 1;
				
				if(C2_C.x > nXm)
				{
					vtr.C.p(i).x++;
					C2_C.x -= (double)c.lb.nX;
					C1.x = C2_C.x - dC.x;
				}
				else if(C2_C.x < -0.5)
				{
					vtr.C.p(i).x--;
					C2_C.x += (double)c.lb.nX;
					C1.x = C2_C.x - dC.x;
				}

				if(C2_C.z < -0.5)
				{
					dist = (-0.5 - C1.z)/dC.z;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 1;
					}
				}

				if(C2_C.y < -0.5)
				{
					dist = (-0.5 - C1.y)/dC.y;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 2;
					}
				}

				if(C2_C.z > nZm)
				{
					dist = (nZm - C1.z)/dC.z;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 3;
					}
				}

				if(C2_C.y > nYm)
				{
					dist = (nYm - C1.y)/dC.y;
					if(dist < dist_use)
					{
						dist_use = dist;
						Tflag = 4;
					}
				}
	
				clVector3Dd C2_Cs = C2_C+0.5;
				ind = C2_Cs.floor();
				ind &= c.nn-1;
				tempt = B_Table(ind.x,ind.y,ind.z).firstlocal();

				while(tempt != NULL)
				{
					bl = tempt->getData();
					if(bcFindIntersectionLinePlane(vCcut,dist,C1,dC,vls.C(c.ls.BL(bl,0)),vls.C(c.ls.BL(bl,1)),vls.C(c.ls.BL(bl,2)),vN) == TRUE)
					{
						if(dist < dist_use)
						{
							dist_use = dist;
							Tflag = 5+Aflag;
							bl_use = bl;
							vCcut_use = vCcut;
							vN_use = vN;
						}
					}

					tempt = tempt->getLink();
				}

				switch(Tflag)
				{
				case 1:
					 {
						C1 += dC*dist_use*0.999;
						C2_C.z = -0.5 - (1.-dist_use)*dC.z;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Z_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					 }
				case 2:
					 {
						C1 += dC*dist_use*0.999;
						C2_C.y = -0.5 - (1.-dist_use)*dC.y;
						dC = C2_C-C1;


						Aflag = 1;
						el(i) ^= Y_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					 }
				case 3:
					{
						C1 += dC*dist_use*0.999;
						C2_C.z = nZm - (1.-dist_use)*dC.z;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Z_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					}
				case 4:
					{
						C1 += dC*dist_use*0.999;
						C2_C.y = nYm - (1.-dist_use)*dC.y;
						dC = C2_C-C1;

						Aflag = 1;
						el(i) ^= Y_SURF;
						Eflag = 0;
						Tflag = 0;
						break;
					}
				case 5:
					 {
						vCvel = bcGetBoundaryVelocity(vCcut_use, vls.C(c.ls.BL(bl_use,0)), vls.C(c.ls.BL(bl_use,1)), vls.C(c.ls.BL(bl_use,2)), bl_use);
						vCvel = (vtr.V(i) + vCvel) / 2.;
						dC_A = vCvel*c.t.dTtr;
						dC = dC_A + dC_D;
						C2_C = C1 + dC;
						

						Aflag = 1;
						Eflag = 0;
						el(i) ^= Surf_Bool(bl_use);
						Tflag = 0;
						break;
					 }
				case 6:
					 {
						C1 = vCcut_use + vN_use*0.001;
						R = U - clTensorS3Dd(vN_use)*2.;
						dC = R * (C2_C - C1);
						C2_C = C1 + dC;

						Eflag = 0;
						el(i) ^= Surf_Bool(bl_use);
						Tflag = 0;
						break;
					}
				}
				dist_use = 100.;
			}
				
			vtr.C(i).x = (double)(vtr.C.p(i).x*c.lb.nX)+C2_C.x;
			vtr.C(i).y = C2_C.y;
			vtr.C(i).z = C2_C.z;
		}	
		
		BOOL testInside(clVector3Dd vd, clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
		{
			clVector3Dd v10 = v1 - v0, v20 = v2 - v0, v21 = v2 - v1, vd0 = vd - v0, vd1 = vd - v1;
			clVector3Dd stot = v10 ^ v20, s0 = vd0 ^ v10, s1 = vd0 ^ v20, s2 = vd1 ^ v21;

			if(fabs(stot.r() - s0.r() - s1.r() - s2.r()) < c.eps)
				return TRUE;
			return FALSE;
		}

		BOOL bcFindIntersectionLinePlane(clVector3Dd &vC, double &dist, clVector3Dd vL0,
			clVector3Dd vLd, clVector3Dd vP0, clVector3Dd vP1, clVector3Dd vP2, clVector3Dd &vN)
		{
			clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
			clVector3Dd vNm = vP10 ^ vP20;
			double den = vNm * vLd;
			if(den == 0.)
				return FALSE;
			vN = vNm.n();

			dist = vNm * (vP0 - vL0) / den;
			vC = vL0 + vLd * dist;
			if(dist >= 0. && dist <= 1.)
				if(testInside(vC,vP0,vP1,vP2) == TRUE)
					return TRUE;

			return FALSE;

		}

		clVector3Dd bcGetBoundaryVelocity(clVector3Dd vC, clVector3Dd vP0, clVector3Dd vP1, clVector3Dd vP2, int bl)
		{
			clVector3Dd vP12 = vP2-vP1, vP0C = vC-vP0, vP10 = vP1-vP0;
			clVector3Dd numer = vP10 ^ vP12, den = vP0C ^ vP12;
			double dist1 = numer.r()/den.r();
			clVector3Dd vI = vP0 + vP0C*dist1, vP1I = vI-vP1;
			double dist2 = vP1I.r();
			clVector3Dd vIv = vls.V(c.ls.BL(bl,1)) - (vls.V(c.ls.BL(bl,1)) - vls.V(c.ls.BL(bl,2))) * dist2/vP12.r();

			return vls.V(c.ls.BL(bl,0)) -  (vls.V(c.ls.BL(bl,0)) - vIv) * vP0C.r()/dist1;
		}

		void update_bead_centers()
		{
			clArray1Dv3d Pos;
			clArray1Di count;
			Pos.zeros(N_BEAD);
			count.zeros(N_BEAD);
			double nXm = (double)c.lb.nX - 0.5;
			Min_vec.z = -0.5;

			for(int i = 0; i < Bead_nBL; i++)
			{
				clVector3Dd C0 = vls.C(Bead_BL(i,0));
				clVector3Dd C1 = vls.C(Bead_BL(i,1));
				clVector3Dd C2 = vls.C(Bead_BL(i,2));

				
				Pos(Bead_Num(i)) += C0 + C1 + C2;
				clVector3Dd Maxval = max3(C0,C1,C2);
				if(Maxval.z > Min_vec.z)
						Min_vec.z = Maxval.z;

				count(Bead_Num(i)) += 3;
			}

			Min_vec.z += 0.2;

			for(int i = 0; i < N_BEAD; i++)
			{
				clVector3Dd BC = Pos(i)/((double)count(i));
				if(BC.x > nXm)
				{
					BC.x = fmod(BC.x+0.5,(double)c.lb.nX) - 0.5;
				}
				Bead_center(i) = BC;
			}
		}


		void  Test_Inside_Beads()
		{
			update_bead_centers();
			double nXm = (double)c.lb.nX-0.5;
			
			int i;
#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
			for(i = 0; i < c.tr.nN; i++)
			{							
				clVector3Dd C1 = vtr.C.s(i);  
				clVector3Dd C2;

				if(C1.z > Min_vec.z)
					continue;
				
				if(C1.x >= (double)c.lb.nX/2.)
					C2 = clVector3Dd(C1.x - (double)c.lb.nX,C1.y,C1.z);
				else
					C2 = clVector3Dd(C1.x + (double)c.lb.nX,C1.y,C1.z);

				
				int bl;
				clNode<int> *tempt;

				clVector3Dd vLd;
				double dist;



				double B_radius = 1.01*R0_NUMBER;

				for(int j = 0; j < N_BEAD; j++)
				{
					clVector3Dd Dist1 = C1 - Bead_center(j);
					clVector3Dd Dist2 = C2 - Bead_center(j);
					if(Dist1.r() <= B_radius || Dist2.r() <= B_radius)
					{
						int Eflag = 0;
						
						double dist_use = 100.;
						clVector3Dd vLd_use;

						while(Eflag == 0)
						{
							Eflag = 1;
							clVector3Dd C1i = C1 + 0.5;
							clVector3Di ind = C1i.floor();


							tempt = B_Table(ind.x,ind.y,ind.z).firstlocal();
							while(tempt != NULL)
							{
								bl = tempt->getData();
								if(TestParInBead(C1,vls.C(c.ls.BL(bl,0)),vls.C(c.ls.BL(bl,1)),vls.C(c.ls.BL(bl,2)),vLd,dist) == TRUE)
								{
									if(dist < dist_use)
									{
										vLd_use = vLd;
										dist_use = dist;	
										Eflag = 2;
									}
									tempt = tempt->getLink();
								}
								else
								{
									tempt = tempt->getLink();
								}
							}

							if(Eflag == 2)
							{
								C1 += vLd_use*2.;

								if(C1.z < -0.5)
								{
									C1.z = -1. - C1.z;
									el(i) ^= Z_SURF;
									Eflag = 0;
								}

							}

							if(C1.x > nXm)
							{
								vtr.C.p(i).x++;
								C1.x -= (double)c.lb.nX;
							}
							if(C1.x < -0.5)
							{
								vtr.C.p(i).x--;
								C1.x += (double)c.lb.nX;
							}

						}
						vtr.C(i).x = (double)(vtr.C.p(i).x*c.lb.nX)+C1.x;
						vtr.C(i).y = C1.y;
						vtr.C(i).z = C1.z;
					}
				}
			}
		}



		BOOL TestParInBead(clVector3Dd vL0, clVector3Dd vP0, clVector3Dd vP1, 
			clVector3Dd vP2, clVector3Dd &vLd, double &dist)
		{
			clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
			clVector3Dd vNm = vP10 ^ vP20;
			clVector3Dd vN = vNm.n();
			clVector3Dd vPL = vP0 - vL0;
			dist = vPL*vN;
			vLd = vN*dist;

			if(dist >= 0. && dist < R0_NUMBER)
			{
			
				clVector3Dd vC = vL0 + vLd;
				if(testInside(vC,vP0,vP1,vP2) == TRUE)
					return TRUE;
			}

			return FALSE;

		}

		void save2file()
		{
			save("tra");
		};
	} A;

	class clU: public clArray3Dv3d
	{
	public:
		void ini()
		{
			clArray3Dv3d::zeros(c.lb.nX+2, c.lb.nY+2, c.lb.nZ+2);
		}
		void update()
		{
			updateLB();
			updateLS();
		}

		void updateLB()
		{
			for(int i = 0; i < c.lb.nX; i++)
				for(int j = 0; j < c.lb.nY; j++)
					for(int m = 0; m < c.lb.nZ; m++)
						if(vlb.M(i,j,m) != LB_SOLID)
							el(i+1, j+1, m+1) = vlb.J(i, j, m) / vlb.Ro(i, j, m);

			if(c.lb.bcPeriodicX == TRUE) bcPeriodicLBx();
			if(c.lb.bcPeriodicY == TRUE) bcPeriodicLBy();
			if(c.lb.bcPeriodicZ == TRUE) bcPeriodicLBz();

			if(c.lb.bcFreeX0 == TRUE) bcFreeFlowLBx0();
			if(c.lb.bcFreeX1 == TRUE) bcFreeFlowLBx1();
			if(c.lb.bcFreeY0 == TRUE) bcFreeFlowLBy0();
			if(c.lb.bcFreeY1 == TRUE) bcFreeFlowLBy1();
			if(c.lb.bcFreeZ0 == TRUE) bcFreeFlowLBz0();
			if(c.lb.bcFreeZ1 == TRUE) bcFreeFlowLBz1();

			if(c.lb.bcWallX0 == TRUE) bcWallLBx0();
			if(c.lb.bcWallX1 == TRUE) bcWallLBx1();
			if(c.lb.bcWallY0 == TRUE) bcWallLBy0();
			if(c.lb.bcWallY1 == TRUE) bcWallLBy1();
			if(c.lb.bcWallZ0 == TRUE) bcWallLBz0();
			if(c.lb.bcWallZ1 == TRUE) bcWallLBz1();


		}
		void bcPeriodicLBx()
		{
			const int nXm0 = c.lb.nX+1;
			const int nXm1 = c.lb.nX;
			for(int j = 0; j < c.lb.nY+2; j++)
				for(int z = 0; z < c.lb.nZ+2; z++)
				{
					el(nXm0,j,z) = el(1,j,z);
					el(0,j,z) = el(nXm1,j,z);
				}
		}
		void bcPeriodicLBy()
		{
			const int nYm0 = c.lb.nY+1;
			const int nYm1 = c.lb.nY;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int z = 0; z < c.lb.nZ+2; z++)
				{
					el(i,nYm0,z) = el(i,1,z);
					el(i,0,z) = el(i,nYm1,z);
				}
		}
		void bcPeriodicLBz()
		{
			const int nZm0 = c.lb.nZ+1;
			const int nZm1 = c.lb.nZ;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int j = 0; j < c.lb.nY+2; j++)
				{
					el(i,j,nZm0) = el(i,j,1);
					el(i,j,0) = el(i,j,nZm1);
				}
		}

		void bcFreeFlowLBx0()
		{
			for(int j = 0; j < c.lb.nY+2; j++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(0,j,z) = el(1,j,z);
		}
		void bcFreeFlowLBx1()
		{
			const int nXm0 = c.lb.nX+1;
			const int nXm1 = c.lb.nX;
			for(int j = 0; j < c.lb.nY+2; j++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(nXm0,j,z) = el(nXm1,j,z);
		}
		void bcFreeFlowLBy0()
		{
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(i,0,z) = el(i,1,z);
		}
		void bcFreeFlowLBy1()
		{
			const int nYm0 = c.lb.nY+1;
			const int nYm1 = c.lb.nY;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(i,nYm0,z) = el(i,nYm1,z);
		}
		void bcFreeFlowLBz0()
		{
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int j = 0; j < c.lb.nY+2; j++)
					el(i,j,0) = el(i,j,1);
		}
		void bcFreeFlowLBz1()
		{
			const int nZm0 = c.lb.nZ+1;
			const int nZm1 = c.lb.nZ;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int j = 0; j < c.lb.nY+2; j++)
					el(i,j,nZm0) = el(i,j,nZm1);
		}

		void bcWallLBx0()
		{
			for(int j = 0; j < c.lb.nY+2; j++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(0,j,z) = c.lb.UX0+c.lb.UX0-el(1,j,z);
		}
		void bcWallLBx1()
		{
			const int nXm0 = c.lb.nX+1;
			for(int j = 0; j < c.lb.nY+2; j++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(nXm0,j,z) = c.lb.UX1+c.lb.UX1-el(c.lb.nX,j,z);
		}
		void bcWallLBy0()
		{
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(i,0,z) = c.lb.UY0+c.lb.UY0-el(i,1,z);
		}
		void bcWallLBy1()
		{
			const int nYm0 = c.lb.nY+1;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int z = 0; z < c.lb.nZ+2; z++)
					el(i,nYm0,z) = c.lb.UY1+c.lb.UY1-el(i,c.lb.nY,z);
		}
		void bcWallLBz0()
		{
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int j = 0; j < c.lb.nY+2; j++)
					el(i,j,0) = c.lb.UZ0+c.lb.UZ0-el(i,j,1);
		}
		void bcWallLBz1()
		{
			const int nZm0 = c.lb.nZ+1;
			for(int i = 0; i < c.lb.nX+2; i++)
				for(int j = 0; j < c.lb.nY+2; j++)
					el(i,j,nZm0) = c.lb.UZ1+c.lb.UZ1-el(i,j,c.lb.nZ);
		}

		void updateLS()
		{
			for(int n = 0; n < vlb.M.nSf; n++)
			{
				if(vlb.M.test_stationary(n) == TRUE)
					continue;

				int n0 = vlb.M.Sf(n,0), n1 = vlb.M.Sf(n,1), n2 = vlb.M.Sf(n,2), n3 = vlb.M.Sf(n,3); 
				if(c.ls.N.NL(n0) == 0 || c.ls.N.NL(n1) == 0 || c.ls.N.NL(n2) == 0 || c.ls.N.NL(n3) == 0)
					continue;

				clVector3Di p0 = vls.C.p(n0), p1 = vls.C.p(n1), p2 = vls.C.p(n2), p3 = vls.C.p(n3);
				if(p0 == p1 && p1 == p2 && p2 == p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					continue;
				}
				if(p0 != p1 && p1 == p2 && p2 == p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n1);
					continue;
				}
				if(p0 == p1 && p1 != p2 && p2 == p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n2);
					continue;
				}
				if(p0 == p1 && p1 == p2 && p2 != p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 != p1 && p1 != p2 && p2 == p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n1);
					updateLS(n0, n1, n2, n3, n2);
					continue;
				}
				if(p0 != p1 && p1 == p2 && p2 != p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n1);
					updateLS(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 == p1 && p1 != p2 && p2 != p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n2);
					updateLS(n0, n1, n2, n3, n3);
					continue;
				}
				if(p0 != p1 && p1 != p2 && p2 != p3)
				{
					updateLS(n0, n1, n2, n3, n0);
					updateLS(n0, n1, n2, n3, n1);
					updateLS(n0, n1, n2, n3, n2);
					updateLS(n0, n1, n2, n3, n3);
					continue;
				}

			}
		};
		void updateLS(int n0, int n1, int n2, int n3, int p)
		{
			clVector3Dd c0 = vls.C.s(n0,p),  c1 = vls.C.s(n1,p), c2 = vls.C.s(n2,p), c3 = vls.C.s(n3,p);
			clVector3Di imin = clVector3Di(clVector3Dd(vlb.M.min(c0, c1, c2, c3)).ceil()).max(clVector3Di(0,0,0));
			clVector3Di imax = clVector3Di(clVector3Dd(vlb.M.max(c0, c1, c2, c3)).floor() + 1).min(c.lb.nn);
			clVector3Di ii;
			for(ii.x = imin.x; ii.x < imax.x; ii.x++)
				for(ii.y = imin.y; ii.y < imax.y; ii.y++)
					for(ii.z = imin.z; ii.z < imax.z; ii.z++)
					{
						clVector3Dd cc = (clVector3Dd)ii;
						if(vlb.M.test_inside(cc, c0, c1, c2, c3) == TRUE)
						{
							clVector3Dd c13 = c1 - c3, c23 = c2 - c3;
							clVector3Dd cc0 = cc - c0;
							clVector3Dd cN0 = c13 ^ c23;
							double den0 = cN0 * cc0;
							if(cc0.abs() < c.eps || fabs(den0) < c.eps)
							{
								el(ii+1) = vls.V(n0);
								continue;
							}

							double dist0 = cN0 * (c3 - c0) / den0;
							clVector3Dd cp = c0 + cc0 * dist0;

							clVector3Dd c03 = c0 - c3;
							clVector3Dd cp1 = cp - c1;
							clVector3Dd cN1 = c03 ^ c23;
							double den1 = cN1 * cp1;
							if(cp1.abs() < c.eps || fabs(den1) < c.eps)
							{
								el(ii+1) = vls.V(n1);
								continue;
							}

							double dist1 = cN1 * (c3 - c1) / den1;
							clVector3Dd cp23 = c1 + cp1 * dist1;

							double d3p = clVector3Dd(c3 - cp23).abs();
							double d2p = clVector3Dd(c2 - cp23).abs();
							//double d23 = clVector3Dd(c2 - c3).abs();
							clVector3Dd vp23 = (vls.V(n3) * d2p + vls.V(n2) * d3p) / (d2p + d3p);

							double d1p = clVector3Dd(c1 - cp).abs();
							double d23p = clVector3Dd(cp23 - cp).abs();
							clVector3Dd vp = (vp23 * d1p + vls.V(n1) * d23p) / (d1p + d23p);

							double dcp = clVector3Dd(cc - cp).abs();
							double d0c = clVector3Dd(c0 - cc).abs();
							clVector3Dd vc = (vp * dcp + vls.V(n0) * d0c) / (dcp + d0c);

							//double Lvc = vc.abs();
							//double Lv0 = vls.V(n0).abs(), Lv1 = vls.V(n1).abs(), Lv2 = vls.V(n2).abs(), Lv3 = vls.V(n3).abs();
							//if(Lvc > Lv0 && Lvc > Lv1 && Lvc > Lv2 && Lvc > Lv3)
							//	printf("Velocity interpolation error!!!!!!!");

							el(ii+1) = vc;
						}
					}

		};

		void save2file()
		{
			savetxt("tru");
		};
	} U;
};

#endif // !defined(AFX_CLVARIABLES_H__INCLUDED_)
