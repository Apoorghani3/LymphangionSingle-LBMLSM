// clConstants.h: interface for the clConstants class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLCONSTANTS_H__INCLUDED_)
#define AFX_CLCONSTANTS_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <stdio.h>
#include "StdAfx.h"
//#include "clArray1D.h"
//#include "clArray2D.h"
//#include "clOutput.h"

class clParemeters
{
protected:
	FILE *outstream;
	BOOL create_out_parameter()
	{
		return create_out_parameter("par.m");
	}

	BOOL create_out_parameter(const char * FileName)
	{
		if((outstream = fopen(FileName, "w")) == NULL)
		{
			printf( "\nUnable to open the output file!\n" );
			return FALSE;
		}
		return TRUE;
	}

	BOOL close_out_parameter()
	{
		if(outstream != NULL)
			fclose(outstream);
		return TRUE;
	}

	void saveparameter(const char *Name, const char *Value)
	{
		fprintf(outstream, "%s = '%s';\n", Name, Value);
	}

	void saveparameter_line()
	{
		fprintf(outstream, "\n");
	}

	void saveparameter(const char *Name, const int Value)
	{
		fprintf(outstream, "%s = %d;\n", Name, Value);
	}

	void saveparameter(const char *Name, const float Value)
	{
		fprintf(outstream, "%s = %f;\n", Name, Value);
	}

	void saveparameter(const char *Name, const double Value)
	{
		fprintf(outstream, "%s = %lg;\n", Name, Value);
	}

	BOOL getparameter(const char * FileName, const char * Format, double *Buf)
	{
		FILE *stream;
		if((stream = fopen(FileName, "r")) == NULL)
			return FALSE;

		long seekpos = 0L;

		fseek(stream, seekpos, SEEK_SET);
		while(TRUE)
		{
			if(fscanf(stream, Format, Buf) != 0)
				break;
			fseek(stream, ++seekpos, SEEK_SET);
		}

		fclose(stream);
		return TRUE;
	}

	BOOL getparameter(const char * FileName, const char * Format, int *Buf)
	{
		FILE *stream;
		if((stream = fopen(FileName, "r")) == NULL)
			return FALSE;

		long seekpos = 0L;

		fseek(stream, seekpos, SEEK_SET);
		while(TRUE)
		{
			if(fscanf(stream, Format, Buf) != 0)
				break;
			fseek(stream, ++seekpos, SEEK_SET);
		}
		fclose(stream);
		return TRUE;
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
	}

};

class clConstants: protected clParemeters  
{
	BOOL loadparameters(const char * FileName)
	{
		getparameter(FileName, "c_Re1 = %lg", &Re1);
		getparameter(FileName, "c_Pe1 = %lg", &Pe1);
		getparameter(FileName, "c_Gx = %lg", &G.x);
		getparameter(FileName, "c_Gy = %lg", &G.y);
		getparameter(FileName, "c_Gz = %lg", &G.z);
		getparameter(FileName, "c_L0 = %lg", &L0);
		getparameter(FileName, "c_Ro0 = %lg", &Ro0);
		getparameter(FileName, "c_nu0 = %lg", &nu0);
		getparameter(FileName, "c_eps = %lg", &eps);

		sL = L0;
		sT = L0 * L0 / nu0;
		sM = Ro0 * L0 * L0 * L0;
		sN = sM * sL / (sT * sT);

		return TRUE;
	}
public:

	const static double Pi, R, Kb;
	static double Re1, Pe1, Angle;
	static double T0, nu0, Ro0, L0;
	static double sL, sT, sM, sN;
	static int nX, nY, nZ;
	static clVector3Dd G;
	static double eps;
	static clVector3Di nn;



	virtual ~clConstants(){};

	void ini()
	{
		srand(1);

		loadparameters("load" SLASH "par.m");
		//if(loadparameters("load" SLASH "par.m") == FALSE || nX == 0 || nY == 0)
		//	iniXY();

		t.ini();
		ls.ini();
		lb.ini();
		tr.ini();

		nX = lb.nX;
		nY = lb.nY;
		nZ = lb.nZ;

		nn = clVector3Di(nX, nY, nZ);
	};

	void iniXY()
	{
		//int x, y;
		//x = counttxt("load" SLASH "cx.txt");
		//if(x > 0) nX = x;

		//y = counttxt("load" SLASH "cy.txt");
		//if(y > 0) nY = y;
	};
	void update()
	{
		t.update();
		lb.update();
		ls.update();
		tr.update();
	}	
	
	void signal()
	{
		static double NextCurrentSignalTime = t.dTlb * (double)t.CurrentSaveTimeStep;
		static double NextDumpSignalTime = t.DumpTimeStep;

		if(NextCurrentSignalTime < t.TimeN)
		{
			stepSignal();
			NextCurrentSignalTime += t.dTlb * t.CurrentSaveTimeStep;
		}

		if(NextDumpSignalTime < t.TimeN)
		{
			fieldSignal();
			if(t.Time < t.DumpStartTime)
				NextDumpSignalTime += t.DumpTimeStep;
			else
				NextDumpSignalTime += t.DumpTimeStepSave;
		}
	}
	
	void stepSignal();
	void fieldSignal();
	void save2file();

	class clTime: protected clParemeters  
	{
		BOOL loadparameters(const char * FileName)
		{
			int fl = 0;
			getparameter(FileName, "c_t_flContinue = %d", &fl);
			getparameter(FileName, "c_t_flOutputDump = %d", &flOutputDump);
			getparameter(FileName, "c_t_dTlb = %lg", &dTlb);
			getparameter(FileName, "c_t_dTls = %lg", &dTls);
			getparameter(FileName, "c_t_dTtr = %lg", &dTtr);
			getparameter(FileName, "c_t_Tstop = %lg", &StopTime);
			getparameter(FileName, "c_t_Ttotal = %lg", &StartTime);
			getparameter(FileName, "c_t_Tdumpstart = %lg", &DumpStartTime);
			getparameter(FileName, "c_t_Tdumpstep = %lg", &DumpTimeStep);
			getparameter(FileName, "c_t_Tdumpstepsave = %lg", &DumpTimeStepSave);
			getparameter(FileName, "c_t_Tcurrentsavestart = %d", &CurrentSaveStartTime);
			getparameter(FileName, "c_t_Tcurrentsavestep = %d", &CurrentSaveTimeStep);
	
			getparameter(FileName, "c_t_oscAmpl1 = %lg", &oscAmpl1);
			getparameter(FileName, "c_t_oscPeriod1 = %lg", &oscPeriod1);
			getparameter(FileName, "c_t_oscAmpl2 = %lg", &oscAmpl2);
			getparameter(FileName, "c_t_oscPeriod2 = %lg", &oscPeriod2);
			getparameter(FileName, "c_t_oscRatio = %lg", &oscRatio);

			if((fl == 0 && StartTime == 0) || fl == 2)
				flContinue = 0;
			else
				flContinue = 1;


			lsSteps = (int)(dTlb / dTls);

			return TRUE;
		}

	public:
		static double Time, TimeN, dTlb, dTls, dTtr, StopTime, StartTime;
		static double DumpTimeStep, DumpTimeStepSave, DumpStartTime; 
		static int CurrentSaveTimeStep, CurrentSaveStartTime, DisplaySignal;
		static int lsSteps, trSteps;
		static double oscAmpl1, oscPeriod1, oscAmpl2, oscPeriod2; 
		static double oscRatio; 
		static BOOL flContinue;
		static BOOL flOutputDump;
		BOOL NoUpdateSignal;

		void ini()
		{
			loadparameters("load" SLASH "par0.m");
			loadparameters("load" SLASH "par.m");
			NoUpdateSignal = FALSE;
			TimeN = dTls;
		};
		void noupdate()
		{
			NoUpdateSignal = TRUE;
		};
		void update();

		double oscillator1()
		{
			return oscillator1(0.);
		}

		double oscillator1sin()
		{
			//return 0.;
			return oscillator1(0.);
		}
		double oscillator1cos()
		{
			//return	1.;
			return oscillator1cos(0.);
		}

		double oscillator_ratio()
		{
			double T = (c.t.Time + c.t.StartTime) / oscPeriod1;
			double T1 = 2. * oscPeriod1 / (1. + 1. / oscRatio);
			double T2 = 2. * oscPeriod1 / (1. + oscRatio);
			T = 2. * (T - floor(T)) * oscPeriod1;
			if(T < T1)
				return oscAmpl1 * cos(c.Pi * T / T1);

			return -oscAmpl1 * cos(c.Pi * (T - T1) / T2);
		}

		//double oscillator1(double phase)
		//{
		//	return oscAmpl1 * sin(2. * c.Pi * (c.t.Time + c.t.StartTime) / oscPeriod1 + phase);
		//}
#ifdef SQUARE_WAVE
		double oscillator1(double phase)
		{
			return 2*oscAmpl1 * atan(sin(2. * c.Pi * (MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD)) / oscPeriod1 + phase)/SQ_BETA)/c.Pi;
		}
		double oscillator2()
		{
			return -2*oscAmpl2 * atan(sin(2. * c.Pi * (MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD) - 0.5 * oscPeriod1) / oscPeriod2)/SQ_BETA)/c.Pi;
		}
#else
		double oscillator1(double phase)
		{
			return oscAmpl1 * sin(2. * c.Pi * (MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD)) / oscPeriod1 + phase);
		}
		//double oscillator2()
		//{
		//	return oscAmpl2 * sin(2. * c.Pi * (c.t.Time + c.t.StartTime) / oscPeriod2);
		//}
		double oscillator2()
		{
			return -oscAmpl2 * sin(2. * c.Pi * (MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD) - 0.5 * oscPeriod1) / oscPeriod2);
		}
#endif // SQUARE_WAVE

#ifdef RAMP
		double ramping()
		{			
			if(fmod(c.t.Time,RAMP_PERIOD)< RAMP_PERIOD*0.1)
			{
				return 1.;
			}
			else if(fmod(c.t.Time,RAMP_PERIOD) >= RAMP_PERIOD*0.1 && fmod(c.t.Time,RAMP_PERIOD) < RAMP_PERIOD*0.4)
			{
				return (RAMP_PERIOD*0.25-fmod(c.t.Time,RAMP_PERIOD))*1/(RAMP_PERIOD*0.15);
			}
			else if(fmod(c.t.Time,RAMP_PERIOD) >= RAMP_PERIOD*0.4 && fmod(c.t.Time,RAMP_PERIOD) < RAMP_PERIOD*0.6)
			{
				return -1. ;
			}
			else if(fmod(c.t.Time,RAMP_PERIOD) >= RAMP_PERIOD*0.6 && fmod(c.t.Time,RAMP_PERIOD) < RAMP_PERIOD*0.9)
			{
				return -(RAMP_PERIOD*0.75-fmod(c.t.Time,RAMP_PERIOD))*1/(RAMP_PERIOD*0.15);
			}
			else
			{
				return 1. ;
			}
		}
#endif

		double oscillator1cos(double phase)
		{
			return oscAmpl1 * cos(2. * c.Pi * (c.t.Time + c.t.StartTime) / oscPeriod1 + phase);
		}
		double random_normal()
		{
			return (double) rand() / (double) RAND_MAX;
		}
		void save2file()
		{
		};
	} t;

	class clLS: protected clParemeters  
	{
		BOOL loadparameters(const char * FileName)
		{
			getparameter(FileName, "c_ls_nN = %d", &nN);
			getparameter(FileName, "c_ls_nL = %d", &nL);
			getparameter(FileName, "c_ls_nAN = %d", &nAN);
//			getparameter(FileName, "c_ls_Lb = %lg", &Lb);
			getparameter(FileName, "c_ls_nPupd = %lg", &nPupd);
			getparameter(FileName, "c_ls_FwDe0 = %lg", &FwDe0);
			getparameter(FileName, "c_ls_FwDe1 = %lg", &FwDe1);
			getparameter(FileName, "c_ls_FwDe2 = %lg", &FwDe2);
			getparameter(FileName, "c_ls_FwDeR = %lg", &FwDeR);
			getparameter(FileName, "c_ls_FwRe = %lg", &FwRe);
			getparameter(FileName, "c_ls_FwKe = %lg", &FwKe);
			getparameter(FileName, "c_ls_FwDeC = %lg", &FwDeC);
			getparameter(FileName, "c_ls_FwReC = %lg", &FwReC);

			getparameter(FileName, "c_ls_Ds0 = %lg", &Ds0);
			getparameter(FileName, "c_ls_Ds1 = %lg", &Ds1);
			getparameter(FileName, "c_ls_Ds2 = %lg", &Ds2);

			getparameter(FileName, "c_ls_Symax0 = %lg", &Symax0);
			getparameter(FileName, "c_ls_Symin0 = %lg", &Symin0);
			getparameter(FileName, "c_ls_Symax1 = %lg", &Symax1);
			getparameter(FileName, "c_ls_Symin1 = %lg", &Symin1);
			getparameter(FileName, "c_ls_Symax2 = %lg", &Symax2);
			getparameter(FileName, "c_ls_Symin2 = %lg", &Symin2);

			getparameter(FileName, "c_ls_FwDeX0 = %lg", &FwDeX0);
			getparameter(FileName, "c_ls_FwDeRX0 = %lg", &FwDeRX0);
			getparameter(FileName, "c_ls_FwDeX1 = %lg", &FwDeX1);
			getparameter(FileName, "c_ls_FwDeRX1 = %lg", &FwDeRX1);
			getparameter(FileName, "c_ls_FwDeY0 = %lg", &FwDeY0);
			getparameter(FileName, "c_ls_FwDeRY0 = %lg", &FwDeRY0);
			getparameter(FileName, "c_ls_FwDeY1 = %lg", &FwDeY1);
			getparameter(FileName, "c_ls_FwDeRY1 = %lg", &FwDeRY1);
			getparameter(FileName, "c_ls_FwDeZ0 = %lg", &FwDeZ0);
			getparameter(FileName, "c_ls_FwDeRZ0 = %lg", &FwDeRZ0);
			getparameter(FileName, "c_ls_FwDeZ1 = %lg", &FwDeZ1);
			getparameter(FileName, "c_ls_FwDeRZ1 = %lg", &FwDeRZ1);

			getparameter(FileName, "c_ls_Cini0x = %lg", &Cini0.x);
			getparameter(FileName, "c_ls_Cini0y = %lg", &Cini0.y);
			getparameter(FileName, "c_ls_Cini0z = %lg", &Cini0.z);
			getparameter(FileName, "c_ls_Cini1x = %lg", &Cini1.x);
			getparameter(FileName, "c_ls_Cini1y = %lg", &Cini1.y);
			getparameter(FileName, "c_ls_Cini1z = %lg", &Cini1.z);
			getparameter(FileName, "c_ls_Cini2x = %lg", &Cini2.x);
			getparameter(FileName, "c_ls_Cini2y = %lg", &Cini2.y);
			getparameter(FileName, "c_ls_Cini2z = %lg", &Cini2.z);

			getparameter(FileName, "c_ls_Kini0 = %lg", &Kini0);
			getparameter(FileName, "c_ls_Kini1 = %lg", &Kini1);
			getparameter(FileName, "c_ls_Kini2 = %lg", &Kini2);

			getparameter(FileName, "c_ls_Mini0 = %lg", &Mini0);
			getparameter(FileName, "c_ls_Mini1 = %lg", &Mini1);
			getparameter(FileName, "c_ls_Mini2 = %lg", &Mini2);

			getparameter(FileName, "c_ls_L0ini0 = %lg", &L0ini0);
			getparameter(FileName, "c_ls_L0ini1 = %lg", &L0ini1);
			getparameter(FileName, "c_ls_L0ini2 = %lg", &L0ini2);

			getparameter(FileName, "c_ls_AKini0 = %lg", &AKini0);
			getparameter(FileName, "c_ls_AKini1 = %lg", &AKini1);
			getparameter(FileName, "c_ls_AKini2 = %lg", &AKini2);

			int n = counttxt("load" SLASH "lsc") / NUMBER_DIMENSIONS;
			if(n > 0) nN = n;

			return TRUE;
		}
	public:
		//static int nN, nL, nBN, nBL, nTa, nTp, nTr, nO, nPN, nBM, nOp, nNd; 
		static int nN, nL, nBN, nBL, nTa, nTr, nTs, nO, nPN, nBM, nOp, nNd, nAN; 
		static double Ds0, Ds1, Ds2;
//		static double Lb;
		static double FwDe0, FwDe1, FwDe2, FwDeR, FwRe, FwKe;
		static double FwDeC, FwReC;
		static double FwDeX0, FwDeRX0, FwDeX1, FwDeRX1;
		static double FwDeY0, FwDeRY0, FwDeY1, FwDeRY1;
		static double FwDeZ0, FwDeRZ0, FwDeZ1, FwDeRZ1;
		static double Kini0, Kini1, Kini2;
		static double AKini0, AKini1, AKini2;
		static double Mini0, Mini1, Mini2;
		static double L0ini0, L0ini1, L0ini2;
		static double Symin0, Symax0, Symin1, Symax1, Symin2, Symax2;
		static clVector3Dd Cini0, Cini1, Cini2;
		static int nBFmax;
		static int nPupd;
		static double valveMinX;  // First position where lso != 0
		static double valveMaxX;  // Last position where lso != 0
		static double waveLength; // Peristaltic wave wavelength = 2 * intervalve (lymphangion) distance

		clArray1Di BN;
		clArray1Di BS;
		clArray1Di BF;
		clArray2Di BL;
		clArray2Dv3d BLFb;
		clArray1Di BNind;
		clArray1Di O;
		clArray1Di Os;
		clArray1Di Op;
		clArray1Di Opstack;
		clArray1Di NO;
		clArray1Dc RO;
		clArray1Dc SO;
		clArray1Dv3c PO;
		clArray1Di T;
		clArray1Di Inout;
		clArray1Di Ta;
		//clArray1Di Tp;
		clArray1Di Tr;
		clArray1Di Ts;
		clArray1Dc Tno_c, Tno_f;
		clArray2Di PN;
		clArray1Dv3d PNC;
		clArray1Di NPN;
		clArray1Dd TaDe;
		clArray2Di AN;
		clArray1Dd AA0;
#ifdef INOUTNORMAL
		clArray1Dv3d InoutNorms;
#endif // INOUTNORMAL

		clArray1Dc DO;

		clArray2Di BLc;//, BLc0;
		clArray2D <clArray1Dsti> BLn;

		void ini()
		{
			loadparameters("load" SLASH "par0.m");
			loadparameters("load" SLASH "par.m");

			N.ini();
			iniA();
			iniO();
			iniB();
			//iniP();
			iniT();
			M.ini();
			K.ini();
			AK.ini();
			Ds.ini();
			//L0.ini();

		};
		void update()
		{
			L0.update();
		};
		void updateO();

		void iniA();
		void iniAA0();
		void iniB();
		void iniBL();
		void iniT();
		void iniO();
		void iniP();

		void setOp(int n, int p);

		void removeB(int n);
		void addB(int bn, int bf, int bs);

		void testBL(int ind, int n0, int n1, int n2);
		void testInOut(int ind);
		void testflatBL(int n0, int n1, int n2);
		void countBL(void);
		void compressarrayBL(void);
		BOOL removecountBL(int nbl);
		void removesingleBL(void);
		void testdoubleBL(void);
		void removedoubleBL(void);
		void removeextraBL(void);
		void removeperiodicBL(void);
		void removenoBL(void);
		BOOL testoverlapBL(int n0, int n1, clVector3Di BL0, clVector3Di BL1);
		clVector3Dd directionBL(int n);
		void sort();

		void removeBL(int n);
		void removeBL(int n1, int n2);
		void addBL(int n0, int n1, int n2);
		class clN: public clArray2Di
		{
		public:
			clArray2Di S, Sind;
			clArray1Di NL;
			int nS, nSf;
			void ini();
			void test_links();
			BOOL test_node(int n);
			void update(){}
			void save2file()
			{
				S.savetxt("lsns");
#ifdef LS_RUPTURE
				clArray2Di::save("lsn");
#else //LS_RUPTURE
				clArray2Di::savetxt("lsn");
#endif //LS_RUPTURE
			};
			void remove_links(int n1, int n2)
			{
				if(n1 == n2 || n1 < 0 || n2 < 0)
					return;

				remove_link(n1, n2);
				remove_link(n2, n1);
			};
			void remove_all_links(int n1)
			{
				if(n1 < 0 || NL(n1) == 0)
					return;

				NL(n1) = 0;
				c.ls.nNd++;
				c.ls.Tno_c(n1) = TRUE;
				c.ls.removeB(n1);
				for(int ind = 0; ind < c.ls.nL; ind++)
					if(c.ls.N(n1,ind) >= 0)
					{
						remove_link(c.ls.N(n1,ind), n1);
						c.ls.N(n1,ind) = -1;
					}
			};
			void remove_link(int n1, int n2)
			{
				if(n1 < 0 || n2 < 0 || NL(n1) == 0)
					return;

				for(int ind = 0; ind < c.ls.nL; ind++)
					if(c.ls.N(n1,ind) == n2)
					{
						el(n1,ind) = -1;
						NL(n1)--;
						if(NL(n1) <= 0)
						{
							NL(n1) = 0;
							c.ls.Tno_c(n1) = TRUE;
							c.ls.removeB(n1);
							c.ls.nNd++;
						}
						break;
					}
			};
			//void remove_link_ind(int n1, int ind)
			//{
			//	if(NL(n1) == 0 || el(n1,ind) < 0)
			//		return;

			//	el(n1,ind) = -1;
			//	NL(n1)--;
			//	if(NL(n1) <= 0)
			//	{
			//		NL(n1) = 0;
			//		c.ls.Tno_c(n1) = TRUE;
			//	}
			//};
			void setNL()
			{
				c.ls.nNd = 0;
				for(int i = 0; i < c.ls.nN; i++)
				{
					NL(i) = count_links(i);
					if(NL(i) == 0)
						c.ls.nNd++;
				}
			};
			void iniS();
			int findS1(int n0, int n1, int ns);
			int findS2(int n0, int n1, int ns);
			int find(int n0, int n1)
			{
				for(int i = 0; i < c.ls.nL; i++)
					if(el(n0, i) == n1)
						return i;
				return -1;
			}
			BOOL find_next(int n0, int n1)
			{
				for(int i = 0; i < c.ls.nL; i++)
				{
					int n = el(n0, i);
					if(n == n1 || (n >= 0 && find(n, n1) >= 0))
						return TRUE;
				}
				return FALSE;
			}
			int count_links(int i)
			{
				int count = 0;
				for(int j = 0; j < c.ls.nL; j++)
					if(el(i,j) >= 0)
						count++;
				return count;
			}

		} N;
		class clM: public clArray1Dd
		{
		public:
			double total, totalp, totala;
			clArray1Dd O;
			void ini();
			void update(){}
			void save2file()
			{
				clArray1Dd::savetxt("lsm");
				O.savetxt("lsmo");
			};
		} M;


		class clK: public clArray2Dd
		{
		public:
			clArray2Dd E;
			void ini();
			void update(){}
			void save2file()
			{
				clArray2Dd::savetxt("lsk");
			};
			void zerosE()
			{
				E.zeros();
			};

		} K;

		class clAK: public clArray1Dd
		{
		public:
			clArray1Dd E;
			void ini();
			void update(){}
			void save2file()
			{
				clArray1Dd::savetxt("lsak");
			};
			void zerosE()
			{
				E.zeros();
			};

		} AK;

		class clDs: public clArray2Dd
		{
		public:
			void ini();
			void update(){}
			void save2file()
			{
				clArray2Dd::savetxt("lsds");
			};
		} Ds;

		class clL0: public clArray2Dd
		{
		public:
			clArray2Dd L0ini;
			void ini();
			void update();
			void save2file()
			{
				clArray2Dd::save("lsl0");
			};
		} L0;
		void save2file()
		{
#ifdef LS_RUPTURE
			N.save2file();
			c.ls.Op.savetxt("lsop");
			c.ls.BNind.savetxt("lsbnind");
			c.ls.BN.save("lsbn");
			c.ls.BS.save("lsbs");
			c.ls.BF.save("lsbf");
			c.ls.BL.save("lsbl");
#endif //LS_RUPTURE

		};

	} ls;

	class clLB: protected clParemeters  
	{
		BOOL loadparameters0(const char * FileName)
		{
			int nFt = 0;
			getparameter(FileName, "c_lb_nX = %d", &nX);
			getparameter(FileName, "c_lb_nY = %d", &nY);
			getparameter(FileName, "c_lb_nZ = %d", &nZ);
			getparameter(FileName, "c_lb_nF = %d", &nFt);
			getparameter(FileName, "c_lb_UX0x = %lg", &UX0.x);
			getparameter(FileName, "c_lb_UX0y = %lg", &UX0.y);
			getparameter(FileName, "c_lb_UX0z = %lg", &UX0.z);
			getparameter(FileName, "c_lb_UX1x = %lg", &UX1.x);
			getparameter(FileName, "c_lb_UX1y = %lg", &UX1.y);
			getparameter(FileName, "c_lb_UX1z = %lg", &UX1.z);
			getparameter(FileName, "c_lb_UY0x = %lg", &UY0.x);
			getparameter(FileName, "c_lb_UY0y = %lg", &UY0.y);
			getparameter(FileName, "c_lb_UY0z = %lg", &UY0.z);
			getparameter(FileName, "c_lb_UY1x = %lg", &UY1.x);
			getparameter(FileName, "c_lb_UY1y = %lg", &UY1.y);
			getparameter(FileName, "c_lb_UY1z = %lg", &UY1.z);
			getparameter(FileName, "c_lb_UZ0x = %lg", &UZ0.x);
			getparameter(FileName, "c_lb_UZ0y = %lg", &UZ0.y);
			getparameter(FileName, "c_lb_UZ0z = %lg", &UZ0.z);
			getparameter(FileName, "c_lb_UZ1x = %lg", &UZ1.x);
			getparameter(FileName, "c_lb_UZ1y = %lg", &UZ1.y);
			getparameter(FileName, "c_lb_UZ1z = %lg", &UZ1.z);

			getparameter(FileName, "c_lb_CoarseStep = %d", &CoarseStep);
			getparameter(FileName, "c_lb_nXc = %d", &nXc);
			getparameter(FileName, "c_lb_nYc = %d", &nYc);
			getparameter(FileName, "c_lb_nZc = %d", &nZc);
			getparameter(FileName, "c_lb_Xc0 = %d", &Xc0);
			getparameter(FileName, "c_lb_Yc0 = %d", &Yc0);
			getparameter(FileName, "c_lb_Zc0 = %d", &Zc0);

			getparameter(FileName, "c_lb_bcWallX0 = %d", &bcWallX0);
			getparameter(FileName, "c_lb_bcWallX1 = %d", &bcWallX1);
			getparameter(FileName, "c_lb_bcWallY0 = %d", &bcWallY0);
			getparameter(FileName, "c_lb_bcWallY1 = %d", &bcWallY1);
			getparameter(FileName, "c_lb_bcWallZ0 = %d", &bcWallZ0);
			getparameter(FileName, "c_lb_bcWallZ1 = %d", &bcWallZ1);
			getparameter(FileName, "c_lb_bcFreeX0 = %d", &bcFreeX0);
			getparameter(FileName, "c_lb_bcFreeX1 = %d", &bcFreeX1);
			getparameter(FileName, "c_lb_bcFreeY0 = %d", &bcFreeY0);
			getparameter(FileName, "c_lb_bcFreeY1 = %d", &bcFreeY1);
			getparameter(FileName, "c_lb_bcFreeZ0 = %d", &bcFreeZ0);
			getparameter(FileName, "c_lb_bcFreeZ1 = %d", &bcFreeZ1);
			getparameter(FileName, "c_lb_bcSymmetryX0 = %d", &bcSymmetryX0);
			getparameter(FileName, "c_lb_bcSymmetryX1 = %d", &bcSymmetryX1);
			getparameter(FileName, "c_lb_bcSymmetryY0 = %d", &bcSymmetryY0);
			getparameter(FileName, "c_lb_bcSymmetryY1 = %d", &bcSymmetryY1);
			getparameter(FileName, "c_lb_bcSymmetryZ0 = %d", &bcSymmetryZ0);
			getparameter(FileName, "c_lb_bcSymmetryZ1 = %d", &bcSymmetryZ1);
			getparameter(FileName, "c_lb_bcCoarseX0 = %d", &bcCoarseX0);
			getparameter(FileName, "c_lb_bcCoarseX1 = %d", &bcCoarseX1);
			getparameter(FileName, "c_lb_bcCoarseY0 = %d", &bcCoarseY0);
			getparameter(FileName, "c_lb_bcCoarseY1 = %d", &bcCoarseY1);
			getparameter(FileName, "c_lb_bcCoarseZ0 = %d", &bcCoarseZ0);
			getparameter(FileName, "c_lb_bcCoarseZ1 = %d", &bcCoarseZ1);
			getparameter(FileName, "c_lb_bcPeriodicX = %d", &bcPeriodicX);
			getparameter(FileName, "c_lb_bcPeriodicY = %d", &bcPeriodicY);
			getparameter(FileName, "c_lb_bcPeriodicZ = %d", &bcPeriodicZ);

			if(nFt > 0)
				nF = nFt + 1;

			if(c.ls.nBFmax + 1 > nF)
				nF = c.ls.nBFmax + 1;

			int n = counttxt("load" SLASH "lbsize");
			if(n == NUMBER_DIMENSIONS)
			{
				clArray1Dd size;
				size.zeros(n);
				size.loadtxt("load" SLASH "lbsize");
				nX = (int)size(0);
				nY = (int)size(1);
				nZ = (int)size(2);
			}

			n = counttxt("load" SLASH "lbsizec");
			if(n == NUMBER_DIMENSIONS * 2)
			{
				clArray1Dd size;
				size.zeros(n);
				size.loadtxt("load" SLASH "lbsizec");
				nXc = (int)size(0);
				nYc = (int)size(1);
				nZc = (int)size(2);
				Xc0 = (int)size(3);
				Yc0 = (int)size(4);
				Zc0 = (int)size(5);
			}

			nn = clVector3Di(nX, nY, nZ);

			return TRUE;
		}
		BOOL loadparameters1(const char * FileName)
		{
			for(int i = 1; i < c.lb.nF; i++)
			{
				char Buf[80];
				sprintf(Buf, "c_lb_lamdab%d = %%lg", i);
				getparameter(FileName, Buf, &lamdab[i]);
				sprintf(Buf, "c_lb_tau%d = %%lg", i);
				getparameter(FileName, Buf, &tau[i]);
				sprintf(Buf, "c_lb_Ro%d = %%lg", i);
				getparameter(FileName, Buf, &Ro[i]);
				sprintf(Buf, "c_lb_Fx%d = %%lg", i);
				getparameter(FileName, Buf, &F[i].x);
				sprintf(Buf, "c_lb_Fy%d = %%lg", i);
				getparameter(FileName, Buf, &F[i].y);
				sprintf(Buf, "c_lb_Fz%d = %%lg", i);
				getparameter(FileName, Buf, &F[i].z);
				lamda[i] = -1. / tau[i];
			}
			return TRUE;
		}
	public:
		const static int nD, nL;
		const static int rev[LB_NUMBER_CONNECTIONS];
		const static int per[LB_NUMBER_CONNECTIONS];
		const static int Cx[LB_NUMBER_CONNECTIONS], Cy[LB_NUMBER_CONNECTIONS], Cz[LB_NUMBER_CONNECTIONS];
		static clVector3Dd C[LB_NUMBER_CONNECTIONS];
		const static double F0[LB_NUMBER_CONNECTIONS];
		const static double M[LB_NUMBER_CONNECTIONS][LB_NUMBER_CONNECTIONS], M1[LB_NUMBER_CONNECTIONS][LB_NUMBER_CONNECTIONS];
		const static int LF[LB_NUMBER_CONNECTIONS];
		static int nX, nY, nZ;
		static int CoarseStep;
		static int nXc, nYc, nZc;
		static int Xc0, Yc0, Zc0;
		static int nF;
		static clVector3Di nn;
		static clVector3Dd UX0, UX1, UY0, UY1, UZ0, UZ1;
		const static double RoX0, RoX1, RoY0, RoY1, RoZ0, RoZ1;
		const static double Cs, Cs2;
		static BOOL bcWallX0, bcWallY0, bcWallZ0, bcWallX1, bcWallY1, bcWallZ1;
		static BOOL bcPressX0, bcPressY0, bcPressZ0, bcPressX1, bcPressY1, bcPressZ1;
		static BOOL bcSymmetryX0, bcSymmetryY0, bcSymmetryZ0, bcSymmetryX1, bcSymmetryY1, bcSymmetryZ1;
		static BOOL bcCoarseX0, bcCoarseY0, bcCoarseZ0, bcCoarseX1, bcCoarseY1, bcCoarseZ1;
		static BOOL bcFreeX0, bcFreeY0, bcFreeZ0, bcFreeX1, bcFreeY1, bcFreeZ1;
		static BOOL bcPeriodicX, bcPeriodicY, bcPeriodicZ;

		clArray1Dd lamda;
		clArray1Dd lamdab;
		clArray1Dd tau;
		clArray1Dd Ro;

		clArray1Dv3d F, G;
		clArray3Dd FM;

		void ini()
		{
			loadparameters0("load" SLASH "par0.m");
			loadparameters0("load" SLASH "par.m");

			lamda.zeros(c.lb.nF);
			lamdab.zeros(c.lb.nF);
			tau.zeros(c.lb.nF);
			Ro.zeros(c.lb.nF);
			F.zeros(c.lb.nF);

			lamda.fill(LB_LAMDA_NUMBER);
			lamdab.fill(LB_LAMDAB_NUMBER);
			tau.fill(LB_TAU_NUMBER);
			Ro.fill(LB_RHO_NUMBER);
			F.fill(clVector3Dd(9999.,9999.,9999.));

			loadparameters1("load" SLASH "par0.m");
			loadparameters1("load" SLASH "par.m");

			lamdab.loadtxt("load" SLASH "lblamdab");
			tau.loadtxt("load" SLASH "lbtau");
			for(int i = 0; i < c.lb.nF; i++)
			{
				if(tau(i) != 0)
					lamda(i) = -1./tau(i);
				else
					lamda(i) = LB_LAMDA_NUMBER;
			}

			Ro.loadtxt("load" SLASH "lbro0");

			if(F.loadtxt("load" SLASH "lbf") == FALSE)
			{
				for(int i = 0; i < c.lb.nF; i++)
				{
					if(F(i) == clVector3Dd(9999.,9999.,9999.))
						F(i) = c.G * Ro(i);
				}
			}

			FM.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);
			if(FM.loadtxt("load" SLASH "lbfm") == FALSE)
				FM.fill(1.);

			save2file();

		};
		void save2file()
		{
			lamda.savetxt("lblamda");
			lamdab.savetxt("lblamdab");
			tau.savetxt("lbtau");
			Ro.savetxt("lbro0");
			F.savetxt("lbf");
			FM.savetxt("lbfm");
		};
		void update()
		{
		};
	} lb;

	class clTr: protected clParemeters  
	{
		BOOL loadparameters(const char * FileName)
		{
			getparameter(FileName, "c_tr_nN = %d", &nN);

			int n = counttxt("load" SLASH "trc") / NUMBER_DIMENSIONS;
			if(n > 0) nN = n;

			return TRUE;
		}
	public:
		static int nN;
		static double D0;
		static double TRUNC;
		static double K;
		static int Stop_dist_x;
		static int remain;

		void ini()
		{
			loadparameters("load" SLASH "par0.m");
			loadparameters("load" SLASH "par.m");
		};
		void save2file()
		{
		};
		void update()
		{
		};
	} tr;

};


#endif // !defined(AFX_CLCONSTANTS_H__INCLUDED_)
