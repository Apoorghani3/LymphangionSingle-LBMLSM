// clOutput.cpp: implementation of the clOutput class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clOutput.h"
#include "clConstants.h"
//#include "clGrid.h"
#include "clMethod.h"
#include "clVariables.h"

void clOutput::ini()
{
	Timeout.zeros(OUTPUT_MAX_LINES);
	LBout.zeros(OUTPUT_MAX_LINES, c.lb.nF, OUTPUT_LB_FILEDS);
	LSout.zeros(OUTPUT_MAX_LINES, c.ls.nO, OUTPUT_LS_FILEDS);
	curLine = 0;

	outDir = OUTPUT_DIR;
	if(c.t.flOutputDump == TRUE)
	{
#if defined(_WIN32)
		_mkdir(outDir);
#else 
		mkdir(outDir, S_IREAD | S_IWRITE | S_IEXEC | S_IRGRP | S_IXGRP | S_IRUSR | S_IXUSR);
#endif //_WIN32
	}

	if(c.t.flContinue == FALSE)
		CleanFile(CURR_OUTPUT_FILE);

	avrCounter = 0;
	avrPeriod = 0;

};

void clOutput::RenameOutputFiles()
{
//	hCurrOut = CopyCurrState(hCurrOut, CURR_OUTPUT_FILE, TMP_OUTPUT_FILE);
	CopyCurrState(CURR_OUTPUT_FILE);

	if(c.t.Time < c.t.DumpStartTime || c.t.flOutputDump == FALSE)
		return;

	countRename = (int)((c.t.Time - c.t.DumpStartTime) / c.t.DumpTimeStepSave);

	MakeDir(countRename);

#ifdef LB_SOLVER
	RenameFile("lbro");
	RenameFile("lbm");
	RenameFile("lbj");
	RenameFile("lbn");
//	RenameFile("lbmsf");
	RenameFile("lbj_cent");	
	RenameFile("lbro_cent");
	RenameFile("lbj_rate");	
	RenameFile("lbj_Qrate");	
#endif //LB_SOLVER

#ifdef LS_SOLVER
	RenameFile("lsc");
	RenameFile("lsv");
	//RenameFile("lsblfb");
	RenameFile("lsfb");

#ifdef LS_RUPTURE
	RenameFile("lsn");
	RenameFile("lsbn");
	RenameFile("lsbl");
//	RenameFile("lsbs");
	RenameFile("lsop");
#endif //LS_RUPTURE
#endif //LS_SOLVER

#ifdef TR_SOLVER
	//RenameFile("trc");
	//RenameFile("trv");
	//RenameFile("tru");
	//RenameFile("tra");
	RenameFile("trbldep");
	RenameFile("trConc");
#endif //TR_SOLVER

//	RenameFile("lstpba");

#ifdef RENAME_BIN_FILES
#ifdef LB_SOLVER
	RenameFileExt("lbf.bin");
#endif //LB_SOLVER

#ifdef LS_SOLVER
	RenameFileExt("lsc.bin");
	RenameFileExt("lsv.bin");
	RenameFileExt("lsf.bin");
#endif //LS_SOLVER

#ifdef TR_SOLVER
	RenameFileExt("trc.bin");
#endif //LS_SOLVER
#endif //RENAME_BIN_FILES

	RenameFileExt("par.m");
}

void clOutput::field()
{
	if(c.t.Time != 0.)
		RenameOutputFiles();

	parameters();

#ifdef LB_SOLVER
	vlb.save2file();
#endif //LB_SOLVER

#ifdef LS_SOLVER
	vls.save2file();
#endif //LS_SOLVER

#ifdef TR_SOLVER
	vtr.save2file();
#endif //TR_SOLVER

	saveBin();
#ifdef SAVE_BIN_IN_LOAD
	CopyFile("par.m", "load" SLASH "par.m");
#endif //SAVE_BIN_IN_LOAD

#ifdef _DEBUG
#endif //_DEBUG
};


void clOutput::parameters()
{
	c.save2file();
}

void clOutput::start()
{
	average();
	bin();
	parameters();
	if(c.t.flContinue == FALSE)
		current();
	field();
}


void clOutput::current()
{
	if(c.t.CurrentSaveStartTime > c.t.Time)
	{
		vls.F.zerosW();
		return;
	}

	if(curLine >= OUTPUT_MAX_LINES)
		CopyCurrState(CURR_OUTPUT_FILE);

	vls.avr.update();
	vlb.updateoutput();
#ifdef BINARY_FLUID
	vlb.N.update();
#endif // BINARY_FLUID
	//vlb.NM.updateoutput();
#ifdef TR_SOLVER
	vtr.updateoutput();
#endif // TR_SOLVER

	//Modified
	clArray1Di Rad_num; Rad_num.zeros(c.lb.nX);
	clArray1Dd Rad_t; Rad_t.zeros(c.lb.nX);
	
	double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;

	//Adds radius value to corresponding x integer coordinate twice (floor + ceiling). Exact values count twice
	for(int i = 0; i < c.ls.nN; i++)
	{
		if(c.ls.O(i)== (NUM_VESSELS-1))
		{
			int flx = (int)(floor(vls.C(i).x)), clx = (int)(ceil(vls.C(i).x));
			if( (flx >= 0) & (flx <= (c.lb.nX-1)) )
			{
				Rad_t(flx) += sqrt( pow((vls.C(i).y-midy),2.0) + pow((vls.C(i).z-midz),2.0) );
				Rad_num(flx)++;				
			}
			if( (clx >= 0) & (clx <= (c.lb.nX-1)) )
			{
				Rad_t(clx) += sqrt( pow((vls.C(i).y-midy),2.0) + pow((vls.C(i).z-midz),2.0) );
				Rad_num(clx)++;				
			}			
		}			
	}
	
	for(int i = 0; i < c.lb.nX; i++)	
		Rad_t(i) /= Rad_num(i);

	clArray3Di VesselIn; VesselIn.zeros(c.lb.nX, c.lb.nY, c.lb.nZ);	
	int lx = SAMPX, hx = c.lb.nX-1-SAMPX;
	int centl=(int)(floor((c.lb.nX-1.0)/2.0));
	int centh=(int)(ceil((c.lb.nX-1.0)/2.0));
	double avgl=0.0, avgm=0.0, avgh=0.0, Qavgl=0.0, Qavgm=0.0, Qavgh=0.0;

	clArray1Dd Qx_t; Qx_t.zeros(c.lb.nX);
	clArray1Dd Jx_t; Jx_t.zeros(c.lb.nX);

	for(int i = 0; i < c.lb.nX; i++)
	{
		for(int j = 0; j < c.lb.nY; j++)
		{
			for(int z =0; z < c.lb.nZ; z++)
			{
				double lbrad = sqrt(pow((j-midy),2.0)+pow((z-midz),2.0));
				VesselIn(i,j,z)=(Rad_t(i) > lbrad);
				int vcrit = VesselIn(i,j,z);

				Jx_t(i) += vlb.J(i,j,z).x * (double)(vcrit);	
				Qx_t(i) += vlb.J(i,j,z).x/vlb.Ro(i,j,z) * (double)(vcrit);

				VesselIn(lx,j,z)=(Rad_t(lx) > lbrad);
				VesselIn(hx,j,z)=(Rad_t(hx) > lbrad);
				VesselIn(centl,j,z)=(Rad_t(centl) > lbrad);
				VesselIn(centh,j,z)=(Rad_t(centh) > lbrad);
			}
		}
	}

	avgl=Jx_t(lx);
	avgm=(Jx_t(centl) + Jx_t(centh))/2.0;
	avgh=Jx_t(hx);

	Qavgl=Qx_t(lx);
	Qavgm=(Qx_t(centl) + Qx_t(centh))/2.0;
	Qavgh=Qx_t(hx);

	//Jx_t.savetxt("lbj_rate");
	//Qx_t.savetxt("lbj_Qrate");

	printf("\n avgl = %f \n",avgl);

	//Modified


	Timeout(curLine) = c.t.Time + c.t.StartTime;
	int ind;
	for(int i = 1; i < c.lb.nF; i++)
	{
		ind = 0;
		LBout(curLine, i, ind++) = vlb.Ro.mass(i); //2
		//LBout(curLine, i, ind++) = vlb.Ro.avr(i); 
		//LBout(curLine, i, ind++) = vlb.J.avr(i).x;
		//LBout(curLine, i, ind++) = vlb.J.avr(i).y;
		//LBout(curLine, i, ind++) = vlb.J.avr(i).z;
		//LBout(curLine, i, ind++) = vlb.J.tot(i).x;
		//LBout(curLine, i, ind++) = vlb.J.tot(i).y;
		//LBout(curLine, i, ind++) = vlb.J.tot(i).z;
		//LBout(curLine, i, ind++) = vtr.A.tot_cilia;
		LBout(curLine, i, ind++) = vlb.NM.S; //3
		//LBout(curLine, i, ind++) = vtr.C.OutputVal.x;
		//LBout(curLine, i, ind++) = vtr.C.OutputVal.y;
		//LBout(curLine, i, ind++) = vtr.C.OutputVal.z;
		//LBout(curLine, i, ind++) = vlb.J.invel;
		//Modified
		LBout(curLine, i, ind++) = Qavgl; //Average flow rate at lower x 4
		LBout(curLine, i, ind++) = Qavgm; //Average flow rate at mid x   5
		LBout(curLine, i, ind++) = Qavgh; //Average flow rate at high x  6
		LBout(curLine, i, ind++) = avgl; //Average mass flow rate at lower x  7
		LBout(curLine, i, ind++) = avgm; //Average mass flow rate at mid x   8
		LBout(curLine, i, ind++) = avgh; //Average mass flow rate at high x   9
		LBout(curLine, i, ind++) = vlb.Ro.max; //Maximum density deviation (for stability check)   10
		LBout(curLine, i, ind++) = vlb.Ro.avgabs; //Average absolute value density deviation (for stability check)  11
		//Modified
	}

	for(int i = 0; i < c.ls.nO; i++)
	{
		ind = 0;
		LSout(curLine, i, ind++) = vls.avr.Fb(i).x;//32-77
		LSout(curLine, i, ind++) = vls.avr.Fb(i).y;
		LSout(curLine, i, ind++) = vls.avr.Fb(i).z; 
		LSout(curLine, i, ind++) = vls.avr.Fc(i).x;//35-80
		LSout(curLine, i, ind++) = vls.avr.Fc(i).y;
		LSout(curLine, i, ind++) = vls.avr.Fc(i).z; 
		LSout(curLine, i, ind++) = vls.avr.Fe(i).x;//38-83
		LSout(curLine, i, ind++) = vls.avr.Fe(i).y;
		LSout(curLine, i, ind++) = vls.avr.Fe(i).z; 
		LSout(curLine, i, ind++) = vls.avr.CG(i).x;//41-86
		LSout(curLine, i, ind++) = vls.avr.CG(i).y;
		LSout(curLine, i, ind++) = vls.avr.CG(i).z; 
		LSout(curLine, i, ind++) = vls.avr.Xmax(i).x;//44-89
		LSout(curLine, i, ind++) = vls.avr.Xmax(i).y;//90
		LSout(curLine, i, ind++) = vls.avr.Xmax(i).z; //91
		LSout(curLine, i, ind++) = vls.avr.Xmin(i).x;//47
		LSout(curLine, i, ind++) = vls.avr.Xmin(i).y;
		LSout(curLine, i, ind++) = vls.avr.Xmin(i).z; 
		LSout(curLine, i, ind++) = vls.avr.V(i).x;
		LSout(curLine, i, ind++) = vls.avr.V(i).y;
		LSout(curLine, i, ind++) = vls.avr.V(i).z;
		LSout(curLine, i, ind++) = vls.avr.Vmax(i); 
		LSout(curLine, i, ind++) = vls.avr.P(i).x;
		LSout(curLine, i, ind++) = vls.avr.P(i).y;
		LSout(curLine, i, ind++) = vls.avr.P(i).z;
		LSout(curLine, i, ind++) = vls.avr.O(i).x;
		LSout(curLine, i, ind++) = vls.avr.O(i).y;
		LSout(curLine, i, ind++) = vls.avr.O(i).z;
		LSout(curLine, i, ind++) = vls.avr.L(i).x;
		LSout(curLine, i, ind++) = vls.avr.L(i).y;
		LSout(curLine, i, ind++) = vls.avr.L(i).z;
		LSout(curLine, i, ind++) = vls.avr.Ravr(i);
		LSout(curLine, i, ind++) = vls.avr.Rmax(i);
		LSout(curLine, i, ind++) = vls.avr.Rmin(i);
		LSout(curLine, i, ind++) = vls.avr.Ek(i);
		LSout(curLine, i, ind++) = vls.avr.Ec(i);
		LSout(curLine, i, ind++) = vls.avr.Ev(i);
		LSout(curLine, i, ind++) = vls.avr.Wb(i);//69
		LSout(curLine, i, ind++) = vls.avr.Wc(i);//70
		LSout(curLine, i, ind++) = vls.avr.We(i);//71
		LSout(curLine, i, ind++) = vls.avr.Smax(i);
		LSout(curLine, i, ind++) = vls.avr.Smin(i);

		LSout(curLine, i, ind++) = vls.R.T(i).x;
		LSout(curLine, i, ind++) = vls.R.T(i).y;
		LSout(curLine, i, ind++) = vls.R.T(i).z;
	}

	curLine++;

	//fprintf(hCurrOut, "%f", c.t.Time + c.t.StartTime);
	//for(int i = 1; i < c.lb.nF; i++)
	//	fprintf(hCurrOut, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g", 
	//		vlb.Ro.mass(i), vlb.Ro.avr(i), 
	//		vlb.J.avr(i).x, vlb.J.avr(i).y, vlb.J.avr(i).z, 
	//		vlb.J.tot(i).x, vlb.J.tot(i).y, vlb.J.tot(i).z);

	//for(int i = 0; i < c.ls.nO; i++)
	//	fprintf(hCurrOut, "\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g",
	//		vls.avr.Fb(i).x, vls.avr.Fb(i).y, vls.avr.Fb(i).z, 
	//		vls.avr.Fc(i).x, vls.avr.Fc(i).y, vls.avr.Fc(i).z, 
	//		vls.avr.Fe(i).x, vls.avr.Fe(i).y, vls.avr.Fe(i).z, 
	//		vls.avr.CG(i).x, vls.avr.CG(i).y, vls.avr.CG(i).z, 
	//		vls.avr.Xmax(i).x, vls.avr.Xmax(i).y, vls.avr.Xmax(i).z, 
	//		vls.avr.Xmin(i).x, vls.avr.Xmin(i).y, vls.avr.Xmin(i).z, 
	//		vls.avr.V(i).x, vls.avr.V(i).y, vls.avr.V(i).z, vls.avr.Vmax(i), 
	//		vls.avr.P(i).x, vls.avr.P(i).y, vls.avr.P(i).z, 
	//		vls.avr.O(i).x, vls.avr.O(i).y, vls.avr.O(i).z, 
	//		vls.avr.L(i).x, vls.avr.L(i).y, vls.avr.L(i).z, 
	//		vls.avr.Ravr(i), vls.avr.Rmax(i), vls.avr.Rmin(i),
	//		vls.avr.Ek(i), vls.avr.Ec(i), vls.avr.Ev(i),
	//		vls.avr.Wb(i), vls.avr.Wc(i), vls.avr.We(i),
	//		vls.avr.Smax(i), vls.avr.Smin(i));

	//fprintf(hCurrOut, "\n");
}

void clOutput::average()
{
}


void clOutput::bin()
{
}

void clOutput::addAverage()
{
}


void clOutput::saveBin()
{
	//////////saveBin(v.P, "p.bin");
	//////////saveBin(v.U, "u.bin");
	//////////saveBin(v.V, "v.bin");
	//////////saveBin(v.T, "t.bin");
	//////////saveBin(v.VoF, "vof.bin");
}
