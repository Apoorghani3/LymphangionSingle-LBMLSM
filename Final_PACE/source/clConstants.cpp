// clConstants.cpp: implementation of the clConstants class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "clConstants.h"
#include "clOutput.h"
#include "clProblem.h"
#include "clVariables.h"

const double clConstants::R = 287.;
const double clConstants::Pi = PI_NUMBER;
const double clConstants::Kb = 1.38066e-23;
double clConstants::Re1 = 1.;
double clConstants::Pe1 = 1.;
//clVector3Dd clConstants::G = clVector3Dd(POISEUILLE_FLAG * 2. * GAMMA_NUMBER * MU_NUMBER / CHANNEL_WIDTH, 0., 0.);
//clVector3Dd clConstants::G = clVector3Dd(POISEUILLE_FLAG * GX_NUMBER, 0., POISEUILLE_FLAG * GZ_NUMBER);
clVector3Dd clConstants::G = clVector3Dd(RE_NUMBER, 0., POISEUILLE_FLAG * GZ_NUMBER);
//clVector3Dd clConstants::G = clVector3Dd(1.e-3, 0., 0.);
//clVector3Dd clConstants::G = clVector3Dd(0., 0.1e-4, 0.);
//clVector3Dd clConstants::G = clVector3Dd(0., 0., 1.e-6);
double clConstants::sL = 1.;
double clConstants::sT = 1.;
double clConstants::sM = 1.;
double clConstants::sN = 1.;
double clConstants::L0 = 2.e-7;
double clConstants::Ro0 = 1000.;
double clConstants::nu0 = 1.e-6 * 6./1.;
double clConstants::T0 = 298.;
int clConstants::nX = 1;
int clConstants::nY = 1;
int clConstants::nZ = 1;
double clConstants::eps = 1.e-5;
clVector3Di clConstants::nn = clVector3Di(clConstants::nX, clConstants::nY, clConstants::nZ);


//////////////////////////////////////////////////////////////////////////////////////
BOOL clConstants::clTime::flContinue = 0;
BOOL clConstants::clTime::flOutputDump = FLAG_DUMP;
double clConstants::clTime::dTlb = LB_STEPS;
double clConstants::clTime::StopTime = STOP_TIME;

int clConstants::clTime::lsSteps = LS_STEPS;
int clConstants::clTime::trSteps = TR_STEPS;
double clConstants::clTime::dTls = clConstants::clTime::dTlb / (double)(clConstants::clTime::lsSteps);
double clConstants::clTime::dTtr = clConstants::clTime::dTlb / (double)(clConstants::clTime::trSteps);
double clConstants::clTime::Time = 0.;
double clConstants::clTime::TimeN = 0.;
double clConstants::clTime::StartTime = 0.;
double clConstants::clTime::DumpStartTime = DUMP_START_TIME;
double clConstants::clTime::DumpTimeStep = DUMP_STEP;
double clConstants::clTime::DumpTimeStepSave = DUMP_STEP_SAVE;
#ifdef _DEBUG
int clConstants::clTime::CurrentSaveTimeStep = 1;
int clConstants::clTime::CurrentSaveStartTime = 0;
int clConstants::clTime::DisplaySignal = 1;
#else //_DEBUG
int clConstants::clTime::CurrentSaveTimeStep = CURRENT_SAVE_STEP;
int clConstants::clTime::CurrentSaveStartTime = CURRENT_SAVE_START_TIME;
#if defined(_WIN32)
int clConstants::clTime::DisplaySignal = ((int)(clConstants::clTime::StopTime / 10000.) < 1 ? 1 : (int)(clConstants::clTime::StopTime / 10000.));
#else //_WIN32
int clConstants::clTime::DisplaySignal = ((int)(clConstants::clTime::StopTime / 100.) < 1 ? 1 : (int)(clConstants::clTime::StopTime / 100.));
#endif //_WIN32
#endif //_DEBUG

double clConstants::clTime::oscRatio = 1;
double clConstants::clTime::oscAmpl1 = OSC_AMPL;
//double clConstants::clTime::oscAmplC2 = OSC_AMPL_2;
double clConstants::clTime::oscPeriod1 = OSC_PERIOD;
double clConstants::clTime::oscAmpl2 = OSC_AMPL_2;
double clConstants::clTime::oscPeriod2 = OSC_PERIOD_2;

//////////////////////////////////////////////////////////////////
int clConstants::clLS::nN = 0;
int clConstants::clLS::nL = LS_NUMBER_CONNECTIONS;
int clConstants::clLS::nBN = 0;
int clConstants::clLS::nBL = 0;
int clConstants::clLS::nBM = 0;
int clConstants::clLS::nTa = 0;
//int clConstants::clLS::nTp = 0;
int clConstants::clLS::nTr = 0;
int clConstants::clLS::nTs = 0;
int clConstants::clLS::nO = 0;
int clConstants::clLS::nOp = 0;
int clConstants::clLS::nPN = 0;
int clConstants::clLS::nNd = 0;
int clConstants::clLS::nAN = 0;
int clConstants::clLS::nBFmax = 0;
int clConstants::clLS::nPupd = 2;

double clConstants::clLS::Ds0 = LS_DISIPATION_0 * LS_DISSIPATION_FLAG;
double clConstants::clLS::Ds1 = LS_DISIPATION_1 * LS_DISSIPATION_FLAG;
double clConstants::clLS::Ds2 = LS_DISIPATION_2 * LS_DISSIPATION_FLAG;
//double clConstants::clLS::Lb = 0.3e-6;

double clConstants::clLS::FwDe0 = LS_DE_0;
double clConstants::clLS::FwDe1 = LS_DE_1;
double clConstants::clLS::FwDe2 = LS_DE_2;
//double clConstants::clLS::FwDeR = clConstants::clLS::FwDe * 2.;
double clConstants::clLS::FwDeR = LS_DE_R_NUMBER;
double clConstants::clLS::FwRe = LS_R_EQ_NUMBER;
double clConstants::clLS::FwKe = LS_K_NUMBER;

double clConstants::clLS::FwDeC = LS_DE_C_NUMBER;
#ifdef LS_RUPTURE
double clConstants::clLS::FwReC = LS_R_EQ_C_NUMBER;
#else //LS_RUPTURE
double clConstants::clLS::FwReC = LS_R_EQ_C_NUMBER * 0.5;
#endif //LS_RUPTURE

double clConstants::clLS::Symax0 = LS_SYMAX_0;
double clConstants::clLS::Symin0 = LS_SYMIN_0;
double clConstants::clLS::Symax1 = LS_SYMAX_1;
double clConstants::clLS::Symin1 = LS_SYMIN_1;
double clConstants::clLS::Symax2 = LS_SYMAX_2;
double clConstants::clLS::Symin2 = LS_SYMIN_2;

// Valve extent (computed during initialization)
double clConstants::clLS::valveMinX = 0.;
double clConstants::clLS::valveMaxX = 0.;
double clConstants::clLS::waveLength = 0.;

double clConstants::clLS::FwDeX0 = 0.;
double clConstants::clLS::FwDeRX0 = clConstants::clLS::FwDeR;
double clConstants::clLS::FwDeX1 = 0.;
double clConstants::clLS::FwDeRX1 = clConstants::clLS::FwDeR;
double clConstants::clLS::FwDeY0 = 0.;
double clConstants::clLS::FwDeRY0 = clConstants::clLS::FwDeR;
double clConstants::clLS::FwDeY1 = 0.;
double clConstants::clLS::FwDeRY1 = clConstants::clLS::FwDeR;
double clConstants::clLS::FwDeZ0 = 0.;
double clConstants::clLS::FwDeRZ0 = clConstants::clLS::FwDeR;
double clConstants::clLS::FwDeZ1 = 0.;
double clConstants::clLS::FwDeRZ1 = clConstants::clLS::FwDeR;

clVector3Dd clConstants::clLS::Cini0 = clVector3Dd(LS_C_INI_X_0, LS_C_INI_Y_0, LS_C_INI_Z_0);
clVector3Dd clConstants::clLS::Cini1 = clVector3Dd(LS_C_INI_X_1, LS_C_INI_Y_1, LS_C_INI_Z_1);
clVector3Dd clConstants::clLS::Cini2 = clVector3Dd(LS_C_INI_X_2, LS_C_INI_Y_2, LS_C_INI_Z_2);

//double clConstants::clLS::Kini0 = LS_STIFFNESS_0 / 1.870;
double clConstants::clLS::Kini0 = LS_STIFFNESS_0;
double clConstants::clLS::Kini1 = LS_STIFFNESS_1;
double clConstants::clLS::Kini2 = LS_STIFFNESS_2;

double clConstants::clLS::AKini0 = LS_ANGLE_STIFFNESS_0;
double clConstants::clLS::AKini1 = LS_ANGLE_STIFFNESS_1;
double clConstants::clLS::AKini2 = LS_ANGLE_STIFFNESS_2;

double clConstants::clLS::Mini0 = LS_MASS_0;
double clConstants::clLS::Mini1 = LS_MASS_1;
double clConstants::clLS::Mini2 = LS_MASS_2;

double clConstants::clLS::L0ini0 = LS_EQ_LENGTH_0;
double clConstants::clLS::L0ini1 = LS_EQ_LENGTH_1;
double clConstants::clLS::L0ini2 = LS_EQ_LENGTH_2;

///////////////////////////////////////////////////////////////////
int clConstants::clLB::nX = 9;
int clConstants::clLB::nY = 9;
int clConstants::clLB::nZ = 9;
int clConstants::clLB::nXc = LB_COARSE_NX;
int clConstants::clLB::nYc = LB_COARSE_NY;
int clConstants::clLB::nZc = LB_COARSE_NZ;
int clConstants::clLB::Xc0 = LB_COARSE_X0;
int clConstants::clLB::Yc0 = LB_COARSE_Y0;
int clConstants::clLB::Zc0 = LB_COARSE_Z0;
int clConstants::clLB::CoarseStep = LB_COARSE_STEP;
clVector3Di clConstants::clLB::nn = clVector3Di(clConstants::clLB::nX, clConstants::clLB::nY, clConstants::clLB::nZ);
const int clConstants::clLB::nD = NUMBER_DIMENSIONS;
const int clConstants::clLB::nL = LB_NUMBER_CONNECTIONS;
int clConstants::clLB::nF = LB_NUMBER_COMPONENTS + 1;

////////////////////////////////////////////////////////// 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
const int clConstants::clLB::rev[LB_NUMBER_CONNECTIONS] = {0,  2,  1,  4,  3,  6,  5,  8,  7, 10,  9, 12, 11, 14, 13, 16, 15, 18, 17};
const int clConstants::clLB::per[LB_NUMBER_CONNECTIONS] = {0,  3,  4,  5,  6,  1,  2,  9, 10,  7,  8, 13, 14, 11, 12, 17, 18, 15, 16};

const int clConstants::clLB::Cx[LB_NUMBER_CONNECTIONS] =  {0,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
const int clConstants::clLB::Cy[LB_NUMBER_CONNECTIONS] =  {0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1};
const int clConstants::clLB::Cz[LB_NUMBER_CONNECTIONS] =  {0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1};
//const int clConstants::clLB::Cx[LB_NUMBER_CONNECTIONS] =  {0,  1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  1, -1, -1,  1,  1, -1,  1, -1};
//const int clConstants::clLB::Cy[LB_NUMBER_CONNECTIONS] =  {0,  0,  0,  1, -1,  0,  0,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1, -1,  1};
//const int clConstants::clLB::Cz[LB_NUMBER_CONNECTIONS] =  {0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1,  1, -1,  1, -1,  0,  0,  0,  0};
clVector3Dd clConstants::clLB::C[LB_NUMBER_CONNECTIONS] = { clVector3Dd(clConstants::clLB::Cx[ 0], clConstants::clLB::Cy[ 0], clConstants::clLB::Cz[ 0]),
															clVector3Dd(clConstants::clLB::Cx[ 1], clConstants::clLB::Cy[ 1], clConstants::clLB::Cz[ 1]),
															clVector3Dd(clConstants::clLB::Cx[ 2], clConstants::clLB::Cy[ 2], clConstants::clLB::Cz[ 2]),
															clVector3Dd(clConstants::clLB::Cx[ 3], clConstants::clLB::Cy[ 3], clConstants::clLB::Cz[ 3]),
															clVector3Dd(clConstants::clLB::Cx[ 4], clConstants::clLB::Cy[ 4], clConstants::clLB::Cz[ 4]),
															clVector3Dd(clConstants::clLB::Cx[ 5], clConstants::clLB::Cy[ 5], clConstants::clLB::Cz[ 5]),
															clVector3Dd(clConstants::clLB::Cx[ 6], clConstants::clLB::Cy[ 6], clConstants::clLB::Cz[ 6]),
															clVector3Dd(clConstants::clLB::Cx[ 7], clConstants::clLB::Cy[ 7], clConstants::clLB::Cz[ 7]),
															clVector3Dd(clConstants::clLB::Cx[ 8], clConstants::clLB::Cy[ 8], clConstants::clLB::Cz[ 8]),
															clVector3Dd(clConstants::clLB::Cx[ 9], clConstants::clLB::Cy[ 9], clConstants::clLB::Cz[ 9]),
															clVector3Dd(clConstants::clLB::Cx[10], clConstants::clLB::Cy[10], clConstants::clLB::Cz[10]),
															clVector3Dd(clConstants::clLB::Cx[11], clConstants::clLB::Cy[11], clConstants::clLB::Cz[11]),
															clVector3Dd(clConstants::clLB::Cx[12], clConstants::clLB::Cy[12], clConstants::clLB::Cz[12]),
															clVector3Dd(clConstants::clLB::Cx[13], clConstants::clLB::Cy[13], clConstants::clLB::Cz[13]),
															clVector3Dd(clConstants::clLB::Cx[14], clConstants::clLB::Cy[14], clConstants::clLB::Cz[14]),
															clVector3Dd(clConstants::clLB::Cx[15], clConstants::clLB::Cy[15], clConstants::clLB::Cz[15]),
															clVector3Dd(clConstants::clLB::Cx[16], clConstants::clLB::Cy[16], clConstants::clLB::Cz[16]),
															clVector3Dd(clConstants::clLB::Cx[17], clConstants::clLB::Cy[17], clConstants::clLB::Cz[17]),
															clVector3Dd(clConstants::clLB::Cx[18], clConstants::clLB::Cy[18], clConstants::clLB::Cz[18])};

const double clConstants::clLB::F0[LB_NUMBER_CONNECTIONS] ={1./3., 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 
															1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 
															1./36., 1./36., 1./36., 1./36., 1./36.};


const int clConstants::clLB::LF[LB_NUMBER_CONNECTIONS] = {	LB_NO_BOUNDARY_LINK, 
															LB_BOUNDARY_LINK_1, 
															LB_BOUNDARY_LINK_2,
															LB_BOUNDARY_LINK_3,
															LB_BOUNDARY_LINK_4,
															LB_BOUNDARY_LINK_5,
															LB_BOUNDARY_LINK_6,
															LB_BOUNDARY_LINK_7,
															LB_BOUNDARY_LINK_8,
															LB_BOUNDARY_LINK_9, 
															LB_BOUNDARY_LINK_10,
															LB_BOUNDARY_LINK_11,
															LB_BOUNDARY_LINK_12,
															LB_BOUNDARY_LINK_13,
															LB_BOUNDARY_LINK_14,
															LB_BOUNDARY_LINK_15,
															LB_BOUNDARY_LINK_16,
															LB_BOUNDARY_LINK_17,
															LB_BOUNDARY_LINK_18};

const double clConstants::clLB::M[LB_NUMBER_CONNECTIONS][LB_NUMBER_CONNECTIONS] = {
															{1,	1,	1,	1,	1,	1,	1,	1,	1,  1,  1,	1,	1,	1,  1,  1,  1,  1,  1},
															{-1, 0, 0,  0,  0,  0,  0,  1,	1,  1,	1,	1,	1,  1,  1,  1,  1,  1,	1}, 
															{1, -2,-2, -2, -2, -2, -2,	1,  1,  1,  1,	1,	1,  1,  1,  1,  1,	1,	1} ,
															{0, 1, -1,	0,	0,	0,	0,	1, -1,	1, -1,	1, -1,  1, -1,  0,  0,	0,	0}, 
															{0,-2,	2,	0,	0,	0,	0,	1, -1,	1, -1,  1, -1,	1, -1,  0,  0,	0,	0}, 
															{0,	0,	0,  1, -1,	0,  0,	1, -1, -1,	1,	0,  0,	0,	0,	1, -1,	1, -1}, 
															{0, 0,	0, -2,	2,	0,	0,	1, -1, -1,	1,  0,  0,	0,	0,  1, -1,	1, -1}, 
															{0, 0,	0,	0,	0,	1, -1,	0,	0,	0,	0,	1, -1, -1,	1,	1, -1, -1,	1}, 
															{0, 0,  0,	0,  0, -2,	2,	0,	0,	0,	0,  1, -1, -1,  1,  1, -1, -1,	1}, 
															{0, 2,  2, -1, -1, -1, -1,	1,	1,	1,	1,  1,  1,	1,	1, -2, -2, -2, -2},
															{0,-2, -2,	1,  1,	1,  1,  1,  1,	1,  1,  1,  1,	1,	1, -2, -2, -2, -2},
															{0,	0,	0,	1,	1, -1, -1,	1,	1,  1,  1, -1, -1, -1, -1,  0,  0,  0,  0},
															{0,	0,  0, -1, -1,	1,	1,	1,	1,	1,  1, -1, -1, -1, -1,  0,	0,  0,  0},
															{0,	0,  0,  0,  0,  0,	0,	1,	1, -1, -1,	0,	0,	0,	0,	0,	0,	0,	0},
															{0,	0,	0,  0,  0,  0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	1, -1, -1},
															{0,	0,	0,  0,  0,  0,  0,  0,  0,  0,	0,  1,  1, -1, -1,  0,	0,  0,  0},
															{0,	0,	0,  0,  0,  0,  0,  1, -1,  1, -1, -1,  1, -1,  1,  0,	0,  0,  0},
															{0,	0,  0,  0,  0,  0,  0, -1,	1,  1, -1,  0,	0,  0,	0,  1, -1,	1, -1},
															{0, 0,  0,  0,  0,	0,  0,  0,  0,  0,  0,  1, -1, -1,	1, -1,  1,	1, -1}};

const double clConstants::clLB::M1[LB_NUMBER_CONNECTIONS][LB_NUMBER_CONNECTIONS] = {
								{1./3., -1./2., 1./6., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., 1./6., -1./6., 0, 0, 0, 0, 1./12., -1./12., 0, 0, 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., -1./6., 1./6., 0, 0, 0, 0, 1./12., -1./12., 0, 0, 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., 0, 0, 1./6., -1./6., 0, 0, -1./24., 1./24., 1./8., -1./8., 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., 0, 0, -1./6., 1./6., 0, 0, -1./24., 1./24., 1./8., -1./8., 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., 0, 0, 0, 0, 1./6., -1./6., -1./24., 1./24., -1./8., 1./8., 0, 0, 0, 0, 0, 0}, 
								{1./18., 0, -1./18., 0, 0, 0, 0, -1./6., 1./6., -1./24., 1./24., -1./8., 1./8., 0, 0, 0, 0, 0, 0}, 
								{1./36., 1./24., 1./72., 1./12., 1./24., 1./12., 1./24., 0, 0, 1./48., 1./48., 1./16., 1./16., 1./4., 0, 0, 1./8., -1./8., 0}, 
								{1./36., 1./24., 1./72., -1./12., -1./24., -1./12., -1./24., 0, 0, 1./48., 1./48., 1./16., 1./16., 1./4., 0, 0, -1./8., 1./8., 0}, 
								{1./36., 1./24., 1./72., 1./12., 1./24., -1./12., -1./24., 0, 0, 1./48., 1./48., 1./16., 1./16., -1./4., 0, 0, 1./8., 1./8., 0}, 
								{1./36., 1./24., 1./72., -1./12., -1./24., 1./12., 1./24., 0, 0, 1./48., 1./48., 1./16., 1./16., -1./4., 0, 0, -1./8., -1./8., 0}, 
								{1./36., 1./24., 1./72., 1./12., 1./24., 0, 0, 1./12., 1./24., 1./48., 1./48., -1./16., -1./16., 0, 0, 1./4., -1./8., 0, 1./8.}, 
								{1./36., 1./24., 1./72., -1./12., -1./24., 0, 0, -1./12., -1./24., 1./48., 1./48., -1./16., -1./16., 0, 0, 1./4., 1./8., 0, -1./8.}, 
								{1./36., 1./24., 1./72., 1./12., 1./24., 0, 0, -1./12., -1./24., 1./48., 1./48., -1./16., -1./16., 0, 0, -1./4., -1./8., 0, -1./8.}, 
								{1./36., 1./24., 1./72., -1./12., -1./24., 0, 0, 1./12., 1./24., 1./48., 1./48., -1./16., -1./16., 0, 0, -1./4., 1./8., 0, 1./8.}, 
								{1./36., 1./24., 1./72., 0, 0, 1./12., 1./24., 1./12., 1./24., -1./24., -1./24., 0, 0, 0, 1./4., 0, 0, 1./8., -1./8.}, 
								{1./36., 1./24., 1./72., 0, 0, -1./12., -1./24., -1./12., -1./24., -1./24., -1./24., 0, 0, 0, 1./4., 0, 0, -1./8., 1./8.}, 
								{1./36., 1./24., 1./72., 0, 0, 1./12., 1./24., -1./12., -1./24., -1./24., -1./24., 0, 0, 0, -1./4., 0, 0, 1./8., 1./8.}, 
								{1./36., 1./24., 1./72., 0, 0, -1./12., -1./24., 1./12., 1./24., -1./24., -1./24., 0, 0, 0, -1./4., 0, 0, -1./8., -1./8.}};

															
															
const double clConstants::clLB::Cs2 = 1./3.;
const double clConstants::clLB::Cs = sqrt(clConstants::clLB::Cs2);

clVector3Dd clConstants::clLB::UX0 = clVector3Dd(LB_UX0_X, LB_UX0_Y, LB_UX0_Z);
clVector3Dd clConstants::clLB::UX1 = clVector3Dd(LB_UX1_X, LB_UX1_Y, LB_UX1_Z);
clVector3Dd clConstants::clLB::UY0 = clVector3Dd(LB_UY0_X, LB_UY0_Y, LB_UY0_Z);
clVector3Dd clConstants::clLB::UY1 = clVector3Dd(LB_UY1_X, LB_UY1_Y, LB_UY1_Z);
clVector3Dd clConstants::clLB::UZ0 = clVector3Dd(LB_UZ0_X, LB_UZ0_Y, LB_UZ0_Z);
clVector3Dd clConstants::clLB::UZ1 = clVector3Dd(LB_UZ1_X, LB_UZ1_Y, LB_UZ1_Z);

const double clConstants::clLB::RoX0 = LB_RO_X0;
const double clConstants::clLB::RoX1 = LB_RO_X1;
const double clConstants::clLB::RoY0 = LB_RO_Y0;
const double clConstants::clLB::RoY1 = LB_RO_Y1;
const double clConstants::clLB::RoZ0 = LB_RO_Z0;
const double clConstants::clLB::RoZ1 = LB_RO_Z1;

BOOL clConstants::clLB::bcWallX0 = LB_BC_WALL_X0;
BOOL clConstants::clLB::bcWallY0 = LB_BC_WALL_Y0;
BOOL clConstants::clLB::bcWallZ0 = LB_BC_WALL_Z0;
BOOL clConstants::clLB::bcWallX1 = LB_BC_WALL_X1;
BOOL clConstants::clLB::bcWallY1 = LB_BC_WALL_Y1;
BOOL clConstants::clLB::bcWallZ1 = LB_BC_WALL_Z1;

BOOL clConstants::clLB::bcPressX0 = LB_BC_PRESS_X0;
BOOL clConstants::clLB::bcPressY0 = LB_BC_PRESS_Y0;
BOOL clConstants::clLB::bcPressZ0 = LB_BC_PRESS_Z0;
BOOL clConstants::clLB::bcPressX1 = LB_BC_PRESS_X1;
BOOL clConstants::clLB::bcPressY1 = LB_BC_PRESS_Y1;
BOOL clConstants::clLB::bcPressZ1 = LB_BC_PRESS_Z1;

BOOL clConstants::clLB::bcFreeX0 = LB_BC_FREE_X0;
BOOL clConstants::clLB::bcFreeY0 = LB_BC_FREE_Y0;
BOOL clConstants::clLB::bcFreeZ0 = LB_BC_FREE_Z0;
BOOL clConstants::clLB::bcFreeX1 = LB_BC_FREE_X1;
BOOL clConstants::clLB::bcFreeY1 = LB_BC_FREE_Y1;
BOOL clConstants::clLB::bcFreeZ1 = LB_BC_FREE_Z1;

BOOL clConstants::clLB::bcSymmetryX0 = LB_BC_SYMMETRY_X0;
BOOL clConstants::clLB::bcSymmetryY0 = LB_BC_SYMMETRY_Y0;
BOOL clConstants::clLB::bcSymmetryZ0 = LB_BC_SYMMETRY_Z0;
BOOL clConstants::clLB::bcSymmetryX1 = LB_BC_SYMMETRY_X1;
BOOL clConstants::clLB::bcSymmetryY1 = LB_BC_SYMMETRY_Y1;
BOOL clConstants::clLB::bcSymmetryZ1 = LB_BC_SYMMETRY_Z1;

BOOL clConstants::clLB::bcCoarseX0 = LB_BC_COARSE_X0;
BOOL clConstants::clLB::bcCoarseY0 = LB_BC_COARSE_Y0;
BOOL clConstants::clLB::bcCoarseZ0 = LB_BC_COARSE_Z0;
BOOL clConstants::clLB::bcCoarseX1 = LB_BC_COARSE_X1;
BOOL clConstants::clLB::bcCoarseY1 = LB_BC_COARSE_Y1;
BOOL clConstants::clLB::bcCoarseZ1 = LB_BC_COARSE_Z1;

BOOL clConstants::clLB::bcPeriodicX = LB_BC_PERIODIC_X;
BOOL clConstants::clLB::bcPeriodicY = LB_BC_PERIODIC_Y;
BOOL clConstants::clLB::bcPeriodicZ = LB_BC_PERIODIC_Z;


int clConstants::clTr::nN = TRC_NUM_TRACERS;
double clConstants::clTr::D0 = TR_D0;
double clConstants::clTr::TRUNC = TR_TRUNC;
double clConstants::clTr::K = K_coeff;
int clConstants::clTr::Stop_dist_x = STOP_DIST_X;
int clConstants::clTr::remain = TRC_NUM_TRACERS;

////////////////////////////////////////////////////////////////////
void clConstants::stepSignal()
{
	clProblem::o.step();	
}

void clConstants::fieldSignal()
{
	clProblem::o.field();	
	clProblem::d.field();	
#if !defined(_WIN32)
#endif //_WIN32
}
void clConstants::save2file()
{
	ls.save2file();
	//lb.save2file();

	if(create_out_parameter() == FALSE)
		return;

	saveparameter("Version", VERSION);
	saveparameter("OpSys", _OPERATION_SYSTEM);
	saveparameter_line();

	saveparameter("c_Re1", Re1);
	saveparameter("c_Pe1", Pe1);
	saveparameter("c_Gx", G.x);
	saveparameter("c_Gy", G.y);
	saveparameter("c_Gz", G.z);
	saveparameter("c_L0", L0);
	saveparameter("c_Ro0", Ro0);
	saveparameter("c_nu0", nu0);
	saveparameter("c_T0", T0);
	saveparameter("c_eps", eps);
	saveparameter_line();

	saveparameter("c_sL", sL);
	saveparameter("c_sM", sM);
	saveparameter("c_sT", sT);
	saveparameter("c_sN", sN);
	saveparameter_line();

	saveparameter("c_t_flContinue", t.flContinue);
	saveparameter("c_t_flOutputDump", t.flOutputDump);
	saveparameter("c_t_dTlb", t.dTlb);
	saveparameter("c_t_dTls", t.dTls);
	saveparameter("c_t_dTtr", t.dTtr);
	saveparameter("c_t_Tstop", t.StopTime);
	saveparameter("c_t_Tcurrent", t.Time);
	saveparameter("c_t_Ttotal", t.StartTime + t.Time);
	saveparameter("c_t_Tdumpstart", t.DumpStartTime);
	saveparameter("c_t_Tdumpstep", t.DumpTimeStep);
	saveparameter("c_t_Tdumpstepsave", t.DumpTimeStepSave);
	saveparameter("c_t_Tcurrentsavestart", t.CurrentSaveStartTime);
	saveparameter("c_t_Tcurrentsavestep", t.CurrentSaveTimeStep);
	saveparameter_line();

	saveparameter("c_t_oscAmpl1", t.oscAmpl1);
	saveparameter("c_t_oscPeriod1", t.oscPeriod1);
	saveparameter("c_t_oscAmpl2", t.oscAmpl2);
	saveparameter("c_t_oscPeriod2", t.oscPeriod2);
	saveparameter("c_t_oscRatio", t.oscRatio);
	saveparameter_line();

	saveparameter("c_tr_nN", tr.nN);
	saveparameter("c_tr_remain", tr.remain);
	saveparameter("c_tr_D0", tr.D0);
	saveparameter_line();

	saveparameter("c_ls_nN", ls.nN);
	saveparameter("c_ls_nL", ls.nL);
	saveparameter("c_ls_nBN", ls.nBN);
	saveparameter("c_ls_nBL", ls.nBL);
	saveparameter("c_ls_nPupd", ls.nPupd);
	saveparameter("c_ls_nAN", ls.nAN);
	//saveparameter("c_ls_Lb", ls.Lb);
	saveparameter_line();

	//saveparameter("c_ls_nTp", ls.nTp);
	saveparameter("c_ls_nTa", ls.nTa);
	saveparameter("c_ls_nTr", ls.nTr);
	saveparameter("c_ls_nTs", ls.nTs);
	saveparameter("c_ls_nO", ls.nO);
	saveparameter("c_ls_nOp", ls.nOp);
	saveparameter("c_ls_nPN", ls.nPN);
	saveparameter("c_ls_nBFmax", ls.nBFmax);
	saveparameter_line();

	saveparameter("c_ls_Ds0", ls.Ds0);
	saveparameter("c_ls_Ds1", ls.Ds0);
	saveparameter("c_ls_Ds2", ls.Ds0);
	saveparameter_line();

	saveparameter("c_ls_FwDe0", ls.FwDe0);
	saveparameter("c_ls_FwDe1", ls.FwDe1);
	saveparameter("c_ls_FwDe2", ls.FwDe2);
	saveparameter("c_ls_FwDeR", ls.FwDeR);
	saveparameter("c_ls_FwRe", ls.FwRe);
	saveparameter("c_ls_FwKe", ls.FwKe);
	saveparameter("c_ls_FwDeC", ls.FwDeC);
	saveparameter("c_ls_FwReC", ls.FwReC);
	saveparameter_line();

	saveparameter("c_ls_Symax0", ls.Symax0);
	saveparameter("c_ls_Symin0", ls.Symin0);
	saveparameter("c_ls_Symax1", ls.Symax1);
	saveparameter("c_ls_Symin1", ls.Symin1);
	saveparameter("c_ls_Symax2", ls.Symax2);
	saveparameter("c_ls_Symin2", ls.Symin2);
	saveparameter_line();

	saveparameter("c_ls_FwDeX0", ls.FwDeX0);
	saveparameter("c_ls_FwDeRX0", ls.FwDeRX0);
	saveparameter("c_ls_FwDeX1", ls.FwDeX1);
	saveparameter("c_ls_FwDeRX1", ls.FwDeRX1);
	saveparameter("c_ls_FwDeY0", ls.FwDeY0);
	saveparameter("c_ls_FwDeRY0", ls.FwDeRY0);
	saveparameter("c_ls_FwDeY1", ls.FwDeY1);
	saveparameter("c_ls_FwDeRY1", ls.FwDeRY1);
	saveparameter("c_ls_FwDeZ0", ls.FwDeZ0);
	saveparameter("c_ls_FwDeRZ0", ls.FwDeRZ0);
	saveparameter("c_ls_FwDeZ1", ls.FwDeZ1);
	saveparameter("c_ls_FwDeRZ1", ls.FwDeRZ1);
	saveparameter_line();

	saveparameter("c_ls_Cini0x", ls.Cini0.x);
	saveparameter("c_ls_Cini0y", ls.Cini0.y);
	saveparameter("c_ls_Cini0z", ls.Cini0.z);
	saveparameter("c_ls_Cini1x", ls.Cini1.x);
	saveparameter("c_ls_Cini1y", ls.Cini1.y);
	saveparameter("c_ls_Cini1z", ls.Cini1.z);
	saveparameter("c_ls_Cini2x", ls.Cini2.x);
	saveparameter("c_ls_Cini2y", ls.Cini2.y);
	saveparameter("c_ls_Cini2z", ls.Cini2.z);
	saveparameter("c_ls_Kini0", ls.Kini0);
	saveparameter("c_ls_Kini1", ls.Kini1);
	saveparameter("c_ls_Kini2", ls.Kini2);
	saveparameter("c_ls_Mini0", ls.Mini0);
	saveparameter("c_ls_Mini1", ls.Mini1);
	saveparameter("c_ls_Mini2", ls.Mini2);
	saveparameter("c_ls_L0ini0", ls.L0ini0);
	saveparameter("c_ls_L0ini1", ls.L0ini1);
	saveparameter("c_ls_L0ini2", ls.L0ini2);
	saveparameter("c_ls_AKini0", ls.AKini0);
	saveparameter("c_ls_AKini1", ls.AKini1);
	saveparameter("c_ls_AKini2", ls.AKini2);
	saveparameter_line();

	for(int i = 0; i < ls.nO; i++)
	{
		char Buf[80];
		sprintf(Buf, "vls_avr_CGx%d", i);
		saveparameter(Buf, vls.avr.CG(i).x);
		sprintf(Buf, "vls_avr_CGy%d", i);
		saveparameter(Buf, vls.avr.CG(i).y);
		sprintf(Buf, "vls_avr_CGz%d", i);
		saveparameter(Buf, vls.avr.CG(i).z);
		sprintf(Buf, "vls_avr_Vx%d", i);
		saveparameter(Buf, vls.avr.V(i).x);
		sprintf(Buf, "vls_avr_Vy%d", i);
		saveparameter(Buf, vls.avr.V(i).y);
		sprintf(Buf, "vls_avr_Vz%d", i);
		saveparameter(Buf, vls.avr.V(i).z);
		sprintf(Buf, "vls_avr_Ox%d", i);
		saveparameter(Buf, vls.avr.O(i).x);
		sprintf(Buf, "vls_avr_Oy%d", i);
		saveparameter(Buf, vls.avr.O(i).y);
		sprintf(Buf, "vls_avr_Oz%d", i);
		saveparameter(Buf, vls.avr.O(i).z);
		sprintf(Buf, "vls_avr_Ek%d", i);
		saveparameter(Buf, vls.avr.Ek(i));
		sprintf(Buf, "vls_avr_Ec%d", i);
		saveparameter(Buf, vls.avr.Ec(i));
		sprintf(Buf, "vls_avr_Ev%d", i);
		saveparameter(Buf, vls.avr.Ev(i));
		sprintf(Buf, "vls_avr_Ravr%d", i);
		saveparameter(Buf, vls.avr.Ravr(i));
		sprintf(Buf, "vls_avr_Rmin%d", i);
		saveparameter(Buf, vls.avr.Rmin(i));
		sprintf(Buf, "vls_avr_Rmax%d", i);
		saveparameter(Buf, vls.avr.Rmax(i));
		sprintf(Buf, "vls_avr_Smin%d", i);
		saveparameter(Buf, vls.avr.Smin(i));
		sprintf(Buf, "vls_avr_Smax%d", i);
		saveparameter(Buf, vls.avr.Smax(i));
		saveparameter_line();
	}

	saveparameter("c_lb_nX", lb.nX);
	saveparameter("c_lb_nY", lb.nY);
	saveparameter("c_lb_nZ", lb.nZ);
	saveparameter("c_lb_nF", lb.nF - 1);
	saveparameter("c_lb_nD", lb.nD);
	saveparameter("c_lb_Cs2", lb.Cs2);
	saveparameter_line();

	saveparameter("c_lb_CoarseStep", lb.CoarseStep);
	saveparameter("c_lb_nXc", lb.nXc);
	saveparameter("c_lb_nYc", lb.nYc);
	saveparameter("c_lb_nZc", lb.nZc);
	saveparameter("c_lb_Xc0", lb.Xc0);
	saveparameter("c_lb_Yc0", lb.Yc0);
	saveparameter("c_lb_Zc0", lb.Zc0);
	saveparameter_line();

	saveparameter("vlb_M_nS", vlb.M.nS);
	saveparameter("vlb_M_nSfull", vlb.M.nSf);
	saveparameter_line();

	saveparameter("c_lb_UX0x", lb.UX0.x);
	saveparameter("c_lb_UX0y", lb.UX0.y);
	saveparameter("c_lb_UX0z", lb.UX0.z);
	saveparameter("c_lb_UX1x", lb.UX1.x);
	saveparameter("c_lb_UX1y", lb.UX1.y);
	saveparameter("c_lb_UX1z", lb.UX1.z);
	saveparameter("c_lb_UY0x", lb.UY0.x);
	saveparameter("c_lb_UY0y", lb.UY0.y);
	saveparameter("c_lb_UY0z", lb.UY0.z);
	saveparameter("c_lb_UY1x", lb.UY1.x);
	saveparameter("c_lb_UY1y", lb.UY1.y);
	saveparameter("c_lb_UY1z", lb.UY1.z);
	saveparameter("c_lb_UZ0x", lb.UZ0.x);
	saveparameter("c_lb_UZ0y", lb.UZ0.y);
	saveparameter("c_lb_UZ0z", lb.UZ0.z);
	saveparameter("c_lb_UZ1x", lb.UZ1.x);
	saveparameter("c_lb_UZ1y", lb.UZ1.y);
	saveparameter("c_lb_UZ1z", lb.UZ1.z);
	saveparameter_line();

	saveparameter("c_lb_bcWallX0", lb.bcWallX0);
	saveparameter("c_lb_bcWallX1", lb.bcWallX1);
	saveparameter("c_lb_bcWallY0", lb.bcWallY0);
	saveparameter("c_lb_bcWallY1", lb.bcWallY1);
	saveparameter("c_lb_bcWallZ0", lb.bcWallZ0);
	saveparameter("c_lb_bcWallZ1", lb.bcWallZ1);
	saveparameter("c_lb_bcFreeX0", lb.bcFreeX0);
	saveparameter("c_lb_bcFreeX1", lb.bcFreeX1);
	saveparameter("c_lb_bcFreeY0", lb.bcFreeY0);
	saveparameter("c_lb_bcFreeY1", lb.bcFreeY1);
	saveparameter("c_lb_bcFreeZ0", lb.bcFreeZ0);
	saveparameter("c_lb_bcFreeZ1", lb.bcFreeZ1);
	saveparameter("c_lb_bcSymmetryX0", lb.bcSymmetryX0);
	saveparameter("c_lb_bcSymmetryX1", lb.bcSymmetryX1);
	saveparameter("c_lb_bcSymmetryY0", lb.bcSymmetryY0);
	saveparameter("c_lb_bcSymmetryY1", lb.bcSymmetryY1);
	saveparameter("c_lb_bcSymmetryZ0", lb.bcSymmetryZ0);
	saveparameter("c_lb_bcSymmetryZ1", lb.bcSymmetryZ1);
	saveparameter("c_lb_bcCoarseX0", lb.bcCoarseX0);
	saveparameter("c_lb_bcCoarseX1", lb.bcCoarseX1);
	saveparameter("c_lb_bcCoarseY0", lb.bcCoarseY0);
	saveparameter("c_lb_bcCoarseY1", lb.bcCoarseY1);
	saveparameter("c_lb_bcCoarseZ0", lb.bcCoarseZ0);
	saveparameter("c_lb_bcCoarseZ1", lb.bcCoarseZ1);
	saveparameter("c_lb_bcPeriodicX", lb.bcPeriodicX);
	saveparameter("c_lb_bcPeriodicY", lb.bcPeriodicY);
	saveparameter("c_lb_bcPeriodicZ", lb.bcPeriodicZ);
	saveparameter_line();


	for(int i = 1; i < c.lb.nF; i++)
	{
		char Buf[80];
		sprintf(Buf, "c_lb_lamda%d", i);
		saveparameter(Buf, lb.lamda[i]);
		sprintf(Buf, "c_lb_lamdab%d", i);
		saveparameter(Buf, lb.lamdab[i]);
		sprintf(Buf, "c_lb_tau%d", i);
		saveparameter(Buf, lb.tau[i]);
		sprintf(Buf, "c_lb_Ro%d", i);
		saveparameter(Buf, lb.Ro[i]);
		sprintf(Buf, "c_lb_Fx%d", i);
		saveparameter(Buf, lb.F[i].x);
		sprintf(Buf, "c_lb_Fy%d", i);
		saveparameter(Buf, lb.F[i].y);
		sprintf(Buf, "c_lb_Fz%d", i);
		saveparameter(Buf, lb.F[i].z);
		saveparameter_line();
	}

	close_out_parameter();
	return;
}


void clConstants::clTime::update()
{
	if(NoUpdateSignal)
	{
        NoUpdateSignal = FALSE;
		return;
	}
	
	Time = TimeN;

#ifdef LS_SOLVER
	TimeN += dTls;
#else
	TimeN += dTlb;
#endif //LS_SOLVER
	//TimeN++;
}


void clConstants::clLS::clN::ini()
{
	zeros(c.ls.nN, c.ls.nL);
	Sind.zeros(c.ls.nN, c.ls.nL);
	NL.zeros(c.ls.nN);
	if(load("load" SLASH "lsn") == FALSE)
	{
//		v.error("LS");
	}
//	fl.zeros(c.ls.nN, c.ls.nL);
	//if(loadtxt("load" SLASH "lsn") == FALSE)
	//{
	//	for(int i = 0; i < c.ls.nN; i++)
	//		for(int j = 0; j < c.ls.nL; j++)
	//			if(c.ls.N(i,j) >= 0)
	//				fl(i,j) = TRUE;
	//			else
	//				fl(i,j) = FALSE;
	
	setNL();

	save2file();
}

void clConstants::clLS::clN::test_links()
{
	for(int i = 0; i < c.ls.nN; i++)
	{
		for(int j = 0; j < c.ls.nL; j++)
		{
			int n = el(i,j);
			if(n >= 0 && !test_node(Sind(i,j)))
				remove_links(i,n);
		}
		if(NL(i) < NUMBER_DIMENSIONS)
			remove_all_links(i);

	}

}
BOOL clConstants::clLS::clN::test_node(int n)
{
	for(int i = 0; i < nSf; i++)
		if(S(n,i) >= 0)
		{
			if(vlb.M.Sfdel(S(n,i)) == FALSE)
				return TRUE;
			else
				S(n,i) = -1;
		}
	return FALSE;
}
//BOOL clConstants::clLS::clN::test_node(int n)
//{
//	BOOL ret = FALSE;
//	for(int i = 0; i < nSf; i++)
//		if(S(n,i) >= 0)
//		{
//			if(vlb.M.Sfdel(S(n,i)) == FALSE)
//				ret = TRUE;
//			else
//				S(n,i) = -1;
//		}
//	return ret;
//}
void clConstants::clLS::clN::iniS()
{
	nS = 0;
	nSf = LS_MAX_SF;
//	if(Sind.load("load" SLASH "lsnsind") != FALSE)
//	{
//		for(int i = 0; i < c.ls.nN; i++)
//			for(int j = 0; j < c.ls.nL; j++)
//				if(Sind(i,j) > nS + 1)
//					nS = Sind(i,j) + 1;
//
//#ifdef LS_RUPTURE
//		Sind.save("lsnsind");
//#else //LS_RUPTURE
//		Sindtxt.save("lsnsind");
//#endif //LS_RUPTURE
//		S.zeros(nS, nSf);
//		if(S.load("load" SLASH "lsns") != FALSE)
//			return;
//	}

	Sind.fill(-1);
	for(int i = 0; i < c.ls.nN; i++)
		for(int j = 0; j < c.ls.nL; j++)
			if(Sind(i,j) < 0 && c.ls.N(i,j) >= 0)
			{
				for(int k = 0; k < c.ls.nL; k++)
					if(c.ls.N(c.ls.N(i,j),k) == i)
					{
						Sind(c.ls.N(i,j),k) = nS;
						break;
					}
				Sind(i,j) = nS++;
			}

	S.zeros(nS, nSf);
	S.fill(-1);
	for(int i = 0; i < c.ls.nN; i++)
		for(int j = 0; j < c.ls.nL; j++)
			if(c.ls.N(i,j) >= 0)
				if(findS1(i, c.ls.N(i,j), Sind(i,j)) == 0)
					findS2(i, c.ls.N(i,j), Sind(i,j));

//#ifdef LS_RUPTURE
//	Sind.save("lsnsind");
//#else //LS_RUPTURE
	Sind.savetxt("lsnsind");
//#endif //LS_RUPTURE
#ifdef LS_DEBUG
	S.savetxt("lsns");
#endif
};
int clConstants::clLS::clN::findS1(int n0, int n1, int ns)
{
	int ind = 0;
	for(int i = 0; i < vlb.M.nSf; i++)
	{
		if((vlb.M.Sf(i,0) == n0 || vlb.M.Sf(i,1) == n0 || vlb.M.Sf(i,2) == n0 || vlb.M.Sf(i,3) == n0)
				&& (vlb.M.Sf(i,0) == n1 || vlb.M.Sf(i,1) == n1 || vlb.M.Sf(i,2) == n1 || vlb.M.Sf(i,3) == n1))
			S(ns, ind++) = i;
		if(ind >= nSf)
			printf("Error!!! Index is larger than LS_MAX_SF\n");
	}
	return ind;
};
int clConstants::clLS::clN::findS2(int n0, int n1, int ns)
{
	int ind0 = 0, ind1 = 0, ind = 0;
	clArray2Di S0, S1;
	clArray1Di Sind0, Sind1;

	S0.zeros(nSf, NUMBER_DIMENSIONS + 1);
	S1.zeros(nSf, NUMBER_DIMENSIONS + 1);
	Sind0.zeros(nSf);
	Sind1.zeros(nSf);
	S0.fill(-1);
	S1.fill(-1);
	Sind0.fill(-1);
	Sind1.fill(-1);

	for(int i = 0; i < vlb.M.nSf; i++)
	{
		if(vlb.M.Sf(i,0) == n0 || vlb.M.Sf(i,1) == n0 || vlb.M.Sf(i,2) == n0 || vlb.M.Sf(i,3) == n0)
		{
			S0(ind0, 0) = vlb.M.Sf(i,0);
			S0(ind0, 1) = vlb.M.Sf(i,1);
			S0(ind0, 2) = vlb.M.Sf(i,2);
			S0(ind0, 3) = vlb.M.Sf(i,3);
			Sind0(ind0) = i;
			ind0++;
			if(ind0 >= nSf)
				printf("Error!!! Index is larger than LS_MAX_SF\n");
		}
		if(vlb.M.Sf(i,0) == n1 || vlb.M.Sf(i,1) == n1 || vlb.M.Sf(i,2) == n1 || vlb.M.Sf(i,3) == n1)
		{
			S1(ind1, 0) = vlb.M.Sf(i,0);
			S1(ind1, 1) = vlb.M.Sf(i,1);
			S1(ind1, 2) = vlb.M.Sf(i,2);
			S1(ind1, 3) = vlb.M.Sf(i,3);
			Sind1(ind1) = i;
			ind1++;
			if(ind1 >= nSf)
				printf("Error!!! Index is larger than LS_MAX_SF\n");
		}
	}
#ifdef LS_DEBUG
	//S0.savetxt("lsns0");
	//Sind0.savetxt("lsnsind0");
	//S1.savetxt("lsns1");
	//Sind1.savetxt("lsnsind1");
#endif

	double contacts;
	for(int i = 0; i < ind0; i++)
	{
		for(int j = 0; j < ind1; j++)
		{
			if(Sind1(j) < 0&& Sind0(i) < 0)
				continue;

			contacts = 0;
			for(int k = 0; k < (NUMBER_DIMENSIONS + 1); k++)
				for(int l = 0; l < (NUMBER_DIMENSIONS + 1); l++)
					if(S0(i,k) == S1(j,l))
						contacts++;
			
			if(contacts >= 3)
			{
				if(Sind0(i) >= 0)
				{
					S(ns, ind++) = Sind0(i);
					Sind0(i) = -1;
				}
				if(Sind1(j) >= 0)
				{
					S(ns, ind++) = Sind1(j);
					Sind1(j) = -1;
				}
				if(ind >= nSf)
					printf("Error!!! Index is larger than LS_MAX_SF\n");
			}

		}
	}

#ifdef LS_DEBUG
	//Sind0.savetxt("lsnsind0p");
	//Sind1.savetxt("lsnsind1p");
	//S.savetxt("lsns");
#endif

	if(ind > 0)
		return ind;

	for(int i = 0; i < ind0; i++)
	{
		for(int j = 0; j < ind1; j++)
		{
			if(Sind1(j) < 0 && Sind0(i) < 0)
				continue;

			contacts = 0;
			for(int k = 0; k < (NUMBER_DIMENSIONS + 1); k++)
				for(int l = 0; l < (NUMBER_DIMENSIONS + 1); l++)
					if(S0(i,k) == S1(j,l))
						contacts++;
			
			if(contacts >= 2)
			{
				if(Sind0(i) >= 0)
				{
					S(ns, ind++) = Sind0(i);
					Sind0(i) = -1;
				}
				if(Sind1(j) >= 0)
				{
					S(ns, ind++) = Sind1(j);
					Sind1(j) = -1;
				}
				if(ind >= nSf)
					printf("Error!!! Index is larger than LS_MAX_SF\n");
			}

		}
	}

#ifdef LS_DEBUG
	//Sind0.savetxt("lsnsind0p");
	//Sind1.savetxt("lsnsind1p");
	//S.savetxt("lsns");
#endif

	if(ind > 0)
		return ind;

	for(int i = 0; i < ind0; i++)
	{
		for(int j = 0; j < ind1; j++)
		{
			if(Sind1(j) < 0)
				continue;

			for(int k = 0; k < (NUMBER_DIMENSIONS + 1); k++)
				for(int l = 0; l < (NUMBER_DIMENSIONS + 1); l++)
					if(S0(i,k) == S1(j,l))
					{
						if(Sind0(i) >= 0)
						{
							S(ns, ind++) = Sind0(i);
							Sind0(i) = -1;
						}
						if(Sind1(j) >= 0)
						{
							S(ns, ind++) = Sind1(j);
							Sind1(j) = -1;
						}
						if(ind >= nSf)
							printf("Error!!! Index is larger than LS_MAX_SF\n");
					}
		}
	}
#ifdef LS_DEBUG
	Sind0.savetxt("lsnsind0p");
	Sind1.savetxt("lsnsind1p");
	S.savetxt("lsns");
#endif
	return ind;
};
void clConstants::clLS::clM::ini()
{
	zeros(c.ls.nN);
	if(loadtxt("load" SLASH "lsm") == FALSE)
	{
//		v.error("Cannot initialze Cx in LS");
	}
	if(c.t.flContinue == FALSE)
	{
		for(int i = 0; i < c.ls.nN; i++)
		{
			if(c.ls.O(i) == 0)
				el(i) *= c.ls.Mini0;
			if(c.ls.O(i) == 1)
				el(i) *= c.ls.Mini1;
			//Modified
			if(c.ls.O(i) >= 2)
			//Modified
				el(i) *= c.ls.Mini2;
		}
	}

	total = 0.;
	for(int i = 0; i < c.ls.nN; i++)
        total += el(i);

	//totalp = 0.;
	//for(int i = 0; i < c.ls.nTp; i++)
	//	totalp += el(c.ls.Tp(i));

	totala = 0.;
	for(int i = 0; i < c.ls.nTa; i++)
        totala += el(c.ls.Ta(i));

	O.zeros(c.ls.nO);
	for(int i = 0; i < c.ls.nN; i++)
		O(c.ls.O(i)) += el(i);

	save2file();
}

void clConstants::clLS::clDs::ini()
{
	zeros(c.ls.nN, c.ls.nL);
#ifdef LS_DISSIPATION
	if(loadtxt("load" SLASH "lsds") == FALSE)
	{
		for(int i = 0; i < c.ls.nN; i++)
			for(int j = 0; j < c.ls.nL; j++)
			{
				el(i,j) = c.ls.Ds0;
				if(c.ls.O(i) == 1)
					el(i,j) = c.ls.Ds1;
				if(c.ls.O(i) >= 2)
					el(i,j) = c.ls.Ds2;
			}
	}
#else
	c.ls.Ds0 = 0.;
#endif //LS_DISSIPATION
	save2file();
}
void clConstants::clLS::iniA()
{
	AN.loadnew("load" SLASH "lsan", 4);

	c.ls.nAN = c.ls.AN.getSize1();
	if(c.ls.nAN == 0)
		return;

	AA0.zeros(c.ls.nAN);

#ifdef LS_DEBUG
#endif //LS_DEBUG
	c.ls.AN.savetxt("lsan");
}
void clConstants::clLS::iniAA0()
{
	if(c.ls.nAN == 0)
		return;
	if(AA0.load("load" SLASH "lsaa0") == FALSE)
	{
		for(int i = 0; i < c.ls.nAN; i++)
		{
			int i0 = c.ls.AN(i,0), i1 = c.ls.AN(i,1), i2 = c.ls.AN(i,2);

			clVector3Dd d1 = vls.C(i0) - vls.C(i1);
			clVector3Dd d2 = vls.C(i2) - vls.C(i1);
			double cs = d1 * d2 / (d1.abs() * d2.abs());
			if(cs > 1.) 
				cs = 1.;
			else if(cs < -1.) 
				cs = -1.;

			AA0(i) = acos(cs);
		}
	}

#ifdef LS_DEBUG
#endif //LS_DEBUG
	c.ls.AA0.savetxt("lsaa0");
}

void clConstants::clLS::clK::ini()
{
	if(c.ls.nN == 0)
		return;

	zeros(c.ls.nN, c.ls.nL);
	E.zeros(c.ls.nN, c.ls.nL);
	if(loadtxt("load" SLASH "lsk") == FALSE)
	{
		v.error("Cannot initialze lsk in LS");
	}
	if(c.t.flContinue == FALSE)
	{
		for(int i = 0; i < c.ls.nN; i++)
			for(int j = 0; j < c.ls.nL; j++)
			{
				if(c.ls.O(i) == 0)
					el(i,j) *= c.ls.Kini0;
				if(c.ls.O(i) == 1)
					el(i,j) *= c.ls.Kini1;
				if(c.ls.O(i) >= 2)
					el(i,j) *= c.ls.Kini2;
			}
	}

	save2file();
}

void clConstants::clLS::clAK::ini()
{
	if(c.ls.nAN == 0)
		return;

	zeros(c.ls.nAN);
	E.zeros(c.ls.nAN);
	if(loadtxt("load" SLASH "lsak") == FALSE)
	{
		v.error("Cannot initialze lsak in LS");
	}
	if(c.t.flContinue == FALSE)
	{
		for(int i = 0; i < c.ls.nAN; i++)
		{
			if(c.ls.O(c.ls.AN(i,1)) == 0)
				el(i) *= c.ls.AKini0;
			if(c.ls.O(c.ls.AN(i,1)) == 1)
				el(i) *= c.ls.AKini1;
			//Modified
			if(c.ls.O(c.ls.AN(i,1)) >= 2)
			//Modified
				el(i) *= c.ls.AKini2;
		}
	}

	save2file();
}

void clConstants::clLS::clL0::ini()
{
	zeros(c.ls.nN, c.ls.nL);
	L0ini.zeros(c.ls.nN, c.ls.nL);
	if(load("load" SLASH "lsl0") == FALSE)
	{
		for(int i = 0; i < c.ls.nN; i++)
		{
			for(int j = 0; j < c.ls.nL; j++)
			{
				int n = c.ls.N(i,j);
				if(n >= 0)
				{
					clVector3Dd d = vls.C(i) - vls.C(n);
					el(i,j) = d.r();
					if(c.ls.O(i) == 0)
						el(i,j) *= c.ls.L0ini0;
					if(c.ls.O(i) == 1)
						el(i,j) *= c.ls.L0ini1;
					//Modified
					if(c.ls.O(i) >= 2)
					//Modified
						el(i,j) *= c.ls.L0ini2;
				}
			}
		}
	}
	L0ini.copy(*this);
	save2file();
}
void clConstants::clLS::clL0::update()
{
	return;
	//for(int i = 0; i < c.ls.nN; i++)
	//	for(int j = 0; j < c.ls.nL; j++)
	//		if(c.ls.N(i,j) >= 0)
	//			el(i,j) *= (1. - 0.05 / STOP_TIME / 10.);
	for(int i = 0; i < c.ls.nN; i++)
	{
//		if(c.ls.O(i) > 0)
//			continue;
		double phase = vls.avr.CG(c.ls.O(i)).x / (double) c.nX * 2. * c.Pi;
		double osc = c.t.oscillator1(phase);
		int Tnode0 = c.ls.T(i);
//		for(int j = 0; j < c.ls.nL; j++)
		for(int j = 0; j < 6; j++)
		{
			int n1 = c.ls.N(i,j);
			if(n1 >= 0)
			{
				int Tnode1 = c.ls.T(n1);
				if(((Tnode0 & LS_CUSTOM_NODE_0) && (Tnode1 & LS_CUSTOM_NODE_1)) || 
					((Tnode0 & LS_CUSTOM_NODE_1) && (Tnode1 & LS_CUSTOM_NODE_0)))
					el(i,j) = L0ini(i,j) * (1. + osc);
				if(((Tnode0 & LS_CUSTOM_NODE_2) && (Tnode1 & LS_CUSTOM_NODE_3)) || 
					((Tnode0 & LS_CUSTOM_NODE_3) && (Tnode1 & LS_CUSTOM_NODE_2)))
					el(i,j) = L0ini(i,j) * (1. - osc);
			}
		}
	}
}

void clConstants::clLS::iniB()
{
	c.ls.BNind.zeros(c.ls.nN);
	c.ls.BNind.fill(-1);

	if(c.ls.BN.loadnew("load" SLASH "lsbn") == TRUE)
	{
		c.ls.nBN = c.ls.BN.getSize(1);

		for(int i = 0; i < c.ls.nBN; i++)
			c.ls.BNind(c.ls.BN(i)) = i;
	}
	else
	{
		int i, j, counter;
		c.ls.nBN = 0;
		for(i = 0; i < c.ls.nN; i++)
		{
			counter = 0;
			for(j = 0; j < c.ls.nL; j++)
				if(c.ls.N(i,j) < 0)
					counter++;

			if(counter > 2)
			{
				c.ls.BNind(i) = c.ls.nBN;
				c.ls.nBN++;
			}
		}

		c.ls.BN.zeros(c.ls.nBN);
		
		for(i = 0; i < c.ls.nN; i++)
			if(c.ls.BNind(i) >= 0)
				c.ls.BN(c.ls.BNind(i)) = i;

	}


	c.ls.BS.zeros(c.ls.nBN);
	c.ls.BS.load("load" SLASH "lsbs");

	c.ls.BF.zeros(c.ls.nBN);
	if(c.ls.BF.load("load" SLASH "lsbf") == FALSE)
		c.ls.BF.fill(1);

	for(int i = 0; i < c.ls.nBN; i++)
		if(c.ls.BF(i) > c.ls.nBFmax)
			c.ls.nBFmax = c.ls.BF(i);

#ifdef LS_RUPTURE
	c.ls.BN.save("lsbn");
	c.ls.BS.save("lsbs");
	c.ls.BF.save("lsbf");
#else //LS_RUPTURE
	c.ls.BN.savetxt("lsbn");
	c.ls.BS.savetxt("lsbs");
	c.ls.BF.savetxt("lsbf");
#endif //LS_RUPTURE
	c.ls.BNind.savetxt("lsbnind");
}
void clConstants::clLS::removeB(int n)
{
	int bn = c.ls.BNind(n);
	if(bn < 0)
		return;

	c.ls.BNind(n) = -1;
	c.ls.T(n) &= !LS_OUTER_NODE;

	for(int i = 0; i < c.ls.nN; i++)
	{
		if(c.ls.BNind(i) > bn)
			c.ls.BNind(i)--;
	}

	for(int i = bn; i < c.ls.nBN-1; i++)
	{
		c.ls.BN(i) = c.ls.BN(i+1);
		c.ls.BF(i) = c.ls.BF(i+1);
		c.ls.BS(i) = c.ls.BS(i+1);
	}

	c.ls.nBN--;
	c.ls.BN.reallocate(c.ls.nBN);
	c.ls.BF.reallocate(c.ls.nBN);
	c.ls.BS.reallocate(c.ls.nBN);

#ifdef LS_DEBUG
	c.ls.BN.savetxt("lsbn");
	c.ls.BS.savetxt("lsbs");
	c.ls.BF.savetxt("lsbf");
	c.ls.BNind.savetxt("lsbnind");
#endif //LS_DEBUG
}
void clConstants::clLS::addB(int bn, int bf, int bs)
{
	int bnind = c.ls.BNind(bn);
	if(bnind >= 0)
		return;

	c.ls.BNind(bn) = c.ls.nBN;
	c.ls.nBN++;
	c.ls.BN.reallocate(c.ls.nBN);
	c.ls.BF.reallocate(c.ls.nBN);
	c.ls.BS.reallocate(c.ls.nBN);

	c.ls.BN(c.ls.nBN-1) = bn;
	c.ls.BF(c.ls.nBN-1) = bf;
	c.ls.BS(c.ls.nBN-1) = bs;

	c.ls.T(bn) |= LS_OUTER_NODE;

#ifdef LS_DEBUG
	c.ls.BN.savetxt("lsbn");
	c.ls.BS.savetxt("lsbs");
	c.ls.BF.savetxt("lsbf");
	c.ls.BNind.savetxt("lsbnind");
#endif //LS_DEBUG
}

void clConstants::clLS::iniBL()
{
	c.ls.BLc.zeros(c.ls.nN, c.ls.nL);
//	c.ls.BLc0.zeros(c.ls.nN, c.ls.nL);
	c.ls.BLn.allocate(c.ls.nN, c.ls.nL);
	c.ls.BLFb.zeros(c.ls.nN, NUMBER_DIMENSIONS);


	if(c.ls.BL.loadnew("load" SLASH "lsbl",  NUMBER_DIMENSIONS) == TRUE)
		c.ls.nBL = c.ls.BL.getSize(1);
	else
	{
		int n = c.ls.nBN;
		c.ls.BL.zeros(n, NUMBER_DIMENSIONS);
		c.ls.nBL = 0;
		for(int i = 0; i < vlb.M.nSf; i++)
		{
			if(n < c.ls.nBL + 4)
			{
				n += c.ls.nBN;
				c.ls.BL.reallocate(n, NUMBER_DIMENSIONS);
			}
			testBL(i, 0, 1, 2);
			testBL(i, 0, 1, 3);
			testBL(i, 0, 2, 3);
			testBL(i, 1, 2, 3);
		}
		removeextraBL();

		for(int i0 = 0; i0 < c.ls.nN; i0++)
			if(c.ls.T(i0) & LS_FLAT_NODE && c.ls.BNind(i0) >= 0)
			{
				int n0 = i0;
				for(int i1 = 0; i1 < c.ls.nL; i1++)
				{
					int n1 = c.ls.N(i0,i1);
					if(n1 > 0 && c.ls.T(n1) & LS_FLAT_NODE && c.ls.BNind(n1) >= 0)
						for(int i2 = i1; i2 < c.ls.nL; i2++)
						{
							int n2 = c.ls.N(i0,i2);
							if(n2 > 0 && c.ls.T(n2) & LS_FLAT_NODE && c.ls.BNind(n2) >= 0 && c.ls.N.find(n1, n2) != -1)
								testflatBL(n0, n1, n2);
						}
				}
			}

		c.ls.BL.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);
	}
	c.ls.BLFb.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);

#ifdef BIO_FLOW
	c.ls.Inout.zeros(c.ls.nBL);
	//c.ls.Inout.loadtxt("load" SLASH "lsinout");
	if(c.ls.Inout.loadtxt("load" SLASH "lsinout") == FALSE)
	{
		//c.ls.Inout.zeros(c.ls.nBL);
		for(int i = 0; i < c.ls.nBL; i++)
		{
			testInOut(i);
		}
	}

	c.ls.Inout.savetxt("lsinout");

#ifdef INOUTNORMAL
	c.ls.InoutNorms.zeros(NUMINOUT);
	c.ls.InoutNorms(1).x = IN_NORM_X;
	c.ls.InoutNorms(1).y = IN_NORM_Y;
	c.ls.InoutNorms(1).z = IN_NORM_Z;
	c.ls.InoutNorms(2).x = OUT2_NORM_X;
	c.ls.InoutNorms(2).y = OUT2_NORM_Y;
	c.ls.InoutNorms(2).z = OUT2_NORM_Z;
	c.ls.InoutNorms(3).x = OUT3_NORM_X;
	c.ls.InoutNorms(3).y = OUT3_NORM_Y;
	c.ls.InoutNorms(3).z = OUT3_NORM_Z;
	c.ls.InoutNorms(4).x = OUT4_NORM_X;
	c.ls.InoutNorms(4).y = OUT4_NORM_Y;
	c.ls.InoutNorms(4).z = OUT4_NORM_Z;
	c.ls.InoutNorms(5).x = OUT5_NORM_X;
	c.ls.InoutNorms(5).y = OUT5_NORM_Y;
	c.ls.InoutNorms(5).z = OUT5_NORM_Z;
#endif // INOUTNORMAL
#endif // BIO_FLOW

#ifdef LS_DEBUG
	c.ls.BLc.savetxt("lsblc");
#endif //LS_DEBUG
#ifdef LS_RUPTURE
	c.ls.BL.save("lsbl");
#else //LS_RUPTURE
	c.ls.BL.savetxt("lsbl");
#endif //LS_RUPTURE
}
void clConstants::clLS::removeextraBL(void)
{
	countBL();
	removedoubleBL();
	removesingleBL();

	removenoBL();
	removeperiodicBL();

	//c.ls.BLc.savetxt("lsblc2");
	compressarrayBL();
	c.ls.BLFb.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);
}
void clConstants::clLS::testflatBL(int n0, int n1, int n2)
{
	clVector3Di BL0 = clVector3Di(n0, n1, n2).sort();
	for(int i0 = 0; i0 < c.ls.nBL; i0++)
	{
		clVector3Di BL1 = clVector3Di(c.ls.BL(i0, 0), c.ls.BL(i0, 1), c.ls.BL(i0, 2)).sort();
		int ident = BL0.compare(BL1);
		if(ident == 7)
			return;
		if(ident == (1 + 2))
			if(testoverlapBL(BL0.x, BL0.y, BL0, BL1) == TRUE)
				return;
		if(ident == (1 + 4))
			if(testoverlapBL(BL0.x, BL0.z, BL0, BL1) == TRUE)
				return;
		if(ident == (2 + 4))
			if(testoverlapBL(BL0.y, BL0.z, BL0, BL1) == TRUE)
				return;
	}
	if(c.ls.BL.getSize1() < c.ls.nBL + 2)
		c.ls.BL.reallocate(c.ls.nBL + ARRAY_INCREASE_SIZE, NUMBER_DIMENSIONS);

	c.ls.BL(c.ls.nBL, 0) = n0;
	c.ls.BL(c.ls.nBL, 1) = n1;
	c.ls.BL(c.ls.nBL, 2) = n2;
	c.ls.nBL++;
	c.ls.BL(c.ls.nBL, 0) = n0;
	c.ls.BL(c.ls.nBL, 1) = n2;
	c.ls.BL(c.ls.nBL, 2) = n1;
	c.ls.nBL++;
}

void clConstants::clLS::testBL(int ind, int i0, int i1, int i2)
{
	int n0 = vlb.M.Sf(ind,i0), n1 = vlb.M.Sf(ind,i1), n2 = vlb.M.Sf(ind,i2);
	int ind0 = c.ls.BNind(n0), ind1 = c.ls.BNind(n1), ind2 = c.ls.BNind(n2);
	if((ind0 >= 0 && ind1 >= 0 && ind2 >= 0 && c.ls.BS(ind0) == c.ls.BS(ind1) && c.ls.BS(ind0) == c.ls.BS(ind2))
				|| ((ind0 == -1 || ind1 == -1 || ind2 == -1)
				&& (c.ls.T(n0) & LS_OUTER_NODE && c.ls.T(n1) & LS_OUTER_NODE && c.ls.T(n2) & LS_OUTER_NODE)))
	{
		clVector3Dd NL = directionBL(n0) + directionBL(n1) + directionBL(n2); 
		clVector3Dd C0 = vls.C(n0), C1 = vls.C(n1), C2 = vls.C(n2);
		clVector3Dd Cn = (C1 - C0) ^ (C2 - C0);
		if(NL * Cn > 0.)
		{
			c.ls.BL(c.ls.nBL, 0) = n0;
			c.ls.BL(c.ls.nBL, 1) = n1;
			c.ls.BL(c.ls.nBL, 2) = n2;
		}
		else
		{
			c.ls.BL(c.ls.nBL, 0) = n0;
			c.ls.BL(c.ls.nBL, 1) = n2;
			c.ls.BL(c.ls.nBL, 2) = n1;
		}
		c.ls.nBL++;
	}
}

void clConstants::clLS::testInOut(int ind)
{
		clVector3Dd Vin, Vout;
		clVector3Dd C0 = vls.C(c.ls.BL(ind,0)), C1 = vls.C(c.ls.BL(ind,1)), C2 = vls.C(c.ls.BL(ind,2));
		clVector3Dd Cn = (C1 - C0) ^ (C2 - C0);
		Vin.x = IN_NORM_X;	Vin.y = IN_NORM_Y;	Vin.z = IN_NORM_Z;
		Vout.x = -OUT2_NORM_X;	Vout.y = -OUT2_NORM_Y;	Vout.z = -OUT2_NORM_Z;

		if((Vin - Cn.n()).r() < 0.01)
		{
			c.ls.Inout(ind) = 1;
		}
		else if((Vout - Cn.n()).r() < 0.01)
		{
			c.ls.Inout(ind) = 2;
		}
}

void clConstants::clLS::countBL(void)
{
	c.ls.BLc.zeros();
	for(int i0 = 0; i0 < c.ls.nN; i0++)
		for(int i1 = 0; i1 < c.ls.nL; i1++)
			c.ls.BLn(i0,i1).reset();

	for(int i0 = 0; i0 < c.ls.nBL; i0++)
	{
		clVector3Di BL0 = clVector3Di(c.ls.BL(i0, 0), c.ls.BL(i0, 1), c.ls.BL(i0, 2)).sort();
		int n0 = BL0.x, n01 = c.ls.N.find(BL0.x, BL0.y), n02 = c.ls.N.find(BL0.x, BL0.z);
		int	n1 = BL0.y, n11 = c.ls.N.find(BL0.y, BL0.z);
		c.ls.BLc(n0, n01)++;
		c.ls.BLc(n0, n02)++;
		c.ls.BLc(n1, n11)++;
		c.ls.BLn(n0, n01).add(i0);
		c.ls.BLn(n0, n02).add(i0);
		c.ls.BLn(n1, n11).add(i0);
	}
//	c.ls.BLc0.copy(c.ls.BLc);
}
void clConstants::clLS::removedoubleBL(void)
{
	//c.ls.BLc.savetxt("lsblc0");
	for(int i0 = 0; i0 < c.ls.nN; i0++)
	{
		for(int i1 = 0; i1 < c.ls.nL; i1++)
		{
			//int ntot = c.ls.BLc(i0,i1);
			int ntot = c.ls.BLn(i0,i1).nN;
			for(int i2 = 0; i2 < ntot - 1; i2++)
			{
				int nbl0 = c.ls.BLn(i0, i1)(i2);
				clVector3Di BL0 = clVector3Di(c.ls.BL(nbl0, 0), c.ls.BL(nbl0, 1), c.ls.BL(nbl0, 2)).sort();
				if(BL0.x == -1)
					continue;

				int ind0 = c.ls.BNind(BL0.x), ind1 = c.ls.BNind(BL0.y), ind2 = c.ls.BNind(BL0.z);
				if(ind0 < 0 || ind1 < 0 || ind2 < 0)
					continue;

				for(int i3 = i2 + 1; i3 < ntot; i3++)
				{
					int nbl1 = c.ls.BLn(i0, i1)(i3);
					clVector3Di BL1 = clVector3Di(c.ls.BL(nbl1, 0), c.ls.BL(nbl1, 1), c.ls.BL(nbl1, 2)).sort();
					if(BL1.x == -1)
						continue;

					if(testoverlapBL(i0, c.ls.N(i0,i1), BL0, BL1) == TRUE)
					{
						removecountBL(nbl0);
						removecountBL(nbl1);
						break;
					}
				}
			}
		}
	}
	//c.ls.BLc.savetxt("lsblc1");
}

BOOL clConstants::clLS::testoverlapBL(int n0, int n1, clVector3Di BL0, clVector3Di BL1)
{

	if(BL0 == BL1)
		return TRUE;

	//if(ntot == 2)
	//	continue; // check if this is correct

	int n02, n12;
	if(BL0.x != n0 && BL0.x != n1)
		n02 = BL0.x;
	else if(BL0.y != n0 && BL0.y != n1)
		n02 = BL0.y;
	else 
		n02 = BL0.z;

	if(BL1.x != n0 && BL1.x != n1)
		n12 = BL1.x;
	else if(BL1.y != n0 && BL1.y != n1)
		n12 = BL1.y;
	else 
		n12 = BL1.z;

	if(c.ls.N.find(n02, n12) == -1)
		return FALSE;

	clVector3Dd vR0 = vls.C(n0), vN0 = clVector3Dd(vls.C(n02) - vR0).n();
	clVector3Dd vR1 = vls.C(n1), vN1 = clVector3Dd(vls.C(n12) - vR1).n();

	double N01 = vN0 * vN1;
	double den = 1. - N01 * N01;

	if(fabs(den) > c.eps)
	{
		clVector3Dd vR01 = vR1 - vR0;
		double R01N0 = vR01 * vN0, R01N1 = vR01 * vN1;
		double L0 = (R01N0 - R01N1 * N01) / den;
		double L1 = -(R01N1 - R01N0 * N01) / den;
		double dist = clVector3Dd(vR01 + vN1 * L1 - vN0 * L0).abs();
		if(dist < 0.1)
			return TRUE;
	}
	else
	{
		vN0 = clVector3Dd(vls.C(n12) - vR0).n();
		vN1 = clVector3Dd(vls.C(n02) - vR1).n();

		N01 = vN0 * vN1;
		den = 1. - N01 * N01;

		if(fabs(den) > c.eps)
		{
			clVector3Dd vR01 = vR1 - vR0;
			double R01N0 = vR01 * vN0, R01N1 = vR01 * vN1;
			double L0 = (R01N0 - R01N1 * N01) / den;
			double L1 = -(R01N1 - R01N0 * N01) / den;
			double dist = clVector3Dd(vR01 + vN1 * L1 - vN0 * L0).abs();
			if(dist < 0.1)
				return TRUE;
		}
	}
	return FALSE;
}
void clConstants::clLS::compressarrayBL(void)
{
	int i0 = 0, i1 = c.ls.nBL - 1;

	for(; i0 < i1; i0++)
	{
		if(c.ls.BL(i0, 0) == -1)
		{
			for(; i1 > i0; i1--)
			{
				if(c.ls.BL(i1, 0) > -1)
				{
					c.ls.BL(i0,0) = c.ls.BL(i1,0);
					c.ls.BL(i0,1) = c.ls.BL(i1,1);
					c.ls.BL(i0,2) = c.ls.BL(i1,2);
					c.ls.BL(i1,0) = -1;
					i1 -= 1;
					break;
				}
			}
		}
	}
	for(; i1 >= 0; i1--)
		if(c.ls.BL(i1, 0) > -1)
			break;

	c.ls.nBL = i1 + 1;
	if(c.ls.nBL < 0) c.ls.nBL = 0;
	c.ls.BL.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);
}
void clConstants::clLS::removesingleBL(void)
{
	for(int i0 = 0; i0 < c.ls.nN; i0++)
		for(int i1 = 0; i1 < c.ls.nL; i1++)
			if(c.ls.BLc(i0,i1) == 1)
			{
				//if(c.ls.NPN(i0) != -1 && c.ls.NPN(c.ls.N(i0,i1)) != -1)
				//{
				//	continue;
				//}
				int ntot = c.ls.BLn(i0,i1).nN;
				//for(int i2 = 0; i2 < c.ls.BLc0(i0,i1); i2++)
				for(int i2 = 0; i2 < ntot; i2++)
				{
					int nbl = c.ls.BLn(i0, i1)(i2);
					clVector3Di bl = clVector3Di(c.ls.BL(nbl, 0), c.ls.BL(nbl, 1), c.ls.BL(nbl, 2)).sort();
					if(bl.x == -1)
						continue;
					int ind0 = c.ls.BNind(bl.x), ind1 = c.ls.BNind(bl.y), ind2 = c.ls.BNind(bl.z);
					//int n0 = bl.x, n01 = c.ls.N.find(bl.x, bl.y), n02 = c.ls.N.find(bl.x, bl.z);
					//int	n1 = bl.y, n11 = c.ls.N.find(bl.y, bl.z);
//					if(c.ls.BLc(n0, n01) > 2 || c.ls.BLc(n0, n02) > 2 || c.ls.BLc(n1, n11) > 2)
					if(ind0 >= 0 && ind1 >= 0 && ind2 >= 0)
						removecountBL(nbl);
//					c.ls.BL(nbl, 0) = -1;

				}
			}
}
void clConstants::clLS::removeperiodicBL(void)
{
	for(int i0 = 0; i0 < c.ls.nBL; i0++)
	{
		int n0 = c.ls.BL(i0,0), n1 = c.ls.BL(i0,1), n2 = c.ls.BL(i0,2); 
		if(n0 == -1 || n1 == -1 || n2 == -1 || (c.ls.NPN(n0) != -1 && c.ls.NPN(n1) != -1 && c.ls.NPN(n2) != -1))
			removecountBL(i0);
	}
}
void clConstants::clLS::removenoBL(void)
{
	for(int i0 = 0; i0 < c.ls.nBL; i0++)
	{
		int n0 = c.ls.BL(i0,0), n1 = c.ls.BL(i0,1), n2 = c.ls.BL(i0,2); 
		if(n0 == -1 || n1 == -1 || n2 == -1 || (c.ls.BNind(n0) == -1 || c.ls.BNind(n1) == -1 || c.ls.BNind(n2) == -1))
			removecountBL(i0);
	}
}
BOOL clConstants::clLS::removecountBL(int nbl)
{
	clVector3Di bl = clVector3Di(c.ls.BL(nbl, 0), c.ls.BL(nbl, 1), c.ls.BL(nbl, 2)).sort();
	if(bl.x == -1)
		return FALSE;
	int n0 = bl.x, n01 = c.ls.N.find(bl.x, bl.y), n02 = c.ls.N.find(bl.x, bl.z);
	int	n1 = bl.y, n11 = c.ls.N.find(bl.y, bl.z);
	c.ls.BLc(n0, n01)--;
	c.ls.BLc(n0, n02)--;
	c.ls.BLc(n1, n11)--;
	c.ls.BL(nbl, 0) = -1;
	return TRUE;
}
//void clConstants::clLS::testdoubleBL(void)
//{
//
//	for(int i0 = 0; i0 < c.ls.nBL; i0++)
//	{
//		clVector3Di BL0 = clVector3Di(c.ls.BL(i0, 0), c.ls.BL(i0, 1), c.ls.BL(i0, 2)).sort();
//		for(int i1 = i0 + 1; i1 < c.ls.nBL; i1++)
//		{
//			clVector3Di BL1 = clVector3Di(c.ls.BL(i1, 0), c.ls.BL(i1, 1), c.ls.BL(i1, 2)).sort();
//			if(BL0 == BL1)
//			{
//				for(int j = i0+1; j < i1; j++)
//				{
//					c.ls.BL(j-1, 0) = c.ls.BL(j, 0);
//					c.ls.BL(j-1, 1) = c.ls.BL(j, 1);
//					c.ls.BL(j-1, 2) = c.ls.BL(j, 2);
//				}
//				for(int j = i1+1; j < c.ls.nBL; j++)
//				{
//					c.ls.BL(j-2, 0) = c.ls.BL(j, 0);
//					c.ls.BL(j-2, 1) = c.ls.BL(j, 1);
//					c.ls.BL(j-2, 2) = c.ls.BL(j, 2);
//				}
//				c.ls.nBL -= 2;
//				i0 -= 1;
//				i1 -= 2;
//				if(i0 >= i1)
//					i1 = i0 + 1;
//			}
//		}
//	}
//}
// old version allows overlaping triangles
//void clConstants::clLS::testdoubleBL(void)
//{
//	for(int i0 = 0; i0 < c.ls.nBL; i0++)
//	{
//		clVector3Di BL0 = clVector3Di(c.ls.BL(i0, 0), c.ls.BL(i0, 1), c.ls.BL(i0, 2)).sort();
//		for(int i1 = i0 + 1; i1 < c.ls.nBL; i1++)
//		{
//			clVector3Di BL1 = clVector3Di(c.ls.BL(i1, 0), c.ls.BL(i1, 1), c.ls.BL(i1, 2)).sort();
//			if(BL0 == BL1)
//			{
//				for(int j = i0+1; j < i1; j++)
//				{
//					c.ls.BL(j-1, 0) = c.ls.BL(j, 0);
//					c.ls.BL(j-1, 1) = c.ls.BL(j, 1);
//					c.ls.BL(j-1, 2) = c.ls.BL(j, 2);
//				}
//				for(int j = i1+1; j < c.ls.nBL; j++)
//				{
//					c.ls.BL(j-2, 0) = c.ls.BL(j, 0);
//					c.ls.BL(j-2, 1) = c.ls.BL(j, 1);
//					c.ls.BL(j-2, 2) = c.ls.BL(j, 2);
//				}
//				c.ls.nBL -= 2;
//				i0 -= 1;
//				i1 -= 2;
//				if(i0 >= i1)
//					i1 = i0 + 1;
//			}
//		}
//	}
//}
//void clConstants::clLS::removedoubleBL(void)
//{
//	//for(int i0 = 0; i0 < c.ls.nBL; i0++)
//	//{
//	//	clVector3Di BL0 = clVector3Di(c.ls.BL(i0, 0), c.ls.BL(i0, 1), c.ls.BL(i0, 2)).sort();
//	//	for(int i1 = i0 + 1; i1 < c.ls.nBL; i1++)
//	//	{
//	//		int n0 = c.ls.BL(i1, 0), n1 = c.ls.BL(i1, 1), n2 = c.ls.BL(i1, 2);
//	//		clVector3Di BL1 = clVector3Di(n0, n1, n2).sort();
//	//		if(BL0 == BL1)
//	//		{
//	//			//int nv0, nv1, nv2;
//	//			//if(c.ls.N(n0,6) >= 0)
//	//			//	nv0 = c.ls.N(n0,6);
//	//			//else
//	//			//	nv0 = c.ls.N(n0,13);
//
//	//			//if(c.ls.N(n1,6) >= 0)
//	//			//	nv1 = c.ls.N(n1,6);
//	//			//else
//	//			//	nv1 = c.ls.N(n1,13);
//
//	//			//if(c.ls.N(n2,6) >= 0)
//	//			//	nv2 = c.ls.N(n2,6);
//	//			//else
//	//			//	nv2 = c.ls.N(n2,13);
//
//	//			//c.ls.N.remove_link(n0, n1);
//	//			//c.ls.N.remove_link(n0, n2);
//	//			//c.ls.N.remove_link(n1, n2);
//	//			//c.ls.N.remove_link(n0, nv1);
//	//			//c.ls.N.remove_link(n0, nv2);
//	//			//c.ls.N.remove_link(n1, nv0);
//	//			//c.ls.N.remove_link(n1, nv2);
//	//			//c.ls.N.remove_link(n2, nv0);
//	//			//c.ls.N.remove_link(n2, nv1);
//
//	//			for(int j = i0+1; j < i1; j++)
//	//			{
//	//				c.ls.BL(j-1, 0) = c.ls.BL(j, 0);
//	//				c.ls.BL(j-1, 1) = c.ls.BL(j, 1);
//	//				c.ls.BL(j-1, 2) = c.ls.BL(j, 2);
//	//			}
//	//			for(int j = i1+1; j < c.ls.nBL; j++)
//	//			{
//	//				c.ls.BL(j-2, 0) = c.ls.BL(j, 0);
//	//				c.ls.BL(j-2, 1) = c.ls.BL(j, 1);
//	//				c.ls.BL(j-2, 2) = c.ls.BL(j, 2);
//	//			}
//	//			c.ls.nBL -= 2;
//	//			i0 -= 1;
//	//			i1 -= 2;
//	//			if(i0 >= i1)
//	//				i1 = i0 + 1;
//	//		}
//	//	}
//	//}
//	testdoubleBL();
//
//	if(c.ls.nBL != c.ls.BL.getSize1())
//		c.ls.BL.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);
//}
clVector3Dd clConstants::clLS::directionBL(int n)
{
	clVector3Dd C0 = vls.C(n); 
	clVector3Dd Ctot; 
	for(int j = 0; j < c.ls.nL; j++)
		if(c.ls.N(n,j) >= 0)
		{
			clVector3Dd dC = (vls.C(c.ls.N(n,j)) - C0); 
			Ctot += dC;
		}

	return -Ctot;
}
void clConstants::clLS::removeBL(int n1, int n2)
{
	for(int i = 0; i < c.ls.nBL; i++)
	{
		int bl0 = c.ls.BL(i, 0), bl1 = c.ls.BL(i, 1), bl2 = c.ls.BL(i, 2);
		if((bl0 == n1 && bl1 == n2) || (bl0 == n2 && bl1 == n1) || 
				(bl0 == n1 && bl2 == n2) || (bl0 == n2 && bl2 == n1) ||
				(bl1 == n1 && bl2 == n2) || (bl1 == n2 && bl2 == n1))
			removeBL(i--);
	}
}
void clConstants::clLS::removeBL(int n)
{
	for(int i = n+1; i < c.ls.nBL; i++)
	{
		c.ls.BL(i-1, 0) = c.ls.BL(i, 0);
		c.ls.BL(i-1, 1) = c.ls.BL(i, 1);
		c.ls.BL(i-1, 2) = c.ls.BL(i, 2);
	}
	c.ls.nBL--;
	//c.ls.BL.reallocate(c.ls.nBL, NUMBER_DIMENSIONS);
}


void clConstants::clLS::addBL(int n0, int n1, int n2)
{
	c.ls.nBL++;
	if(c.ls.nBL > c.ls.BL.getSize1())
		c.ls.BL.reallocate(c.ls.nBL + LS_ARRAY_INCREASE_STEP, NUMBER_DIMENSIONS);
	c.ls.BL(c.ls.nBL-1, 0) = n0;
	c.ls.BL(c.ls.nBL-1, 1) = n1;
	c.ls.BL(c.ls.nBL-1, 2) = n2;
}
void clConstants::clLS::iniT()
{
	//c.ls.nTp = 0;
	c.ls.nTa = 0;
	c.ls.nTr = 0;
	c.ls.nTs = 0;

	c.ls.T.zeros(c.ls.nN);
	c.ls.Ts.zeros(c.ls.nN);
	c.ls.Tno_c.zeros(c.ls.nN);
	c.ls.Tno_f.zeros(c.ls.nN);
	c.ls.T.loadtxt("load" SLASH "lst");

	c.ls.TaDe.zeros(c.ls.nN);
	if(c.ls.TaDe.loadtxt("load" SLASH "lstade") == FALSE)
	{
		for(int i = 0; i < c.ls.nN; i++)
			{
				if(c.ls.O(i) == 0)
					c.ls.TaDe(i) = c.ls.FwDe0;
				if(c.ls.O(i) == 1)
					c.ls.TaDe(i) = c.ls.FwDe1;
				//Modified
				if(c.ls.O(i) >= 2)
				//Modified
					c.ls.TaDe(i) = c.ls.FwDe2;
			}
	}

	for(int i = 0; i < c.ls.nN; i++)
	{
		if(c.ls.RO(c.ls.O(i)) == TRUE && c.ls.T(i) & LS_STATIC_NODE)
			c.ls.SO(c.ls.O(i)) = TRUE;
	}

	for(int i = 0; i < c.ls.nN; i++)
	{
		if(c.ls.RO(c.ls.O(i)) == TRUE)
			c.ls.T(i) = c.ls.T(i) | LS_RIGID_NODE;

		if(c.ls.SO(c.ls.O(i)) == TRUE)
			c.ls.T(i) = c.ls.T(i) | LS_STATIC_NODE;

		//if(c.ls.T(i) & LS_PASSIVE_NODE)
		//	c.ls.nTp++;

		//if((c.ls.T(i) & LS_ACTIVE_NODE) || (c.ls.T(i) & LS_ACTIVE_REP_NODE))
		if(c.ls.T(i) & LS_ACTIVE_NODE)
			c.ls.nTa++;

		if(c.ls.T(i) & LS_STATIC_NODE)
			c.ls.Ts(i) = TRUE;

		//if(c.ls.T(i) & LS_STATIC_NODE)
		//	c.ls.Tno_f(i) = TRUE;

		if(c.ls.T(i) & LS_STATIC_NODE || c.ls.T(i) & LS_PERIODIC_NODE)
			c.ls.Tno_c(i) = TRUE;

		if(c.ls.N.NL(i) == 0)
			c.ls.Tno_c(i) = TRUE;
	}

	for(int i = 0; i < c.ls.nN; i++)
	{
		if(c.ls.T(i) & LS_RIGID_NODE)
		{
//			c.ls.RO(c.ls.O(i)) = TRUE;
			c.ls.nTr++;
		}
		if(c.ls.Ts(i) == TRUE)
		{
			c.ls.nTs++;
		}
	}

	for(int i = 0; i < c.ls.nN; i++)
	{
		BOOL flS = TRUE;
		for(int j = 0; j < c.ls.nL; j++)
		{
			if(c.ls.N(i, j) >= 0 && !(c.ls.T(c.ls.N(i, j)) & LS_STATIC_NODE))
			{
				flS = FALSE;
				break;
			}
		}
		//c.ls.Tno_f(i) = flS;
		////if(c.ls.T(i) & LS_STATIC_NODE)
		////	c.ls.Tno_f(i) = TRUE;
		////else
		////	c.ls.Tno_f(i) = FALSE;
	}
//	c.ls.Ta.zeros(c.ls.nTa,c.lb.nD);
	c.ls.Ta.zeros(c.ls.nTa);
//	c.ls.Tp.zeros(c.ls.nTp);

	c.ls.Tr.zeros(c.ls.nTr);

	int ip = 0, ia = 0, ir = 0;
	for(int i = 0; i < c.ls.nN; i++)
	{
		//if(c.ls.T(i) & LS_PASSIVE_NODE)
		//	c.ls.Tp(ip++) = i;

		//if((c.ls.T(i) & LS_ACTIVE_NODE) || (c.ls.T(i) & LS_ACTIVE_REP_NODE))
		if(c.ls.T(i) & LS_ACTIVE_NODE)
			c.ls.Ta(ia++) = i;

		if(c.ls.T(i) & LS_RIGID_NODE)
			c.ls.Tr(ir++) = i;
	}

	c.ls.T.savetxt("lst");
	c.ls.Ts.savetxt("lsts");
	//c.ls.Tp.savetxt("lstp");
	c.ls.Ta.savetxt("lsta");
	c.ls.TaDe.savetxt("lstade");
	c.ls.Tr.savetxt("lstr");
}

void clConstants::clLS::iniP()
{
	c.ls.NPN.zeros(c.ls.nN);
	c.ls.NPN.fill(-1);
	if(c.ls.PN.loadtxtnew("load" SLASH "lspn", 2) == FALSE)
		return;

	c.ls.nPN = c.ls.PN.getSize(1);
	c.ls.PNC.zeros(c.ls.nPN);

	if(c.ls.PNC.loadtxt("load" SLASH "lspnc") == FALSE)
	{
		for(int i = 0; i < c.ls.nPN; i++)
			c.ls.PNC(i) = vls.C.c0(c.ls.PN(i,0)) -  vls.C.c0(c.ls.PN(i,1));
	}

	for(int i = 0; i < c.ls.nPN; i++)
	{
		c.ls.NPN(c.ls.PN(i,0)) = c.ls.PN(i,1);
		c.ls.NPN(c.ls.PN(i,1)) = c.ls.PN(i,0);
	}


	c.ls.PN.savetxt("lspn");
	c.ls.PNC.savetxt("lspnc");
}
void clConstants::clLS::updateO()
{
	nOp = 0;
	Op.fill(-1);
	for(int i = 0; i < c.ls.nN; i++)
	{
		int nStack = 0, nCur = -1;
		int n0 = i;
		if(Op(n0) >= 0)
			continue;

		while(TRUE)
		{
			if(Op(n0) < 0)
			{
				Op(n0) = nOp;
				for(int j = 0; j < c.ls.nL; j++)
				{
					int n1 = N(n0,j);
					if(n1 >= 0 && Op(n1) < 0)
						Opstack(nStack++) = n1;
				}
			}
			if(++nCur >= nStack)
				break;
			n0 = Opstack(nCur);
		}
		nOp++;
	}
}
//void clConstants::clLS::updateO()
//{
//	nOp = 0;
//	Op.fill(-1);
//	for(int i = 0; i < c.ls.nN; i++)
//		if(Op(i) < 0)
//			setOp(i, nOp++);
//}
//void clConstants::clLS::setOp(int n, int p)
//{
//	if(Op(n) == p)
//		return;
//
//	if(Op(n) >= 0)
//	{
//		printf("\nError in Op, double index");
//		return;
//	}
//
//	Op(n) = p;
//	for(int i = 0; i < c.ls.nL; i++)
//		if(N(n,i) >= 0)
//			setOp(N(n,i), p);
//}
void clConstants::clLS::iniO()
{
	c.ls.O.zeros(c.ls.nN);
	c.ls.O.loadtxt("load" SLASH "lso");

	c.ls.Os.zeros(c.ls.nN);
	c.ls.Os.loadtxt("load" SLASH "lsos");

	c.ls.Op.zeros(c.ls.nN);
	c.ls.Opstack.zeros(c.ls.nN * c.ls.nL);

	nO = 0;
	for(int i = 0; i < c.ls.nN; i++)
		if(c.ls.O(i) > c.ls.nO)
			c.ls.nO = c.ls.O(i);

	if(c.ls.nN > 0)
		c.ls.nO++;

	c.ls.NO.zeros(c.ls.nO);
	for(int i = 0; i < c.ls.nN; i++)
		c.ls.NO(c.ls.O(i))++;

	updateO();

	c.ls.O.savetxt("lso");
	c.ls.Op.savetxt("lsop");
	c.ls.NO.savetxt("lsno");
//////////////////////////////////////////////////////////////////////
	c.ls.RO.zeros(c.ls.nO);
	if(c.ls.RO.loadtxt("load" SLASH "lsro") == FALSE)
	{
		c.ls.RO.fill(FALSE);
	}
	//c.ls.RO(0) = TRUE;
	//c.ls.RO(1) = TRUE;
	c.ls.RO.savetxt("lsro");
//////////////////////////////////////////////////////////////////////
	c.ls.SO.zeros(c.ls.nO);
	if(c.ls.SO.loadtxt("load" SLASH "lsso") == FALSE)
	{
		c.ls.SO.fill(FALSE);
	}
//	c.ls.SO(0) = TRUE;
	c.ls.SO.savetxt("lsso");
//////////////////////////////////////////////////////////////////////
	c.ls.PO.zeros(c.ls.nO);
	if(c.ls.PO.loadtxt("load" SLASH "lspo") == FALSE)
	{
		c.ls.PO.fill(clVector3Dc(TRUE, TRUE, TRUE));
	}
//	c.ls.PO(0) = clVector3Dc(FALSE, FALSE, FALSE);
	c.ls.PO.savetxt("lspo");
//////////////////////////////////////////////////////////////////////
	c.ls.DO.zeros(c.ls.nO);
	if(c.ls.DO.loadtxt("load" SLASH "lsdo") == FALSE)
	{
		c.ls.DO.fill(TRUE);
	}
//	c.ls.DO(0) = TRUE;
	c.ls.DO.savetxt("lsdo");
}

