// ConstDef.h : include file for constants
// (c) Alexander Alexeev, 2007 
//

#if !defined(__CONSTDEF_H__INCLUDED_)
#define __CONSTDEF_H__INCLUDED_

//#define NUM_THREADS 0

#define LB_USE_OPENMP
#define LS_USE_OPENMP
#define BC_USE_OPENMP			// parallel boundary-link search (clMethod.h bcFindBoundaryNodes)
#define WALL_LAW_PAPER			// r = R - A/2 (1-cos(2 pi x/L)) sin(wt) between valve roots (clMethod.h wallMod)
#define LS_HINGE_BENDING
#define LS_HINGE_WALL_FACTOR	0.2		// hinge stiffness x0.2 at the vessel wall ...
#define LS_HINGE_WALL_RAMP	3.0		// ... rising linearly to full over 3 lattice units
#define LS_STRETCH_EDGE_FACTOR	0.1		// leaflet in-plane springs x0.1 at the free edge ...
#define LS_STRETCH_EDGE_RAMP	24.0	// ... rising linearly to full over 24 lattice units (whole leaflet)
#define LS_COMMISSURE_TIE_K	5.0		// tie outermost free-edge nodes of opposite leaflets (2026-09-26)
#define LS_COMMISSURE_TIE_D	5.0
#define LS_COMMISSURE_TIE_OPEN	0.5		// tie only while the valve centre is shut (2026-09-26)
#define TR_USE_OPENMP


#define LB_SOLVER
#define LS_SOLVER
//#define TR_SOLVER
#define LB_USE_LBGK
//#define LB_USE_MRTLB

#define TR_USE_DIFFUSION
//#define TR_SOLVER_2D
//#define TR_ABSORBING_CILIA

//#define CREATE_BIN_FILES
#define LB_TEST_BC
#define LB_THIN_WALL
//#define LB_COARSE
#define LS_DISSIPATION
//#define LS_POTENCIAL
//#define LS_RUPTURE
//#define RENAME_BIN_FILES
#define RENAME_DUMP_FILES
//#define SAVE_BIN_IN_LOAD
//#define LS_FAST_FORCE_CALC
#define LB_NO_FLUX_CORRECTION
#define LB_NO_MASS_CORRECTION
//#define LBN_CORRECTION
//#define BINARY_FLUID
//#define BINARY_SCALAR
//#define SET_ORBIT		// Define to prescribe bead orbit
#ifndef SET_ORBIT
//#define LS_BEADS		// Define if using beads with magnetic field
#endif // SET_ORBIT
//#define LS_BEADS_FRICTION
//#define BEADS_MIXING
//#define STATIC_LS
#define LEAFLETS
//#define STATIC_VESSEL
//#define TRANSLATE_LS		(-0.01)
//#define MIXSTOP
//#define FIND_JUNCTIONS	// Comment this out if there is no vessel or more than one vessel object
#define LB_HW_BB				// half-way bounce-back: the interpolated (Bouzidi) scheme leaked ~35% of the stroke through the wall, mostly at the leaflet-wall junctions (2026-09-24)
#define BIO_FLOW
//#define SYMMETRIC_VESSEL
#define CONSTRAIN_MIDLINE
#define SYMMETRIC_GAP_DIST	(0.02)
#define NUM_VESSELS			(1)	// Number of vessel objects, for determining which nodes to constrain to not cross the symmetry plane (flaps only)

//#define OSC_FLOW
#define SQUARE_WAVE
#define SQ_BETA					(0.015)	// Beta determines roundedness of square wave
//#define INVEL_BC
#define INLET_P_BC
#define OUTLET_P_BC
//#define OUTLET_INTERP
//#define RAMP
#define OSCILLATING_VESSEL
#define OSC_AMPLITUDE 3
//#define REACTIVE_PRESS
#define DOMAIN_FACT				(2.5)
#define INTERVAL				(150)
//Modified
#define SAMPX					(20) //X-coordinate of the sample to measure Q from
//Modified

#define INLET_VEL_MAG			(0.001)
#define CYCLE_PERIOD 188496
#define RAMP_PER_CYCLE			(int)(2)
#define	RAMP_PERIOD				(int)(CYCLE_PERIOD/RAMP_PER_CYCLE)
#define PHASE_DIFFERENCE 0
#define INHALE					(int)(CYCLE_PERIOD / 2.)	// Sine function input for Womersley model
#define EXHALE					(int)(CYCLE_PERIOD - INHALE)

#define OSC_PERIOD				(int)(2.*INHALE)
#define OSC_PERIOD_2			(int)(2.*EXHALE)
#ifdef INVEL_BC
#define OSC_AMPL				(INLET_VEL_MAG)
#define OSC_AMPL_2				(INLET_VEL_MAG / 2.)
#else
#define OSC_AMPL				(1.)
#define OSC_AMPL_2				(1.)
#endif // INVEL_BC
#define DELTA_RO 0
#define DELTA_RO_WALL			DELTA_RO	// Density difference divided by six (direction gives resulting flow direction) - Zhou He

#define INOUTNORMAL				// Forces flow at inlet/outlets to be normal to inlet/outlet (INLET/OUTLET NORMALS MUST BE CONSTANT THROUGHOUT RUN AS WRITTEN)
#define NUMINOUT				(6)		// Total number of inlets and outlets PLUS ONE
#define IN_NORM_X				(1.)	// Unit normal INTO VESSEL for inout 1 (inlet)
#define IN_NORM_Y				(0.)
#define IN_NORM_Z				(0.)
#define OUT2_NORM_X				(1.)	// Unit normals OUT OF VESSEL for inout 2-5 (outlets)
#define OUT2_NORM_Y				(0.)
#define OUT2_NORM_Z				(0.)
#define OUT3_NORM_X				(0.)
#define OUT3_NORM_Y				(0.)
#define OUT3_NORM_Z				(0.)
#define OUT4_NORM_X				(0.)
#define OUT4_NORM_Y				(0.)
#define OUT4_NORM_Z				(0.)
#define OUT5_NORM_X				(0.)
#define OUT5_NORM_Y				(0.)
#define OUT5_NORM_Z				(0.)


#define TRC_NUM_TRACERS		(int)(1e6)
#define PECLET_NUMBER_BEAD	(500.)
#define TR_D0				(2 * R0_NUMBER * BEAD_VELOCITY / PECLET_NUMBER_BEAD)  //Diffusion coefficient for tracers
//#define TR_D0				0.0005  //Diffusion coefficient for tracers
#define BEADS_UPDATE_TIME	100 //time between updating list of bead position



//Set these defines to match the geometry used in simulation
//#define TR_RIDGES
#define TR_HIDDEN_DISKS   //If pillars are below wall keep this defined (method is only working right now with this defined)
#define TR_BEADS

///////Adsorbing vs Reflective Boundaries (1 for adsorbing 0 for reflective)
#define	BEAD_SURF					0
#define	PILLAR_SURF					0
#define RIDGE_SURF					0
#define	Z_SURF						0	
#define Y_SURF						0



//#define INI_SHEAR
//#define INI_POISEUILLE
//#define THETA_NUMBER				(800.)		//1000~1e4
#define THETA_NUMBER				(12.8)		//1000~1e4  -  FOR MIXING
#define LAMBDA_NUMBER				(25e5)		//
#define PHI_NUMBER					(0.)		
#define PILLAR_H1					(-0.5)
#define PILLAR_H2					(0.1)
#define PILLAR_RADIUS				(3.5)
#define N_DISK					4		// Set to zero if not using beads
#define N_BEAD					8		// Set to zero if not using beads
#define N_ROWS_DISK				4		// Set to zero if not using beads
#define N_BEADS_DISK			2		// Number of beads per disk
#define CX0_DISK					(10.)
#define CY0_DISK					(17.5)
#define SPACING					(21.)

#define R0_NUMBER				(3.5)		//bead radius
#define N0_NUMBER				(324)

#define ORBIT_FACTOR			(2.6 / 1.4)
#define ORBIT_RADIUS			(R0_NUMBER * ORBIT_FACTOR)
#define BEAD_VELOCITY			(0.0125)
#define BEAD_OMEGA				(BEAD_VELOCITY / ORBIT_RADIUS)
#define BEAD_OMEGA_SQ			(BEAD_OMEGA * BEAD_OMEGA)

//#ifdef SET_ORBIT
//#define OSC_PERIOD					(4000)//2 * PI_NUMBER / BEAD_OMEGA)
//#else
//#define OSC_PERIOD					((2. * R0_NUMBER) * (2. * R0_NUMBER) * THETA_NUMBER / MU_NUMBER)
//#endif // SET_ORBIT
#define MAGNETIC_FORCE				(LAMBDA_NUMBER * 3. * PI_NUMBER * MU_NUMBER * MU_NUMBER / THETA_NUMBER)
//#define B_BIAS					(0.18/MU0_MS)	//ratio of external magnetic field to mu0*Ms - FOR TRAJECTORY RUNS
#define B_BIAS					(0.1)	//ratio of external magnetic field to mu0*Ms - FOR MIXING RUNS
//#define MU0_MS					(0.83)	//tesla
#define MU0_MS					(100*0.83)	//tesla  -  FOR MIXING
#define MSP						(15e3)	//A/m, saturation magnetization of bead 
#define FORCE_CONVERSION		((1./4. * 1e-6)*(1./4. * 1e-6) / 9e-9)	//from SI(A*T/m) to LB(A * T/m)
//#define OSC_AMPL				(1.)
//#define PILLAR2CAPSULE			(12.* PI_NUMBER * PI_NUMBER * MU_NUMBER * MU_NUMBER / THETA_NUMBER)

#define RATIO					(0.)	//repulsive force ratio
#define TANG_FORCE				0.0		// Constant tangential force adjustment amplitude
#define PI_NUMBER				3.1415926535897931

// â”€â”€ IN-PLANE (EXTENSIONAL) STIFFNESS - TEMPLATE v2 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
// These are now the ONLY place the solid's in-plane stiffness is set.
// load/lsk.txt carries topology weights in [0,1] (1 = full-strength link,
// 0.5 = shared edge, 0 = no link) and clConstants.cpp multiplies each node's
// weight by the constant for its group. Written by editConstDef.m.
//
// GROUPS are the values in load/lso.txt, and the mapping in clConstants.cpp
// (clK::ini) is NOT one-per-leaflet:
//     _0  ->  O == 0   the vessel wall            (2640 nodes)
//     _1  ->  O == 1   inlet leaflet 1 ONLY       (183 nodes)
//     _2  ->  O >= 2   inlet leaflet 2 AND both outlet leaflets (549 nodes)
// So _1 and _2 must be kept EQUAL unless you deliberately want one half of
// the inlet valve stiffer than the other half. buildConfig.m writes them
// equal and warns if they diverge.
#define LS_STIFFNESS_0 2.598076211
#define LS_STIFFNESS_1 0.2501851166	// x0.1 of 2.501851166 (2026-09-26)
#define LS_STIFFNESS_2 0.2501851166	// x0.1 of 2.501851166 (2026-09-26)

#define LS_C_INI_X_0		0.
#define LS_C_INI_Y_0		0.
#define LS_C_INI_Z_0		0.

#define LS_DX_0					(4./3.)
#define LS_DX_1					(4./3.)
#define LS_DX_2					(4./3.)

#define LS_CRIT_F1				0.

#define LS_INI_VELOCITY			0.
#define LS_INI_FORCE			0.
#define LS_DAMPING_CONSTANT		0.1
////////////////////////////////////////////////////////////////////////////////////////
//#define SP_NUMBER				5.0
//#define DELTA_NUMBER			1.
////////////////////////////////////////////////////////////////////////////////////////
//#define MAX_TIME				(5.*100000.)
//#define MAX_PERIODS				(30)
//#define STOP_PERIODS			((int)(MAX_TIME / (double)OSC_PERIOD) < MAX_PERIODS ? (int)(MAX_TIME / (double)OSC_PERIOD) : MAX_PERIODS)
//#define STOP_TIME				(STOP_PERIODS * OSC_PERIOD)
//#define STOP_TIME				(100.)

#define STOP_TIME				(int)(1.*CYCLE_PERIOD)
#define DUMP_TOTAL				(10)
//#define SAVE_TOTAL				(1000)
#define DUMP_PER_PERIOD 200
#define PERIOD_PER_DUMP			(1)
#ifdef DUMP_PER_PERIOD
#define DUMP_STEP_SAVE			(CYCLE_PERIOD/DUMP_PER_PERIOD)
#else
#define DUMP_STEP_SAVE			(CYCLE_PERIOD*PERIOD_PER_DUMP)
#endif // DUMP_PER_PERIOD
#define DUMP_STEP				(DUMP_STEP_SAVE)	//(STOP_TIME/DUMP_TOTAL)
#define DUMP_START_TIME			(STOP_TIME - 1.*CYCLE_PERIOD)
//#define DUMP_START_TIME         (1)

#define CURRENT_SAVE_START_TIME	(1)
#define CURRENT_SAVE_STEP		(STOP_TIME < 5000 ? 1: (int)((STOP_TIME) / 5000))
//#define CURRENT_SAVE_STEP        (LS_STEPS)

//#define PARTICLE_RELEASE_PERIOD	(10.)
//#define PARTICLE_RELEASE_TIME	(PARTICLE_RELEASE_PERIOD*OSC_PERIOD)
#define PARTICLE_RELEASE_TIME	(0.)

#define LB_STEPS				(1)
#define LS_STEPS				(4)
#define TR_STEPS				(1)


// Binary Fluid Config
// LBB_A, LBB_B, and LBB_KAPPA determine interfacial tension
#define LBB_A					-0.001		// negative for homogenous state or positive for two phases
#define LBB_B					0.		// must be positive
#define LBB_KAPPA				0.00		// Interface Width for Species A
#define LBB_MOBILITY			0.08333		// Mobility Number
#define LBB_GAMMA_NUMBER		2.*LBB_MOBILITY		// For miscible fluids, D = -LBB_GAMMA_NUMBER*LBB_A/2 if LBB_B=0
#define LBB_SIGMA				0.00083		// surface tension


//#define LBB_LAMBDA_NUMBER		0.014		// Interaction Strength
//#define LBB_THETA_NUMBER		0.11	    // Reduced Temperature


// Changed 2026-08-13 from 0.6791635268 to 1.0.
//
// tau = 1.0 gives nu = cs^2 (tau - 1/2) = 1/6, which is the viscosity the
// paper states and the one buildConfig.m has always assumed when deriving the
// period and the valve stiffness. At 0.6791635268 the solver was actually
// running at nu = 0.0597 - a factor of 2.79 - so every case came out at
//   Gamma_actual = 1.67 * Gamma_requested   and   K_actual = 2.79 * K_requested
// i.e. the valve leaflets were ~2.8x STIFFER than intended, which is the
// leading explanation for the measured valve opening being 3.5-5x smaller
// than the paper's Figure 3a. See SESSION_NOTES.md, 2026-08-13.
//
// tau = 1.0 is also the numerically best-behaved BGK relaxation time.
#define LB_TAU_NUMBER			(1.0)
#define LB_TAUB_NUMBER			(1.)		// Relaxation for Species B
#ifndef BINARY_SCALAR
	#define LB_TAUN_NUMBER			(1.)		// Relaxation for Species A, 1/tau_N
#endif // BINARY_SCALAR

#ifdef BINARY_SCALAR
#define BINARY_DIFFUSION		(0.0005)
#define LB_TAUN_NUMBER			(1./(3*BINARY_DIFFUSION + 0.5))		// Relaxation for Species A, 1/tau_N
#endif // BINARY_SCALAR


#define WIDTH_NODES				4.
#define ASPECT_RATIO			10.
#define H0_NUMBER				((WIDTH_NODES - 1.) * LS_DX_1)
#define H_NUMBER				(H0_NUMBER + 0. * LS_DX_1 * 0.15)
#define L_NUMBER				(ASPECT_RATIO * H0_NUMBER)
//#define ANGLE_NUMBER			(1. * PI_NUMBER / 4.)
#define NUMBER_ACTIVE_NODES		(WIDTH_NODES * WIDTH_NODES)
#define CHANNEL_WIDTH			(H0_NUMBER * 2.)
//#define CHANNEL_HIGHT			(R0_NUMBER * 4.)
#define CHANNEL_HIGHT			(26.)
//#define R0_NUMBER				10.


#define EI_NUMBER				(LS_STIFFNESS_1 * H_NUMBER * H_NUMBER * H_NUMBER * H_NUMBER / 12.)
#define MU_NUMBER				(1./6.)
#define KSI_NUMBER				(4. * PI_NUMBER * MU_NUMBER)
//#define SPL_NUMBER				(SP_NUMBER / L_NUMBER)
//#define OMEGA_NUMBER			((SPL_NUMBER * SPL_NUMBER * SPL_NUMBER * SPL_NUMBER) * EI_NUMBER / KSI_NUMBER)
//#define OSC_PERIOD				100.
//#define OSC_PERIOD				(int)(0. * 2. * PI_NUMBER / OMEGA_NUMBER)
//#define OSC_AMPL				(0.01 / NUMBER_ACTIVE_NODES)
//#define FORCE_NUMBER			((DELTA_NUMBER * 3. * EI_NUMBER / (L_NUMBER * L_NUMBER) / cos(ANGLE_NUMBER)) / NUMBER_ACTIVE_NODES)
//#define FORCE_NUMBER			(0. * (DELTA_NUMBER * 3. * EI_NUMBER / (L_NUMBER * L_NUMBER)) / NUMBER_ACTIVE_NODES)
//#define OSC_AMPL				FORCE_NUMBER
//#define OSC_AMPL				0.


//#define CA_NUMBER				(0. * 1.e-3 / LS_STIFFNESS_0)
#define CA_NUMBER				5.
//#define EH_NUMBER				(LS_STIFFNESS_0 * LS_DX_0)
//#define FI_NUMBER				(0. * 10./ LS_STIFFNESS_0)
//#define FI_NUMBER				0.25
#define K_NUMBER				1.
#define R_EQ_NUMBER				0.75
//#define NC_NUMBER				(3. * PI_NUMBER * K_NUMBER * (3./4. * K_NUMBER + R_EQ_NUMBER) / (LS_DX_1 * LS_DX_1))
//#define N0_NUMBER				(4. * PI_NUMBER * R0_NUMBER * R0_NUMBER * NC_NUMBER / (LS_DX_0 * LS_DX_0))

#define RE_NUMBER				(DELTA_RO/((double)(INTERVAL*DOMAIN_FACT)))		// Test Re
//#define RE_NUMBER				(1. * 37.9583 * 0.5)		// For channel height 14 with 2 beads per disc - Ux/Vb
//#define RE_NUMBER				(0.2 * 44.0649)		// For channel height 14 with 2 beads per disc (offset discs) - Ux/Vb
//#define RE_NUMBER				(0.2 * 37.8315)		// For channel height 14 with beads - Ux/Vb
//#define RE_NUMBER				(0.4 * 3.0244)		// UPDATEFor channel height 30 with beads
//#define RE_NUMBER				(0.2 * 3.1193)		// UPDATEFor channel height 40 with ridges (no beads)
//#define RE_NUMBER				(0.2 * 3.2135)		// UPDATEFor channel height 45 with ridges and beads
#define GAMMA_NUMBER				(4. * RE_NUMBER * MU_NUMBER / CHANNEL_HIGHT / CHANNEL_HIGHT)
//#define GAMMA_NUMBER			(CA_NUMBER * EH_NUMBER / MU_NUMBER / R0_NUMBER)
//#define GAMMA_NUMBER			CA_NUMBER * LBB_SIGMA / MU_NUMBER / 10.
#define POISEUILLE_FLAG			0.
#define SHEAR_FLAG				0.
//#define GX_NUMBER				2. * GAMMA_NUMBER * MU_NUMBER / CHANNEL_HIGHT
#define GX_NUMBER				RE_NUMBER
#define GZ_NUMBER				0.

//#define LS_DE_0				(1. * EH_NUMBER * K_NUMBER * K_NUMBER * FI_NUMBER / N0_NUMBER)
#define LS_DE_0				0.
#define LS_DE_1				0.
#define LS_DE_2				0.

#define LS_DISIPATION_0		0.3
#define LS_DISIPATION_1		LS_DISIPATION_0
#define LS_DISIPATION_2		LS_DISIPATION_0
#define LS_DE_R_NUMBER		(0.005 * 50.)//(0.005 * 1000.)
#define LS_K_NUMBER			K_NUMBER
#define LS_R_EQ_NUMBER		R_EQ_NUMBER
#define LS_DE_C_NUMBER		LS_DE_R_NUMBER
#define LS_R_EQ_C_NUMBER	0.5

#define LS_SYMAX_0				1.
#define LS_SYMAX_1				1.
#define LS_SYMAX_2				1.
#define LS_SYMIN_0				-1.
#define LS_SYMIN_1				-1.
#define LS_SYMIN_2				-1.
#define LS_SY_NOISE				(0.*0.01)

#define LS_MASS_0				1.
#define LS_MASS_1				1.
#define LS_MASS_2				1.

#define LS_EQ_LENGTH_0			1.
#define LS_EQ_LENGTH_1			1.
#define LS_EQ_LENGTH_2			1.

#define LS_C_INI_X_1		0.
#define LS_C_INI_Y_1		0.
#define LS_C_INI_Z_1		0.
#define LS_C_INI_X_2		0.
#define LS_C_INI_Y_2		0.
#define LS_C_INI_Z_2		0.

// â”€â”€ BENDING STIFFNESS - TEMPLATE v2 â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
// The valve-opening knob. load/lsak.txt carries topology weights in [0,1]
// (1 = full-strength triplet, 0.5 = shared, 0.25 = double-length) and
// clConstants.cpp (clAK::ini) multiplies by the constant for the group of
// the triplet's CENTRE node, AN(i,1). Same group mapping as LS_STIFFNESS_*:
//     _0  ->  wall,  _1  ->  inlet leaflet 1,  _2  ->  the other 3 leaflets.
// Keep _1 == _2. Written by editConstDef.m from buildConfig.m's paper.K.
//
// Reference points, from the runs that exist:
//   4.0029618  = the 2026-08-13 run (Gamma 0.2, K 0.065, valveSoften 2)
//   3.9229026  = SimulationPaper    (Gamma 0.14, K 0.065) - measured
//                delta/DeltaV^(1/3) = 0.0825, i.e. 5x short of the paper.
// Those two differ by 2%, which is the whole reason v2 exists: in v1 the
// number was buried in the mesh where nobody could read it off.
#define LS_ANGLE_STIFFNESS_0 15
#define LS_ANGLE_STIFFNESS_1 0.08339503888	// x0.1 of 0.8339503888 (2026-09-26)
#define LS_ANGLE_STIFFNESS_2 0.08339503888	// x0.1 of 0.8339503888 (2026-09-26)

#define LB_COARSE_STEP		2
#define LB_COARSE_X0		0
#define LB_COARSE_Y0		-5
#define LB_COARSE_Z0		-5
#define LB_COARSE_NX		20
#define LB_COARSE_NY		20
#define LB_COARSE_NZ		20

//FOR TRACER PARTICLES
//#define TR_D0				0.0005
#define D0_TAU				((0.08128 * 4. * CHAN_HALF_HEIGHT * CHAN_HALF_HEIGHT) / TR_D0)
#define TR_TRUNC			6.
#define K_coeff				0.01
#define STOP_DIST_X			0.

#define LB_WALL_SHIFT			0

#define TOPWALLREFLECT
#define BOTTOMWALLREFLECT

//Boundary conditions

#define LB_BC_WALL_X0		FALSE
#define LB_BC_WALL_Y0		TRUE
#define LB_BC_WALL_Z0		TRUE
#define LB_BC_WALL_X1		LB_BC_WALL_X0
#define LB_BC_WALL_Y1		LB_BC_WALL_Y0
#ifdef SYMMETRIC_VESSEL
#define LB_BC_WALL_Z1		FALSE
#else
#define LB_BC_WALL_Z1		LB_BC_WALL_Z0
#endif // SYMMETRIC_VESSEL

#define LB_BC_PRESS_X0		TRUE
#define LB_BC_PRESS_Y0		FALSE
#define LB_BC_PRESS_Z0		FALSE
#define LB_BC_PRESS_X1		LB_BC_PRESS_X0
#define LB_BC_PRESS_Y1		LB_BC_PRESS_Y0
#define LB_BC_PRESS_Z1		LB_BC_PRESS_Z0

#define LB_BC_SYMMETRY_X0		FALSE
#define LB_BC_SYMMETRY_Y0		FALSE
#define LB_BC_SYMMETRY_Z0		FALSE
#define LB_BC_SYMMETRY_X1		FALSE
#define LB_BC_SYMMETRY_Y1		FALSE
#ifdef SYMMETRIC_VESSEL
#define LB_BC_SYMMETRY_Z1		TRUE
#else
#define LB_BC_SYMMETRY_Z1		FALSE
#endif // SYMMETRIC VESSEL

//#define LB_BC_COARSE_X0		FALSE
//#define LB_BC_COARSE_Y0		TRUE
//#define LB_BC_COARSE_Z0		TRUE
//#define LB_BC_COARSE_X1		!LB_BC_COARSE_X0
//#define LB_BC_COARSE_Y1		LB_BC_COARSE_Y0
//#define LB_BC_COARSE_Z1		LB_BC_COARSE_Z0

#define LB_BC_FREE_X0		FALSE
#define LB_BC_FREE_Y0		FALSE
#define LB_BC_FREE_Z0		FALSE
#define LB_BC_FREE_X1		LB_BC_FREE_X0
#define LB_BC_FREE_Y1		LB_BC_FREE_Y0
#define LB_BC_FREE_Z1		LB_BC_FREE_Z0

#define LB_BC_Out_X0		FALSE
#define LB_BC_Out_Y0		FALSE
#define LB_BC_Out_Z0		FALSE
#define LB_BC_Out_X1		FALSE			
#define LB_BC_Out_Y1		LB_BC_Out_Y0
#define LB_BC_Out_Z1		LB_BC_Out_Z0

#define LB_BC_PERIODIC_X	(!LB_BC_WALL_X0 && !LB_BC_PRESS_X0 && !LB_BC_FREE_X0 && !LB_BC_SYMMETRY_X0 && !LB_BC_Out_X0)
#define LB_BC_PERIODIC_Y	(!LB_BC_WALL_Y0 && !LB_BC_PRESS_Y0 && !LB_BC_FREE_Y0 && !LB_BC_SYMMETRY_Y0 && !LB_BC_Out_Y0)
#define LB_BC_PERIODIC_Z	(!LB_BC_WALL_Z0 && !LB_BC_PRESS_Z0 && !LB_BC_FREE_Z0 && !LB_BC_SYMMETRY_Z0 && !LB_BC_Out_Z0)


//Wall velocity
#define LB_UX0_X			0.
#define LB_UX0_Y			0.
#define LB_UX0_Z			0.
#define LB_UX1_X			0.
#define LB_UX1_Y			0.
#define LB_UX1_Z			0.

#define LB_UY0_X			(-SHEAR_FLAG * GAMMA_NUMBER * CHANNEL_HIGHT)
#define LB_UY0_Y			0.
#define LB_UY0_Z			0.
#define LB_UY1_X			(SHEAR_FLAG * GAMMA_NUMBER * CHANNEL_HIGHT)
#define LB_UY1_Y			0.
#define LB_UY1_Z			0.

#define LB_UZ0_X			0.
#define LB_UZ0_Y			0.
#define LB_UZ0_Z			0.
#define LB_UZ1_X			0.
#define LB_UZ1_Y			0.
#define LB_UZ1_Z			0.

// Inlet / Outlet pressure
#define LB_P_X0				(1./3. + DELTA_RO_WALL)
#define LB_P_X1				(1./3. - DELTA_RO_WALL)
#define LB_P_Y0				(1./3.)
#define LB_P_Y1				(1./3.)
#define LB_P_Z0				(1./3.)
#define LB_P_Z1				(1./3.)

// Inlet / Outlet density
#define LB_RO_X0a			(3. * LB_P_X0)
#define LB_RO_X1a			(3. * LB_P_X1)
#define LB_RO_Y0a			(3. * LB_P_Y0)
#define LB_RO_Y1a			(3. * LB_P_Y1)
#define LB_RO_Z0a			(3. * LB_P_Z0)
#define LB_RO_Z1a			(3. * LB_P_Z1)

// Inlet / Outlet density
#define LB_RO_X0			(LB_RO_X0a - 1.)
#define LB_RO_X1			(LB_RO_X1a - 1.)
#define LB_RO_Y0			(LB_RO_Y0a - 1.)
#define LB_RO_Y1			(LB_RO_Y1a - 1.)
#define LB_RO_Z0			(LB_RO_Z0a - 1.)
#define LB_RO_Z1			(LB_RO_Z1a - 1.)

#ifdef LS_DISSIPATION
#define LS_DISSIPATION_FLAG		1.
#else
#define LS_DISSIPATION_FLAG		0.
#endif //LS_DISSIPATION

#define LB_LAMDA_NUMBER			(-1./LB_TAU_NUMBER)
#define LB_LAMDAB_NUMBER		(-1./LB_TAUB_NUMBER)
#define LB_RHO_NUMBER			1.
#define LB_FLUID_COARSE			1


#define ARRAY_INCREASE_SIZE		100
#define NUMBER_DIMENSIONS		3
#define LB_NUMBER_COMPONENTS	3
#define LB_NUMBER_CONNECTIONS	19
#define LS_NUMBER_CONNECTIONS	20
#define LS_MAX_SF				LS_NUMBER_CONNECTIONS
#define LS_ARRAY_INCREASE_STEP	3
#define LS_SEARCH_ACTIVE_SIZE	3
#define LS_SEARCH_ACTIVE_STEP	3
#define LS_SEARCH_REPULSIVE_STEP	2
#define LS_SEARCH_REPULSIVE_SIZE	2

//#define OUTPUT_LB_FILEDS		10
//Modified
//#define OUTPUT_LB_FILEDS		15
#define OUTPUT_LB_FILEDS		10
//Modified
#define OUTPUT_LS_FILEDS		45
#define OUTPUT_MAX_LINES		10000


#define LB_NO_BOUNDARY			0
#define LB_BOUNDARY_LIQUID		1
#define LB_BOUNDARY_SOLID		2

#define LB_SOLID				0
#define LB_LIQUID1				1
#define LB_LIQUID2				2

#define LB_NO_BOUNDARY_LINK		0x0
#define LB_BOUNDARY_LINK_1		0x1
#define LB_BOUNDARY_LINK_2		0x2
#define LB_BOUNDARY_LINK_3		0x4
#define LB_BOUNDARY_LINK_4		0x8
#define LB_BOUNDARY_LINK_5		0x10
#define LB_BOUNDARY_LINK_6		0x20
#define LB_BOUNDARY_LINK_7		0x40
#define LB_BOUNDARY_LINK_8		0x80
#define LB_BOUNDARY_LINK_9		0x100
#define LB_BOUNDARY_LINK_10		0x200
#define LB_BOUNDARY_LINK_11		0x400
#define LB_BOUNDARY_LINK_12		0x800
#define LB_BOUNDARY_LINK_13		0x1000
#define LB_BOUNDARY_LINK_14		0x2000
#define LB_BOUNDARY_LINK_15		0x4000
#define LB_BOUNDARY_LINK_16		0x8000
#define LB_BOUNDARY_LINK_17		0x10000
#define LB_BOUNDARY_LINK_18		0x20000

#define LS_REGULAR_NODE			0x0
#define LS_VOLUME_NODE			0x1
#define LS_FLAT_NODE			0x2
#define LS_PERIODIC_NODE		0x4
#define LS_ACTIVE_NODE			0x8
#define LS_STATIC_NODE			0x10
#define LS_RIGID_NODE			0x20
#define LS_WALL_INTER_NODE		0x40
#define LS_OUTER_NODE			0x80
#define LS_CUSTOM_NODE_0		0x100
#define LS_CUSTOM_NODE_1		0x200
#define LS_CUSTOM_NODE_2		0x400
#define LS_CUSTOM_NODE_3		0x800
#define LS_CUSTOM_NODE_4		0x1000

#ifdef RENAME_DUMP_FILES
#define FLAG_DUMP	TRUE
#else //RENAME_DUMP_FILES
#define FLAG_DUMP	FALSE
#endif //RENAME_DUMP_FILES

#ifdef SAVE_BIN_IN_LOAD
#define BIN_SAVE_DIR	"load" SLASH
#else //SAVE_BIN_IN_LOAD
#define BIN_SAVE_DIR	""
#endif //SAVE_BIN_IN_LOAD

#ifdef LB_USE_OPENMP
#define USE_OPENMP
#endif //LB_USE_OPENMP
#ifdef LS_USE_OPENMP
#define USE_OPENMP
#endif //LS_USE_OPENMP
#ifdef TR_USE_OPENMP
#define USE_OPENMP
#endif //TR_USE_OPENMP

#ifndef LS_SOLVER
#undef LS_STEPS
#define LS_STEPS				(1)
#endif //LS_SOLVER

#ifndef LB_COARSE
#undef LB_BC_COARSE_X0
#undef LB_BC_COARSE_Y0
#undef LB_BC_COARSE_Z0
#undef LB_BC_COARSE_X1
#undef LB_BC_COARSE_Y1		
#undef LB_BC_COARSE_Z1
#define LB_BC_COARSE_X0		FALSE
#define LB_BC_COARSE_Y0		FALSE
#define LB_BC_COARSE_Z0		FALSE
#define LB_BC_COARSE_X1		FALSE
#define LB_BC_COARSE_Y1		FALSE
#define LB_BC_COARSE_Z1		FALSE
#endif //LB_COARSE

#define OUTPUT_DIR				"results"
#define CURR_OUTPUT_FILE		"tout.txt"
#define UNIX_OUTPUT_FILE		"stdout.txt"
#define TMP_OUTPUT_FILE			"tout_tmp.txt"

#endif // !defined(__CONSTDEF_H__INCLUDED_)
