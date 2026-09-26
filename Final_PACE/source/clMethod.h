// clMethod.h: interface for the clMethod class.
//
// (c) Alexander Alexeev, 2006 
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CLMETHOD_H__INCLUDED_)
#define AFX_CLMETHOD_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "StdAfx.h"
#include <vector>
#include <map>
#include <algorithm>
#include "clProblem.h"
//#include "clSolver.h"
//#include "clGrid.h"
#include "clVariables.h"



class clMethod  
{
private:

public:

	double FtotalX, FtotalY;

	void ini()
	{

		v.update();
		vls.update();       ////////////////////////////////////????????????????
//		vlb.update();
		p.update();

//		c.Mass0 = Mass = v.mass();
	};

	void update()
	{

//		c.update();
//		g.update();

		v.update();

//		vls.update();
//		vlb.update();


		p.update();

//		Mass = v.mass();
		c.signal();
	};

	void run()
	{
		p.start();

		while(1)
		{
#ifdef _DEBUG
			p.debug();
#endif //_DEBUG

			solve();

			update();

			p.step();

			if(testFinish() == TRUE)	
			{
				p.finish();
				break;
			}
		}
	};

	BOOL testFinish()
	{
		
//#ifdef TR_SOLVER
//		if((double)vtr.C.OutputVal.x == 0.)
//			return TRUE;
//#endif
		
		//static double CG0yini = vls.avr.CG(0).y;
#ifdef MIXSTOP
		if((double)vtr.C.OutputVal.x <= 0.29 && c.t.Time > c.t.DumpTimeStepSave)
			return TRUE;
#endif // MIXSTOP
		if(c.t.Time >= c.t.StopTime)
			return TRUE;
//		return FALSE;

//		vls.avr.calcCG();
//		if(vls.avr.CG(0).x < c.lb.nX - 0.5)
		//{
		//	double dCG0y = vls.avr.CG(0).y - CG0yini;
		//	for(int i = 0; i < c.ls.nN; i++)
		//		if(c.ls.O(i) == 0)
		//			vls.C(i).y -= dCG0y;
		//}

		//if(vls.avr.CG(0).x - vls.avr.CG(0).y * 2. >= c.lb.nX)
		//	return TRUE;
		//static double TimeV0 = -1.;
		//if(c.t.Time >= 10000.)
		//{
		//	if(vls.avr.V(0).abs() < 1.e-4)
		//	{
		//		if(TimeV0 < 0)
		//			TimeV0 = c.t.Time;
		//		if(c.t.Time - TimeV0 > 2000.)
		//			return TRUE;
		//	}
		//	else
		//	{
		//		TimeV0 = -1.;
		//	}
		//}

		return FALSE;
	};

	void solve()
	{
#ifdef LB_SOLVER
		solveLB();
		vlb.update();
#endif //LB_SOLVER

#ifdef LS_SOLVER
		vls.V.avr.zeros();

		for(int i = 0; i < c.t.lsSteps; i++)
		{
			solveLS(i);
			vls.update();
			c.update();

			for(int j = 0; j < c.ls.nN; j++)
				vls.V.avr(j) += vls.V(j);
		}

//		vls.D.update();
		for(int i = 0; i < c.ls.nN; i++)
			vls.V.avr(i) /= (double)c.t.lsSteps;
#else
		c.update();
#endif //LS_SOLVER
//
//#ifdef TR_SOLVER
//		if(c.t.Time >= (PARTICLE_RELEASE_TIME))
//		{	
//			vtr.update1();
//			for(int i = 0; i < c.t.trSteps; i++)
//			{
//				vtr.update2();
//				solveTr();
//				vtr.A.update();
//			}
//		}
//#endif //TR_SOLVER

		
#ifdef TR_SOLVER

		vtr.update1();
		for(int i = 0; i < c.t.trSteps; i++)
		{
			vtr.update2();
			solveTr();
			vtr.A.update();
		}

#endif //TR_SOLVER

	}

////////////////////////////////////////////////////////////////////////////
////////////////////////Tr solver///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	void solveTr()
	{
		double dTtr = c.t.dTtr;
		double D = c.tr.D0;
		double diff = sqrt(2.*D*dTtr);

		int i;
#ifdef TR_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //TR_USE_OPENMP
		for(i = 0; i < c.tr.remain; i++)
		{
			clVector3Dd d_w;
			clVector3Dd d_u;

#ifdef LB_SOLVER
			d_u = vtr.V(i) * dTtr;
#else
			d_u = (0,0,0);
#endif

#ifdef TR_USE_DIFFUSION
			d_w = vtr.W(i)*diff;
#else
			d_w = (0,0,0);
#endif

			vtr.C(i) += d_u + d_w;
			vtr.C.d(i) = d_u + d_w;
		}
	}


//	void solveTr()
//	{
//		double dTtr = c.t.dTtr;
//		double D = c.tr.D0;
//		double diff = sqrt(2.*D*dTtr);
//
//		int it;
//#ifdef TR_USE_OPENMP
//		#pragma omp parallel for default(shared) private(it)
//#endif //TR_USE_OPENMP
//		for(it = 0; it < c.tr.remain; it++)
//		{
//			int i = vtr.S(it);
//			
//			if(vtr.A(i) != 0)
//				continue;
//
//			clVector3Dd d_w;
//			clVector3Dd d_u;
//
//#ifdef LB_SOLVER
//			d_u = vtr.V(it) * dTtr;
//#else
//			d_u = (0,0,0);
//#endif
//
//#ifdef TR_USE_DIFFUSION
//			d_w = vtr.W(it)*diff;
//#else
//			d_w = (0,0,0);
//#endif
//
//			vtr.C(i) += d_u + d_w;
//			vtr.C.d(i) = d_u + d_w;
//		}
//	}


////////////////////////////////////////////////////////////////////////////
////////////////////////LS solver///////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	void solveLS(int lsStep)
	{
		solveVerlet1();

		vls.C.shift();
        
		vls.F.updateW();

		vls.C.zerosLS();
		c.ls.K.zerosE();
		c.ls.AK.zerosE();
		vls.F.zerosP();
		vls.F.zerosC();
		vls.F.zerosS();
		

		solveLSspring();
		solveLSangle();
		solveLSties();
#ifdef LS_BEADS
		vls.F.updateE();
#endif // LS_BEADS
		vls.P.updateF();
		vls.F.updateA(lsStep);

#ifdef LS_DEBUG
		vls.F.save2file();
#endif //LS_DEBUG

		int i;
#ifdef LS_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //LS_USE_OPENMP
		for(i = 0; i < c.ls.nN; i++)
		{
			//vls.F.b(i).y = 0.;
#ifdef LS_NO_COUPLING
			vls.F(i) = vls.F.s(i) + vls.F.c(i) + vls.F.e(i);
#else //LS_NO_COUPLING
			vls.F(i) = vls.F.s(i) + vls.F.b(i) + vls.F.c(i) + vls.F.e(i);
#endif //LS_NO_COUPLING
		}

		solveRigidBody();
		solveRigidBodyCorrection();

		solveVerlet2();

#ifdef LS_DEBUG
		vls.save2file();
#endif //LS_DEBUG

	}
	// Prescribed radial wall displacement at axial position localX (measured from the
	// inlet valve root, valveMinX) and its time derivative.
	// WALL_LAW_PAPER: the published model, r = R - A/2 (1 - cos(2 pi x/L)) sin(w t),
	// applied only between the valve roots (0 <= x <= L); the outlet valve sits on a
	// stationary extension and the inlet valve where the wall barely moves. tau = t/T
	// with contraction first.
	// Otherwise the legacy law A sin(2 pi x/lambda) sin(w t), lambda = 2L, which moves
	// the outlet-valve region in anti-phase and starts with expansion.
	double wallMod(double localX)
	{
		double w = 2.0 * PI_NUMBER / RAMP_PERIOD;
#ifdef WALL_LAW_PAPER
		double L = 0.5 * c.ls.waveLength;
		if(localX < 0. || localX > L)
			return 0.;
		return -0.5 * OSC_AMPLITUDE * (1. - cos(2.0 * PI_NUMBER * localX / L)) * sin(w * c.t.Time);
#else
		return OSC_AMPLITUDE * sin(2.0 * PI_NUMBER * (localX / c.ls.waveLength)) * sin(w * c.t.Time);
#endif
	}
	double wallRate(double localX)
	{
		double w = 2.0 * PI_NUMBER / RAMP_PERIOD;
#ifdef WALL_LAW_PAPER
		double L = 0.5 * c.ls.waveLength;
		if(localX < 0. || localX > L)
			return 0.;
		return -0.5 * OSC_AMPLITUDE * (1. - cos(2.0 * PI_NUMBER * localX / L)) * w * cos(w * c.t.Time);
#else
		return OSC_AMPLITUDE * sin(2.0 * PI_NUMBER * (localX / c.ls.waveLength)) * w * cos(w * c.t.Time);
#endif
	}

	void solveVerlet1()
	{
		int i;
#ifdef SET_ORBIT
		int o;
#endif // SET_ORBIT
		double dTls = c.t.dTls; 
#ifdef LS_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //LS_USE_OPENMP
		for(i = 0; i < c.ls.nN; i++)
		{
			if(c.ls.Tno_c(i) == TRUE)
//			if(::c.ls.T(i) & LS_CUSTOM_NODE_0)
			{
//				vls.V(i) = 0.;
				continue;
			}
#ifdef SET_ORBIT
			o = c.ls.O(i);
			if(o < N_BEAD)
			{
				double oscsin = ::c.t.oscillator1sin();
				double osccos = ::c.t.oscillator1cos();
				if(o % N_BEADS_DISK == 0)
				{
					vls.V(i).x = BEAD_VELOCITY * oscsin;
					vls.V(i).y = -1. * BEAD_VELOCITY * osccos;
					vls.V(i).z = 0;
				}
				else
				{
					vls.V(i).x = -BEAD_VELOCITY * oscsin;
					vls.V(i).y = 1. * BEAD_VELOCITY * osccos;
					vls.V(i).z = 0;
				}
			}
			else
			{
				vls.V(i).x = 0.;
				vls.V(i).y = 0.;
				vls.V(i).z = 0.;
			}

			//int o = c.ls.O(i);	// Bead object number
			//int od = o + N_BEAD;	// Determine object number of corresponding disk
			////clVector3Dd R_bd = vls.avr.CG(o) - vls.avr.CG(od);	// Unit vector from disc center to bead center
			//clVector3Dd R_bd = vls.avr.CG(o) - c.ls.Cdisc(o);	// Unit vector from disc center to bead center
			//R_bd.z = 0;
			//R_bd = R_bd / R_bd.abs();	// Unit vector from disc center to bead center
			//double velmag = BEAD_OMEGA * ORBIT_RADIUS;
			//vls.V(i).x = -R_bd.y * velmag;
			//vls.V(i).y = R_bd.x * velmag;
			//vls.V(i).z = 0;

			//int o = c.ls.O(i);	// Bead object number
			//int od = o + N_BEAD;	// Determine object number of corresponding disk
			//clVector3Dd R_bd = vls.C(i) - vls.avr.CG(od);	// Vector from disc center to bead center
			//vls.V(i).x = R_bd.y * BEAD_OMEGA;
			//vls.V(i).y = -R_bd.x * BEAD_OMEGA;
			//vls.V(i).z = 0;
#else
#ifdef STATIC_LS
			vls.V(i).x = 0.;
			vls.V(i).y = 0.;
			vls.V(i).z = 0.;

#elif	defined TRANSLATE_LS
			vls.V(i).x = TRANSLATE_LS;
			vls.V(i).y = 0.;
			vls.V(i).z = 0.;

#elif defined LEAFLETS
#ifdef STATIC_VESSEL
			if(c.ls.O(i) == 0 || ::c.ls.T(i) & LS_CUSTOM_NODE_0 || ::c.ls.T(i) & LS_CUSTOM_NODE_1)
#else
			if(::c.ls.T(i) & LS_CUSTOM_NODE_0 || ::c.ls.T(i) & LS_CUSTOM_NODE_1)
#endif // STATIC_VESSEL
			{
				//vls.V(i) = 0.;
				double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;

				// Peristaltic modulation as a function of location and time, evaluated
				// everywhere (no X-range gating): staticness is handled entirely by
				// LS_STATIC_NODE, set at geometry time on entrance/exit vessel wall
				// nodes outside the valve footprint (GeometryCreatorAmir.m's
				// entrance_static_flag block). The previous X-range clamp
				// (isOutsideValveRegion / x<valveMinX / x>=valveMaxX) zeroed
				// modulation for any node at or beyond valveMaxX - which is exactly
				// where the SECOND valve's leaflets live (lso==3,4 span
				// [valveMaxX, valveMaxX+leaflet length]), leaving them structurally
				// undriven while the first valve's leaflets (inside [valveMinX,
				// valveMaxX)) worked normally.
				double localX = vls.C(i).x - c.ls.valveMinX;
				double radiusModulation = wallMod(localX);
				if(::c.ls.T(i) & LS_STATIC_NODE)
				{
					radiusModulation = 0.0;
				}

				// Calculate current radial position
				double currentY = vls.C(i).y - midy;
				double currentZ = vls.C(i).z - midz;
				double currentRadius = sqrt(currentY*currentY + currentZ*currentZ);

				// Target radius: node's own initial (reference) radius + modulation from wave.
				// Using each node's own baseline (instead of one vessel-wide radius) preserves
				// non-circular anchor shapes (e.g. leaflet root crescent) instead of snapping
				// them onto a circle.
				double initY = vls.C.c0(i).y - midy;
				double initZ = vls.C.c0(i).z - midz;
				double baseRadius = sqrt(initY*initY + initZ*initZ);
				double targetRadius = baseRadius + radiusModulation;

				// Avoid division by zero
				if(currentRadius < 1.e-6) {
					vls.V(i).x = 0.0;
					vls.V(i).y = 0.0;
					vls.V(i).z = 0.0;
				}
				else {
					// Scale factor: how much to stretch/shrink from current to target radius
					double scaleFactor = targetRadius / currentRadius;

					// Target positions in Y/Z
					double targetY = midy + currentY * scaleFactor;
					double targetZ = midz + currentZ * scaleFactor;

					// Velocity = analytic wall velocity (what the fluid must see, via
					// vls.V.avr) + a position-tracking correction that keeps the node on
					// its target. Using only (target-current)/dt gave the fluid a zero
					// wall velocity, because solveVerlet2 then found the node already on
					// target (see WALL_VELOCITY note in solveVerlet2).
					double drdt = (::c.ls.T(i) & LS_STATIC_NODE) ? 0.0 : wallRate(localX);
					vls.V(i).x = 0.0;
					vls.V(i).y = drdt * currentY / currentRadius + (targetY - vls.C(i).y) / dTls;
					vls.V(i).z = drdt * currentZ / currentRadius + (targetZ - vls.C(i).z) / dTls;
				}				
			}
			else
			{
				vls.V(i) += vls.F(i) * 0.5 * dTls / c.ls.M(i);
#ifdef SYMMETRIC_VESSEL
				if(c.ls.O(i) >= NUM_VESSELS && vls.C(i).z >= c.lb.nZ - 0.5 - SYMMETRIC_GAP_DIST && vls.V(i).z > 0.){
					vls.V(i).z = 0.;
				}
#endif // SYMMETRIC_VESSEL
#ifdef CONSTRAIN_MIDLINE
				double Cz = vls.C(i).z;
				double midline = c.lb.nZ/2. - 0.5;
				double dist_midline = fabs(midline - Cz);
				if(c.ls.O(i) >= NUM_VESSELS && dist_midline <= SYMMETRIC_GAP_DIST){
					if((Cz <= c.lb.nZ/2. - 0.5 && vls.V(i).z > 0.) || (Cz >= c.lb.nZ/2. - 0.5 && vls.V(i).z < 0.)){
						vls.V(i).z = 0.;
					}
				}
#endif // CONSTRAIN_MIDLINE

#ifdef OSCILLATING_VESSEL

				// Position-based (matches the LEAFLETS branch above): the old version
				// prescribed a radial VELOCITY ("speed") each step and let it integrate
				// into position. OSC_AMPLITUDE is a radius amplitude, not a velocity, so
				// without a frequency-scaled integration that drifted the radius outward
				// without bound (grew quadratically with time - verified against
				// Simulation/pressureDiff_0/results/*, extent 20->380 over 23 dumps).
				// Prescribing the target radius directly and driving velocity from
				// (target-current)/dt, like the leaflet anchors, keeps it bounded.
				// No X-range gating here either (see the LEAFLETS branch above) -
				// LS_STATIC_NODE (set on entrance/exit vessel wall nodes in
				// GeometryCreatorAmir.m's entrance_static_flag block) is what keeps
				// those nodes stationary now.
				double localX_osc = vls.C(i).x - c.ls.valveMinX;
				double radiusModulation_osc = wallMod(localX_osc);
				if(::c.ls.T(i) & LS_STATIC_NODE)
				{
					radiusModulation_osc = 0.0;
				}

				if(c.ls.O(i) == (NUM_VESSELS-1))
				{
					double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;

					double currentY_osc = vls.C(i).y - midy;
					double currentZ_osc = vls.C(i).z - midz;
					double currentRadius_osc = sqrt(currentY_osc*currentY_osc + currentZ_osc*currentZ_osc);

					double initY_osc = vls.C.c0(i).y - midy;
					double initZ_osc = vls.C.c0(i).z - midz;
					double baseRadius_osc = sqrt(initY_osc*initY_osc + initZ_osc*initZ_osc);
					double targetRadius_osc = baseRadius_osc + radiusModulation_osc;

					if(currentRadius_osc < 1.e-6) {
						vls.V(i).x = 0.0;
						vls.V(i).y = 0.0;
						vls.V(i).z = 0.0;
					}
					else {
						double scaleFactor_osc = targetRadius_osc / currentRadius_osc;
						double targetY_osc = midy + currentY_osc * scaleFactor_osc;
						double targetZ_osc = midz + currentZ_osc * scaleFactor_osc;

						// analytic wall velocity + position-tracking correction (see LEAFLETS branch)
						double drdt_osc = (::c.ls.T(i) & LS_STATIC_NODE) ? 0.0 : wallRate(localX_osc);
						vls.V(i).x = 0.0;
						vls.V(i).y = drdt_osc * currentY_osc / currentRadius_osc + (targetY_osc - vls.C(i).y) / dTls;
						vls.V(i).z = drdt_osc * currentZ_osc / currentRadius_osc + (targetZ_osc - vls.C(i).z) / dTls;
					}
				}
#endif
			}
#else
			vls.V(i) += vls.F(i) * 0.5 * dTls / c.ls.M(i);

#endif // STATIC_LS
#endif // SET_ORBIT
			//int o = c.ls.O(i);
			//if (o < N_BEAD && vls.C(i).z < -0.49)
			//{
			//	vls.V(i) = clVector3Dd(0., 0., 0.);
			//}

			clVector3Dd dC = vls.V(i) * dTls;
			vls.C.d(i) = dC;
			vls.C(i) += dC;
		}

		//clVector3Dd dCtot = clVector3Dd(0., ::c.t.oscillator_ratio(::c.t.Time + (double)lsStep * ::c.t.dTls) - ::c.t.oscAmpl1, 0.);
		//for(int i = 0; i < ::c.ls.nN; i++)
		//{
		//	if(::c.ls.T(i) & LS_CUSTOM_NODE_0)
		//	{
		//		clVector3Dd C1 = vls.C.c0(i) + dCtot;
		//		clVector3Dd dC = C1 - vls.C(i);
		//		//clVector3Dd dV = dC / ::c.t.dTls - vls.V(i);
		//		vls.V(i) = dC / ::c.t.dTls;
		//		vls.C.d(i) = dC;
		//		vls.C(i) = C1;
		//	}
		//}

		//double dy = c.t.oscillator1(0.5 * c.Pi) - c.t.oscAmpl1;
//		double dy = c.t.oscillator_ratio() - c.t.oscAmpl1;
//		for(i = 0; i < c.ls.nN; i++)
//		{
//			if(c.ls.Ts(i) == TRUE)
//			{
//				//vls.V(i).y = vls.C.d(i).y / dTls;
//				//vls.C(i).y += vls.C.d(i).y;
//
//				double y1 = vls.C.c0(i).y + dy;
//				double dy = y1 - vls.C(i).y;
//				double Vy1 = dy / dTls, dVy = Vy1 - vls.V(i).y;
//				clVector3Dd Ftot;
////				if(vls.C(i).x > -0.5)
//					Ftot = clVector3Dd(0., dVy, 0.) * c.ls.M(i) / dTls;
//				vls.F.e(i) = Ftot - vls.F.b_m1(i) - vls.F.c(i);
//				vls.F(i) = Ftot;
//
//				//vls.V(i).y = dy / dTls;
////				vls.C.d(i).y = dy;
//				//vls.C(i).y = y1;
//			}
//		}

		for(i = 0; i < c.ls.nPN; i++)
		{
			int n0 = c.ls.PN(i,0), n1 = c.ls.PN(i,1);
			vls.F.b(n1) += vls.F.b(n0);
			vls.F.b(n0) = 0.;
			vls.C(n0) = vls.C(n1) + c.ls.PNC(i);
			vls.V(n0) = vls.V(n1);
		}
	}
	void solveVerlet2()
	{
		int i;
#ifdef SET_ORBIT
		int o;
#endif // SET_ORBIT
		double dTls = c.t.dTls; 
		double maxV = 0., V;
#ifdef LS_USE_OPENMP
		#pragma omp parallel for default(shared) private(i, V)	// maxV: serial pass below (MSVC = OpenMP 2.0, no max reduction)

#endif //LS_USE_OPENMP
		for(i = 0; i < c.ls.nN; i++)
		{
			if(c.ls.Tno_c(i) == TRUE)
//			if(::c.ls.T(i) & LS_CUSTOM_NODE_0)
				continue;
#ifdef SET_ORBIT
			o = c.ls.O(i);
			if(o < N_BEAD)
			{
				double oscsin = ::c.t.oscillator1sin();
				double osccos = ::c.t.oscillator1cos();
				if(o % N_BEADS_DISK == 0)
				{
					vls.V(i).x = BEAD_VELOCITY * oscsin;
					vls.V(i).y = -1. * BEAD_VELOCITY * osccos;
					vls.V(i).z = 0;
				}
				else
				{
					vls.V(i).x = -BEAD_VELOCITY * oscsin;
					vls.V(i).y = 1. * BEAD_VELOCITY * osccos;
					vls.V(i).z = 0;
				}
			}
			else
			{
				vls.V(i).x = 0.;
				vls.V(i).y = 0.;
				vls.V(i).z = 0.;
			}


			//int o = c.ls.O(i);	// Bead object number
			//int od = o + N_BEAD;	// Determine object number of corresponding disk
			////clVector3Dd R_bd = vls.avr.CG(o) - vls.avr.CG(od);	// Unit vector from disc center to bead center
			//clVector3Dd R_bd = vls.avr.CG(o) - c.ls.Cdisc(o);	// Unit vector from disc center to bead center
			//R_bd.z = 0;
			//R_bd = R_bd / R_bd.abs();	// Unit vector from disc center to bead center
			//double velmag = BEAD_OMEGA * ORBIT_RADIUS;
			//vls.V(i).x = -R_bd.y * velmag;
			//vls.V(i).y = R_bd.x * velmag;
			//vls.V(i).z = 0;

			//int o = c.ls.O(i);	// Bead object number
			//int od = o + N_BEAD;	// Determine object number of corresponding disk
			//clVector3Dd R_bd = vls.C(i) - vls.avr.CG(od);	// Vector from disc center to bead center
			//vls.V(i).x = R_bd.y * BEAD_OMEGA;
			//vls.V(i).y = -R_bd.x * BEAD_OMEGA;
			//vls.V(i).z = 0;
#else
#ifdef STATIC_LS
			vls.V(i).x = 0.;
			vls.V(i).y = 0.;
			vls.V(i).z = 0.;

#elif	defined TRANSLATE_LS
			vls.V(i).x = TRANSLATE_LS;
			vls.V(i).y = 0.;
			vls.V(i).z = 0.;

#elif defined LEAFLETS
#ifdef STATIC_VESSEL
			if(c.ls.O(i) == 0 || ::c.ls.T(i) & LS_CUSTOM_NODE_0 || ::c.ls.T(i) & LS_CUSTOM_NODE_1)
#else
			if(::c.ls.T(i) & LS_CUSTOM_NODE_0 || ::c.ls.T(i) & LS_CUSTOM_NODE_1)
#endif // STATIC_VESSEL
			{
				//vls.V(i) = 0.;
				double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;

				double localX_v2 = vls.C(i).x - c.ls.valveMinX;
				double radiusModulation = wallMod(localX_v2);
				if(::c.ls.T(i) & LS_STATIC_NODE)
				{
					radiusModulation = 0.0;
				}

				// Calculate current radial position
				double currentY_v2 = vls.C(i).y - midy;
				double currentZ_v2 = vls.C(i).z - midz;
				double currentRadius_v2 = sqrt(currentY_v2*currentY_v2 + currentZ_v2*currentZ_v2);

				// Target radius: node's own initial (reference) radius + modulation from wave.
				// Using each node's own baseline (instead of one vessel-wide radius) preserves
				// non-circular anchor shapes (e.g. leaflet root crescent) instead of snapping
				// them onto a circle.
				double initY_v2 = vls.C.c0(i).y - midy;
				double initZ_v2 = vls.C.c0(i).z - midz;
				double baseRadius_v2 = sqrt(initY_v2*initY_v2 + initZ_v2*initZ_v2);
				double targetRadius_v2 = baseRadius_v2 + radiusModulation;

				// Avoid division by zero
				if(currentRadius_v2 < 1.e-6) {
					vls.V(i).x = 0.0;
					vls.V(i).y = 0.0;
					vls.V(i).z = 0.0;
				}
				else {
					// Scale factor: how much to stretch/shrink from current to target radius
					double scaleFactor_v2 = targetRadius_v2 / currentRadius_v2;

					// Target positions in Y/Z
					double targetY_v2 = midy + currentY_v2 * scaleFactor_v2;
					double targetZ_v2 = midz + currentZ_v2 * scaleFactor_v2;

					// WALL_VELOCITY: this is the velocity accumulated into vls.V.avr and
					// handed to the fluid. Position was already set in solveVerlet1, so
					// (target-current)/dt is ~0 here and the fluid saw a motionless wall.
					// Use the analytic rate of the prescribed motion instead.
					double drdt_v2 = (::c.ls.T(i) & LS_STATIC_NODE) ? 0.0 : wallRate(localX_v2);
					vls.V(i).x = 0.0;
					vls.V(i).y = drdt_v2 * currentY_v2 / currentRadius_v2;
					vls.V(i).z = drdt_v2 * currentZ_v2 / currentRadius_v2;
				}	
			}
			else
			{
				vls.V(i) += vls.F(i) * 0.5 * dTls / c.ls.M(i);
#ifdef SYMMETRIC_VESSEL
				if(c.ls.O(i) >= NUM_VESSELS && vls.C(i).z >= c.lb.nZ - 0.5 - SYMMETRIC_GAP_DIST && vls.V(i).z > 0.){
					vls.V(i).z = 0.;
				}
#endif // SYMMETRIC_VESSEL
#ifdef CONSTRAIN_MIDLINE
				double Cz = vls.C(i).z;
				double midline = c.lb.nZ/2. - 0.5;
				double dist_midline = fabs(midline - Cz);
				if(c.ls.O(i) >= NUM_VESSELS && dist_midline <= SYMMETRIC_GAP_DIST){
					if((Cz <= c.lb.nZ/2. - 0.5 && vls.V(i).z > 0.) || (Cz >= c.lb.nZ/2. - 0.5 && vls.V(i).z < 0.)){
						vls.V(i).z = 0.;
					}
				}
#endif // CONSTRAIN_MIDLINE

#ifdef OSCILLATING_VESSEL

				// Position-based (matches solveVerlet1's version and the LEAFLETS
				// branch): prescribing target radius directly and driving velocity
				// from (target-current)/dt keeps the radius bounded, unlike the old
				// speed/velocity prescription which let the radius drift outward
				// without bound (quadratic growth in time - verified against
				// Simulation/pressureDiff_0/results/*, extent 20->380 over 23 dumps).
				// No X-range gating (see solveVerlet1) - LS_STATIC_NODE handles it.
				double localX_osc2 = vls.C(i).x - c.ls.valveMinX;
				double radiusModulation_osc2 = wallMod(localX_osc2);
				if(::c.ls.T(i) & LS_STATIC_NODE)
				{
					radiusModulation_osc2 = 0.0;
				}

				if(c.ls.O(i) == (NUM_VESSELS-1))
				{
					double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;

					double currentY_osc2 = vls.C(i).y - midy;
					double currentZ_osc2 = vls.C(i).z - midz;
					double currentRadius_osc2 = sqrt(currentY_osc2*currentY_osc2 + currentZ_osc2*currentZ_osc2);

					double initY_osc2 = vls.C.c0(i).y - midy;
					double initZ_osc2 = vls.C.c0(i).z - midz;
					double baseRadius_osc2 = sqrt(initY_osc2*initY_osc2 + initZ_osc2*initZ_osc2);
					double targetRadius_osc2 = baseRadius_osc2 + radiusModulation_osc2;

					if(currentRadius_osc2 < 1.e-6) {
						vls.V(i).x = 0.0;
						vls.V(i).y = 0.0;
						vls.V(i).z = 0.0;
					}
					else {
						double scaleFactor_osc2 = targetRadius_osc2 / currentRadius_osc2;
						double targetY_osc2 = midy + currentY_osc2 * scaleFactor_osc2;
						double targetZ_osc2 = midz + currentZ_osc2 * scaleFactor_osc2;

						// WALL_VELOCITY: analytic rate, see LEAFLETS branch above
						double drdt_osc2 = (::c.ls.T(i) & LS_STATIC_NODE) ? 0.0 : wallRate(localX_osc2);
						vls.V(i).x = 0.0;
						vls.V(i).y = drdt_osc2 * currentY_osc2 / currentRadius_osc2;
						vls.V(i).z = drdt_osc2 * currentZ_osc2 / currentRadius_osc2;
					}
				}
#endif
			}
#else
			vls.V(i) += vls.F(i) * 0.5 * dTls / c.ls.M(i);

#endif // STATIC_LS
#endif // SET_ORBIT

			V = vls.V(i).abs();
			//int o = c.ls.O(i);
			//if (o < N_BEAD && vls.C(i).z < -0.49)
			//{
			//	vls.V(i) = clVector3Dd(0., 0., 0.);
			//}
#ifndef LS_USE_OPENMP
			if(V > maxV) maxV = V;
#endif //LS_USE_OPENMP
		}
#ifdef LS_USE_OPENMP
		for(i = 0; i < c.ls.nN; i++)
			if(c.ls.Tno_c(i) == FALSE && vls.V(i).abs() > maxV)
				maxV = vls.V(i).abs();
#endif //LS_USE_OPENMP
		vls.V.max = maxV;

		for(i = 0; i < c.ls.nPN; i++)
		{
			int n0 = c.ls.PN(i,0), n1 = c.ls.PN(i,1);
			vls.V(n0) = vls.V(n1);
		}
	}
#ifndef LS_STRETCH_EDGE_FACTOR
#define LS_STRETCH_EDGE_FACTOR	1.0
#endif
#ifndef LS_STRETCH_EDGE_RAMP
#define LS_STRETCH_EDGE_RAMP	24.0
#endif
	// ---------------------------------------------------------------------------
	// Free-edge softening of the leaflet in-plane springs (2026-09-26, Amir's request).
	// Every leaflet spring's stiffness K is multiplied once, at the first spring solve, by
	//     f = LS_STRETCH_EDGE_FACTOR + (1 - LS_STRETCH_EDGE_FACTOR) * min(de / LS_STRETCH_EDGE_RAMP, 1)
	// de = distance of the spring midpoint from the leaflet's free edge in the reference mesh.
	// Free edge = boundary edges of the leaflet's triangles with at least one unanchored end
	// (anchored = LS_CUSTOM_NODE_0/1, the clamp band), i.e. the 12 segments / 13 nodes that
	// CompareOldData.py uses. Damping Ds is left unchanged. FACTOR = 1 (default) -> no change.
	// Per-spring factors are written to lsedgesoft.txt (i, j, de, f) for checking.
	// ---------------------------------------------------------------------------
	bool esReady = false;

	void buildEdgeSoftening()
	{
		esReady = true;
		if(LS_STRETCH_EDGE_FACTOR == 1.0) return;
		std::map< std::vector<int>, int > triSeen;
		std::map< std::pair<int,int>, int > edgeCnt;
		for(int b = 0; b < c.ls.nBL; b++)
		{
			int n[3] = { c.ls.BL(b,0), c.ls.BL(b,1), c.ls.BL(b,2) };
			if(n[0] < 0 || n[1] < 0 || n[2] < 0) continue;
			int o = c.ls.O(n[0]);
			if(o < 1 || c.ls.O(n[1]) != o || c.ls.O(n[2]) != o) continue;
			std::vector<int> key(n, n + 3); std::sort(key.begin(), key.end());
			if(triSeen.count(key)) continue;
			triSeen[key] = 1;
			for(int e = 0; e < 3; e++)
				edgeCnt[std::pair<int,int>(std::min(n[e], n[(e + 1) % 3]), std::max(n[e], n[(e + 1) % 3]))]++;
		}
		std::vector<int> segA, segB;                                     // free-edge segments
		for(std::map< std::pair<int,int>, int >::iterator it = edgeCnt.begin(); it != edgeCnt.end(); ++it)
		{
			if(it->second != 1) continue;
			int a = it->first.first, bb = it->first.second;
			bool ancA = (c.ls.T(a) & (LS_CUSTOM_NODE_0 | LS_CUSTOM_NODE_1)) != 0;
			bool ancB = (c.ls.T(bb) & (LS_CUSTOM_NODE_0 | LS_CUSTOM_NODE_1)) != 0;
			if(ancA && ancB) continue;
			segA.push_back(a); segB.push_back(bb);
		}
		FILE *fo = fopen("lsedgesoft.txt", "w");
		double fMin = 1., fSum = 0.; int ns = 0;
		for(int k0 = 0; k0 < c.ls.nN; k0++)
		{
			int o = c.ls.O(k0);
			if(o < 1 || c.ls.N.NL(k0) == 0) continue;
			for(int jk = 0; jk < c.ls.nL; jk++)
			{
				int k1 = c.ls.N(k0,jk);
				if(k1 < 0 || k1 == k0 || c.ls.O(k1) != o) continue;
				clVector3Dd mid = (vls.C.c0(k0) + vls.C.c0(k1)) * 0.5;
				double de = 1.e30;
				for(size_t s = 0; s < segA.size(); s++)
				{
					if(c.ls.O(segA[s]) != o) continue;
					clVector3Dd a = vls.C.c0(segA[s]), ab = vls.C.c0(segB[s]) - a;
					double t = ((mid - a) * ab) / (ab * ab);
					if(t < 0.) t = 0.;
					if(t > 1.) t = 1.;
					double dd = (a + ab * t - mid).abs();
					if(dd < de) de = dd;
				}
				double ramp = (LS_STRETCH_EDGE_RAMP > 0.) ? de / LS_STRETCH_EDGE_RAMP : 1.;
				if(ramp > 1.) ramp = 1.;
				double f = LS_STRETCH_EDGE_FACTOR + (1. - LS_STRETCH_EDGE_FACTOR) * ramp;
				c.ls.K(k0,jk) *= f;                                          // both rows of the spring are scaled
				if(k0 < k1)
				{
					if(f < fMin) fMin = f;
					fSum += f; ns++;
					if(fo) fprintf(fo, "%d %d %.4f %.4f\n", k0, k1, de, f);
				}
			}
		}
		if(fo) fclose(fo);
		printf("\nLS_STRETCH_EDGE: %d free-edge segments, %d leaflet springs softened x%.3g at the free edge over %.3g lu -> factor min %.3f, mean %.3f\n",
			(int)segA.size(), ns, (double)LS_STRETCH_EDGE_FACTOR, (double)LS_STRETCH_EDGE_RAMP, fMin, ns ? fSum / ns : 1.);
	}

	void solveLSspring()
	{
		if(!esReady) buildEdgeSoftening();
		double Smax = 0., Smin = 0.;
		int i0;
#ifdef LS_USE_OPENMP
		// Gather form: every node sums the forces of all its own springs, so no two
		// threads ever write the same node (the old OpenMP branch here had its force
		// accumulation commented out, i.e. it silently dropped all spring forces).
		// Each spring is evaluated from its lower-index end's row, exactly as the
		// serial scatter loop below does, because v2 stores stiffness per node row.
		#pragma omp parallel default(shared)
		{
			double Smax_t = 0., Smin_t = 0.;
			int k0;
			#pragma omp for
			for(k0 = 0; k0 < c.ls.nN; k0++)
			{
				if(c.ls.N.NL(k0) == 0)
					continue;
				BOOL Tno_fk = c.ls.Tno_f(k0);
				clVector3Dd Fsum(0., 0., 0.);
				for(int jk = 0; jk < c.ls.nL; jk++)
				{
					int k1 = c.ls.N(k0,jk);
					if(k1 < 0 || k1 == k0 || (Tno_fk && c.ls.Tno_f(k1)))
						continue;
					int ra = k0, rb = jk;						// row that owns this spring
					if(k1 < k0) { ra = k1; rb = c.ls.N.find(k1, k0); }
					clVector3Dd dC = vls.C(k0) - vls.C(k1);
					double L = dC.abs(), L0 = c.ls.L0(ra,rb);
					clVector3Dd nn = dC/L;
					double dL = L - L0, S = dL / L0;
					vls.C.L(k0,jk) = L;
					vls.C.S(k0,jk) = S;
					if(S > Smax_t) Smax_t = S;
					else if(S < Smin_t) Smin_t = S;
					double K = c.ls.K(ra,rb);
					c.ls.K.E(k0,jk) = 0.5 * K * dL * dL;
					double Fd = 0.;
#ifdef LS_DISSIPATION
#ifndef LS_DISSIPATION_VELOCITY
					clVector3Dd dV = vls.V(k0) - vls.V(k1);
					Fd = c.ls.Ds(ra,rb) * (dV * nn);
#endif //LS_DISSIPATION_VELOCITY
#endif //LS_DISSIPATION
					Fsum += nn * (-(K * dL + Fd));
#ifdef LS_DISSIPATION_VELOCITY
					Fsum -= vls.V(k0) * c.ls.Ds(ra,rb);
#endif //LS_DISSIPATION_VELOCITY
				}
				vls.F.s(k0) += Fsum;
			}
			#pragma omp critical
			{
				if(Smax_t > Smax) Smax = Smax_t;
				if(Smin_t < Smin) Smin = Smin_t;
			}
		}
#else //LS_USE_OPENMP
		for(i0 = 0; i0 < c.ls.nN; i0++)
		{
			if(c.ls.N.NL(i0) == 0)
				continue;

			BOOL Tno_f0 = c.ls.Tno_f(i0);

			for(int j0 = 0; j0 < c.ls.nL; j0++)
			{
				int i1 = c.ls.N(i0,j0);
				if(i1 < i0 || (Tno_f0 && c.ls.Tno_f(i1)))
					continue;

				clVector3Dd dC = vls.C(i0) - vls.C(i1);
				double L = dC.abs(), L0 = c.ls.L0(i0,j0);
				clVector3Dd nn = dC/L;
				int j1 = c.ls.N.find(i1, i0);
				double dL = L - L0, S = dL / L0;
				vls.C.L(i0,j0) = L;
				vls.C.L(i1,j1) = L;
				vls.C.S(i0,j0) = S;
				vls.C.S(i1,j1) = S;

//#ifdef LS_USE_OPENMP
//				#pragma omp critical
//#endif //LS_USE_OPENMP
				{
					if(S > Smax)
						Smax = S;
					else if(S < Smin)
						Smin = S;
				}
				
				double K = c.ls.K(i0,j0);
				double Ek = 0.5 * K * dL * dL;
				c.ls.K.E(i0,j0) = Ek;
				c.ls.K.E(i1,j1) = Ek;

				double F, Fd = 0., Fk;
				Fk = K * (L - L0);
#ifdef LS_DISSIPATION
#ifndef LS_DISSIPATION_VELOCITY
				clVector3Dd dV = vls.V(i0) - vls.V(i1);
				Fd = c.ls.Ds(i0,j0) * (dV * nn);
#endif //LS_DISSIPATION_VELOCITY
#endif //LS_DISSIPATION
				F = -(Fk + Fd);

				clVector3Dd Ftot = nn * F;
#ifdef LS_USE_OPENMP
//#pragma omp atomic
				//vls.F(i0).x += Ftot.x;
//#pragma omp atomic
				//vls.F(i0).y += Ftot.y;
//#pragma omp atomic
				//vls.F(i0).z += Ftot.z;
//#pragma omp atomic
				//vls.F(i1).x -= Ftot.x;
//#pragma omp atomic
				//vls.F(i1).y -= Ftot.y;
//#pragma omp atomic
				//vls.F(i1).z -= Ftot.z;
#else 
				vls.F.s(i0) += Ftot;
				vls.F.s(i1) -= Ftot;
#ifdef LS_DISSIPATION_VELOCITY
				vls.F.s(i0) -= vls.V(i0) * c.ls.Ds(i0,j0);
				vls.F.s(i1) -= vls.V(i1) * c.ls.Ds(i0,j0);
#endif //LS_DISSIPATION_VELOCITY
#endif //LS_USE_OPENMP
			}
		}
#endif //LS_USE_OPENMP

		vls.C.Smax = Smax;
		vls.C.Smin = Smin;
	}

#ifdef LS_HINGE_BENDING
#ifndef LS_HINGE_WALL_FACTOR
#define LS_HINGE_WALL_FACTOR	1.0
#endif
#ifndef LS_HINGE_WALL_RAMP
#define LS_HINGE_WALL_RAMP	3.0
#endif
	// ---------------------------------------------------------------------------
	// Dihedral-hinge bending for the leaflets (2026-09-26).
	// Every interior edge shared by two leaflet triangles is a hinge (i, j, k, l):
	// j-k is the shared edge, i and l are the tips of the two triangles. theta is the
	// signed dihedral angle between the two triangle planes and
	//     E = 0.5 * kh * (theta - theta0)^2 ,
	// theta0 taken from the reference (undeformed) mesh. For a triangulated sheet the
	// continuum bending rigidity is D_b = (sqrt(3)/2) kh (Seung & Nelson 1988), and the
	// straight-triplet model it replaces had D_b = (3 sqrt(3)/4) k_b, so kh = 1.5 k_b keeps
	// the same D_b and the same meaning of K; k_b is LS_ANGLE_STIFFNESS_1/2 as before.
	// Checked in Template_v2/HingeBendingCheck.py: gradient vs finite differences (1e-9)
	// and, on the real leaflet mesh bent to a cylinder, hinge energy / plate energy 0.94-0.99.
	// ---------------------------------------------------------------------------
	std::vector<int> hgN;          // 4 node indices per hinge: i, j, k, l
	std::vector<double> hgTh0;     // rest dihedral angle
	std::vector<double> hgK;       // hinge stiffness kh
	bool hgReady = false;

	static double hgAngle(clVector3Dd xi, clVector3Dd xj, clVector3Dd xk, clVector3Dd xl)
	{
		clVector3Dd b1 = xj - xi, b2 = xk - xj, b3 = xl - xk;
		clVector3Dd m = b1 ^ b2, n = b2 ^ b3;
		return atan2(b2.abs() * (b1 * n), m * n);
	}

	void buildHinges()
	{
		hgN.clear(); hgTh0.clear(); hgK.clear();
		std::map< std::vector<int>, int > triSeen;                 // dedupe the two windings
		std::map< std::pair<int,int>, std::vector<int> > edgeOpp;  // edge -> opposite vertices
		for(int b = 0; b < c.ls.nBL; b++)
		{
			int n[3] = { c.ls.BL(b,0), c.ls.BL(b,1), c.ls.BL(b,2) };
			if(n[0] < 0 || n[1] < 0 || n[2] < 0) continue;
			int o = c.ls.O(n[0]);
			if(o < 1 || c.ls.O(n[1]) != o || c.ls.O(n[2]) != o) continue;      // leaflet triangles only
			std::vector<int> key(n, n + 3); std::sort(key.begin(), key.end());
			if(triSeen.count(key)) continue;
			triSeen[key] = 1;
			for(int e = 0; e < 3; e++)
			{
				int a = n[e], bb = n[(e + 1) % 3], opp = n[(e + 2) % 3];
				std::pair<int,int> ed(std::min(a, bb), std::max(a, bb));
				edgeOpp[ed].push_back(opp);
			}
		}
		// Optional softening near the vessel wall (root / commissures), 2026-09-26:
		//   kh -> kh * f,  f = LS_HINGE_WALL_FACTOR + (1 - LS_HINGE_WALL_FACTOR) * min(dw / LS_HINGE_WALL_RAMP, 1)
		// dw = distance of the hinge edge's midpoint from the wall in the reference mesh.
		// LS_HINGE_WALL_FACTOR = 1 (default) -> uniform stiffness, i.e. no softening.
		double hgYc = (c.lb.nY - 1.0) / 2.0, hgZc = (c.lb.nZ - 1.0) / 2.0, hgRw = 0.; int hgNw = 0;
		for(int q = 0; q < c.ls.nN; q++)
			if(c.ls.O(q) == 0)
			{
				double dy = vls.C.c0(q).y - hgYc, dz = vls.C.c0(q).z - hgZc;
				hgRw += sqrt(dy*dy + dz*dz); hgNw++;
			}
		hgRw = (hgNw > 0) ? hgRw / hgNw : 1.e30;
		double fMin = 1., fSum = 0.;
		for(std::map< std::pair<int,int>, std::vector<int> >::iterator it = edgeOpp.begin(); it != edgeOpp.end(); ++it)
		{
			if(it->second.size() != 2) continue;                          // boundary edge: no hinge
			int i = it->second[0], j = it->first.first, k = it->first.second, l = it->second[1];
			hgN.push_back(i); hgN.push_back(j); hgN.push_back(k); hgN.push_back(l);
			hgTh0.push_back(hgAngle(vls.C.c0(i), vls.C.c0(j), vls.C.c0(k), vls.C.c0(l)));
			double kb = (c.ls.O(j) == 1) ? c.ls.AKini1 : c.ls.AKini2;
			clVector3Dd mid = (vls.C.c0(j) + vls.C.c0(k)) * 0.5;
			double dy = mid.y - hgYc, dz = mid.z - hgZc;
			double dw = hgRw - sqrt(dy*dy + dz*dz);
			double ramp = (LS_HINGE_WALL_RAMP > 0.) ? dw / LS_HINGE_WALL_RAMP : 1.;
			if(ramp < 0.) ramp = 0.;
			if(ramp > 1.) ramp = 1.;
			double fw = LS_HINGE_WALL_FACTOR + (1. - LS_HINGE_WALL_FACTOR) * ramp;
			if(fw < fMin) fMin = fw;
			fSum += fw;
			hgK.push_back(1.5 * kb * fw);
		}
		hgReady = true;
		printf("\nLS_HINGE_BENDING: %d leaflet hinges built (kh = 1.5 * LS_ANGLE_STIFFNESS); wall softening %.3g over %.3g lu -> factor min %.3f, mean %.3f\n",
			(int)hgTh0.size(), (double)LS_HINGE_WALL_FACTOR, (double)LS_HINGE_WALL_RAMP, fMin, hgTh0.size() ? fSum / hgTh0.size() : 1.);
	}

	void solveLShinges()
	{
		if(!hgReady) buildHinges();
		int nh = (int)hgTh0.size();
		for(int h = 0; h < nh; h++)
		{
			int i = hgN[4*h], j = hgN[4*h+1], k = hgN[4*h+2], l = hgN[4*h+3];
			if(c.ls.Tno_f(i) == TRUE && c.ls.Tno_f(j) == TRUE && c.ls.Tno_f(k) == TRUE && c.ls.Tno_f(l) == TRUE)
				continue;
			clVector3Dd xi = vls.C(i), xj = vls.C(j), xk = vls.C(k), xl = vls.C(l);
			clVector3Dd b1 = xj - xi, b2 = xk - xj, b3 = xl - xk;
			clVector3Dd m = b1 ^ b2, n = b2 ^ b3;
			double mm = m * m, nn = n * n, b2l = b2.abs();
			if(mm < 1.e-20 || nn < 1.e-20 || b2l < 1.e-12) continue;      // degenerate hinge
			double th = atan2(b2l * (b1 * n), m * n);
			double d = th - hgTh0[h];
			while(d >  PI_NUMBER) d -= 2.0 * PI_NUMBER;
			while(d < -PI_NUMBER) d += 2.0 * PI_NUMBER;
			// d(theta)/dx for the four nodes (b1 = xj - xi convention)
			clVector3Dd gi = m * (-b2l / mm);
			clVector3Dd gl = n * ( b2l / nn);
			double p = (b1 * b2) / (b2l * b2l), q = (b3 * b2) / (b2l * b2l);
			clVector3Dd gj = gi * (-(1.0 + p)) + gl * q;
			clVector3Dd gk = gl * (-(1.0 + q)) + gi * p;
			double f = -hgK[h] * d;                                        // F = -dE/dx = -kh (theta - theta0) dtheta/dx
			vls.F.s(i) += gi * f;
			vls.F.s(j) += gj * f;
			vls.F.s(k) += gk * f;
			vls.F.s(l) += gl * f;
		}
	}
#endif //LS_HINGE_BENDING

#ifndef LS_COMMISSURE_TIE_K
#define LS_COMMISSURE_TIE_K	0.0
#endif
#ifndef LS_COMMISSURE_TIE_D
#define LS_COMMISSURE_TIE_D	5.0
#endif
#ifndef LS_COMMISSURE_TIE_OPEN
#define LS_COMMISSURE_TIE_OPEN	0.0
#endif
	// ---------------------------------------------------------------------------
	// Commissure ties (2026-09-26, Amir's request): the closed valves kept two thin slits next to
	// the commissures, where the free edges stop ~0.7 lu short of the mid-plane. Real leaflets are
	// joined at the commissures, so the outermost UNCLAMPED free-edge node of each leaflet is tied to
	// its mirror node on the opposite leaflet of the same valve (1-2 inlet, 3-4 outlet) by a zero-length
	// spring with damping on the relative velocity:
	//     F_a = K (x_b - x_a) + D (v_b - v_a),  F_b = -F_a
	// Free edge = leaflet boundary edges, nodes outside the clamp band (LS_CUSTOM_NODE_0/1) - the same
	// definition as buildEdgeSoftening / CompareOldData.py. K = LS_COMMISSURE_TIE_K (0 = off, default).
	// Closing-only ties (LS_COMMISSURE_TIE_OPEN > 0, 2026-09-26): a permanent tie pinched the OPEN orifice
	// (outlet open area -30%), so each tie is scaled by w = clamp((OPEN - dc) / (OPEN/2), 0, 1), where dc is
	// the current distance between the two CENTRAL free-edge nodes of that valve: full strength once the valve
	// centre has shut (dc <= OPEN/2), zero as soon as it opens past OPEN. OPEN = 0 -> always on.
	// ---------------------------------------------------------------------------
	std::vector<int> ctA, ctB, ctCa, ctCb;       // tie pair, and that valve's central free-edge node pair
	bool ctReady = false;

	void buildCommissureTies()
	{
		ctReady = true;
		if(LS_COMMISSURE_TIE_K == 0.0) return;
		std::map< std::vector<int>, int > triSeen;
		std::map< std::pair<int,int>, int > edgeCnt;
		for(int b = 0; b < c.ls.nBL; b++)
		{
			int n[3] = { c.ls.BL(b,0), c.ls.BL(b,1), c.ls.BL(b,2) };
			if(n[0] < 0 || n[1] < 0 || n[2] < 0) continue;
			int o = c.ls.O(n[0]);
			if(o < 1 || c.ls.O(n[1]) != o || c.ls.O(n[2]) != o) continue;
			std::vector<int> key(n, n + 3); std::sort(key.begin(), key.end());
			if(triSeen.count(key)) continue;
			triSeen[key] = 1;
			for(int e = 0; e < 3; e++)
				edgeCnt[std::pair<int,int>(std::min(n[e], n[(e + 1) % 3]), std::max(n[e], n[(e + 1) % 3]))]++;
		}
		std::vector<char> isFree(c.ls.nN, 0);
		for(std::map< std::pair<int,int>, int >::iterator it = edgeCnt.begin(); it != edgeCnt.end(); ++it)
		{
			if(it->second != 1) continue;
			int ab[2] = { it->first.first, it->first.second };
			for(int q = 0; q < 2; q++)
				if((c.ls.T(ab[q]) & (LS_CUSTOM_NODE_0 | LS_CUSTOM_NODE_1)) == 0) isFree[ab[q]] = 1;
		}
		double yc = (c.lb.nY - 1.0) / 2.0;
		int pairs[2][2] = { {1, 2}, {3, 4} };
		for(int p = 0; p < 2; p++)
			for(int s = -1; s <= 1; s += 2)
			{
				int best[2] = { -1, -1 };
				for(int q = 0; q < 2; q++)
				{
					double bmax = -1.e30;
					for(int i = 0; i < c.ls.nN; i++)
						if(isFree[i] && c.ls.O(i) == pairs[p][q] && s * (vls.C.c0(i).y - yc) > bmax)
						{ bmax = s * (vls.C.c0(i).y - yc); best[q] = i; }
				}
				if(best[0] < 0 || best[1] < 0) continue;
				int cen[2] = { -1, -1 };
				for(int q = 0; q < 2; q++)
				{
					double bmin = 1.e30;
					for(int i = 0; i < c.ls.nN; i++)
						if(isFree[i] && c.ls.O(i) == pairs[p][q] && fabs(vls.C.c0(i).y - yc) < bmin)
						{ bmin = fabs(vls.C.c0(i).y - yc); cen[q] = i; }
				}
				ctA.push_back(best[0]); ctB.push_back(best[1]); ctCa.push_back(cen[0]); ctCb.push_back(cen[1]);
				clVector3Dd d0 = vls.C.c0(best[1]) - vls.C.c0(best[0]);
				printf("\nLS_COMMISSURE_TIE: leaflets %d-%d side %+d: nodes %d-%d at y-yc %+.2f, rest separation %.3f",
					pairs[p][0], pairs[p][1], s, best[0], best[1], vls.C.c0(best[0]).y - yc, d0.abs());
			}
		printf("\nLS_COMMISSURE_TIE: %d ties, K = %g, D = %g, closing-only threshold %g lu (0 = always on)\n", (int)ctA.size(),
			(double)LS_COMMISSURE_TIE_K, (double)LS_COMMISSURE_TIE_D, (double)LS_COMMISSURE_TIE_OPEN);
	}

	void solveLSties()
	{
		if(!ctReady) buildCommissureTies();
		for(size_t t = 0; t < ctA.size(); t++)
		{
			int a = ctA[t], b = ctB[t];
			if(c.ls.Tno_f(a) == TRUE && c.ls.Tno_f(b) == TRUE) continue;
			double w = 1.;
			if(LS_COMMISSURE_TIE_OPEN > 0.)
			{
				double dc = (vls.C(ctCb[t]) - vls.C(ctCa[t])).abs();
				w = (LS_COMMISSURE_TIE_OPEN - dc) / (0.5 * LS_COMMISSURE_TIE_OPEN);
				if(w <= 0.) continue;
				if(w > 1.) w = 1.;
			}
			clVector3Dd F = ((vls.C(b) - vls.C(a)) * LS_COMMISSURE_TIE_K + (vls.V(b) - vls.V(a)) * LS_COMMISSURE_TIE_D) * w;
			vls.F.s(a) += F;
			vls.F.s(b) -= F;
		}
	}

	void solveLSangle()
	{
		int i;
#ifdef LS_USE_OPENMP
		#pragma omp parallel for default(shared) private(i)
#endif //LS_USE_OPENMP
		for(i = 0; i < c.ls.nAN; i++)
		{
			int i0 = c.ls.AN(i,0), i1 = c.ls.AN(i,1), i2 = c.ls.AN(i,2);

			if(c.ls.Tno_f(i0) == TRUE && c.ls.Tno_f(i1) == TRUE && c.ls.Tno_f(i2) == TRUE)
				continue;
#ifdef LS_HINGE_BENDING
			if(c.ls.O(i1) >= 1)          // leaflets use dihedral hinges instead (solveLShinges)
				continue;
#endif //LS_HINGE_BENDING

			clVector3Dd d1 = vls.C(i0) - vls.C(i1);
			clVector3Dd d2 = vls.C(i2) - vls.C(i1);
			double r1 = d1.abs(), r2 = d2.abs();
			double cs = d1 * d2 / (r1 * r2); //cos
			if(cs > 1.) 
				cs = 1.;
			else if(cs < -1.) 
				cs = -1.;

			//double sn = sqrt(1. - cs * cs);
			double k = c.ls.AK(i);

			c.ls.AK.E(i) += k * (1. + cs);
        
			clVector3Dd f0 = -(d2 / (r1 * r2) - d1 * cs / (r1 * r1)) * k;
			clVector3Dd f2 = -(d1 / (r1 * r2) - d2 * cs / (r2 * r2)) * k;
			clVector3Dd f1 = -(f0 + f2);

#ifdef LS_USE_OPENMP
			#pragma omp atomic
			vls.F.s(i0).x +=  f0.x;
			#pragma omp atomic
			vls.F.s(i0).y +=  f0.y;
			#pragma omp atomic
			vls.F.s(i0).z +=  f0.z;
			#pragma omp atomic
			vls.F.s(i1).x +=  f1.x;
			#pragma omp atomic
			vls.F.s(i1).y +=  f1.y;
			#pragma omp atomic
			vls.F.s(i1).z +=  f1.z;
			#pragma omp atomic
			vls.F.s(i2).x +=  f2.x;
			#pragma omp atomic
			vls.F.s(i2).y +=  f2.y;
			#pragma omp atomic
			vls.F.s(i2).z +=  f2.z;
#else
			vls.F.s(i0) +=  f0;
			vls.F.s(i1) +=  f1;
			vls.F.s(i2) +=  f2;
#endif //LS_USE_OPENMP
		}
#ifdef LS_HINGE_BENDING
		solveLShinges();
#endif //LS_HINGE_BENDING
	}

	void solveRigidBody()
	{
		if(c.ls.nTr == 0)
			return;

		vls.R.update();
		//return;

		for(int i = 0; i < c.ls.nTr; i++)
		{
			int n0 = c.ls.Tr(i);
			int o = c.ls.O(n0);

			double M = c.ls.M(n0);
			clVector3Dd R = vls.C(n0) - vls.R.CG(o);
			//clVector3Dd Rt = vls.R.On(o) * (R * vls.R.On(o));
			//clVector3Dd Rn = R - Rt;

			clVector3Dd O = vls.R.O(o);
			//if (vls.R.CG(o).z < 1.01 * (R0_NUMBER - 0.5))	//rolling, in contact with wall
			//{
			//		O = clVector3Dd(0., 0., 1) ^ vls.R.V(o) / R0_NUMBER;
			//}
			clVector3Dd Fcentrefug = O ^ (O ^ R) * M;
			//clVector3Dd Vr = vls.V(n0) - vls.R.V(o);
			//clVector3Dd Vo = O ^ R;
			//clVector3Dd Ucoriolis = Vr - Vo;
			//clVector3Dd Fcoriolis = O ^ Ucoriolis * 2. * M;

			clVector3Dd Odot = vls.R.II1(o) * vls.R.T(o);
			clVector3Dd AO = Odot ^ R;

			vls.F(n0) = vls.R.F(o) *M/c.ls.M.O(o) + Fcentrefug + AO * M;

///////////////////////////////////////////////////////////////////////////////////
////				vls.F(i).x += (FA.x + MA * C.y + Fomega * C.n().x) * Mi;
////				vls.F(i).y += (FA.y - MA * C.x + Fomega * C.n().y) * Mi;
///////////////////////////////////////////////////////////////////////////////////
		}
	}

	void solveRigidBodyCorrection()
	{
		for(int i = 0; i < c.ls.nTr; i++)
		{
			double damp_mu = LS_DAMPING_CONSTANT;
			int n0 = c.ls.Tr(i);
			for(int j = 0; j < c.ls.nL; j++)
			{
				int n = c.ls.N(n0,j);
				if(n < 0 || !(c.ls.T(n) & LS_RIGID_NODE))
					continue;
				clVector3Dd d = vls.C(n0) - vls.C(n);
				double L = d.abs();
				double F = -c.ls.K(n0,j) * (L - c.ls.L0(n0,j));
				vls.F(n0) +=  d * F / L;
				// Damping
				clVector3Dd dv = vls.V(n0) - vls.V(n);
				double dvdotd = dv * d;
				vls.F(n0) -= d * damp_mu * dvdotd / (L * L);
			}
		}
		for(int i = 0; i < c.ls.nAN; i++)
		{
			int i0 = c.ls.AN(i,0), i1 = c.ls.AN(i,1), i2 = c.ls.AN(i,2);

			if(c.ls.Tno_f(i1) == TRUE || !(c.ls.T(i1) & LS_RIGID_NODE))
				continue;

			clVector3Dd d1 = vls.C(i0) - vls.C(i1);
			clVector3Dd d2 = vls.C(i2) - vls.C(i1);
			double r1 = d1.abs(), r2 = d2.abs();
			double cs = d1 * d2 / (r1 * r2); //cos
			if(cs > 1.) 
				cs = 1.;
			else if(cs < -1.) 
				cs = -1.;

			double sn = sqrt(1. - cs * cs);
			if (sn < c.eps) sn = c.eps;
			sn = 1. / sn;

			//double dcos = cs - c.ls.AA0(i);
			double dx = acos(cs) - c.ls.AA0(i);
			double kdx = dx * c.ls.AK(i) * sn;
			double en = dx * dx * c.ls.AK(i);
	
			double a11 = kdx * cs / (r1*r1), a12 = -kdx / (r1*r2), a22 = kdx * cs / (r2*r2);
        
			clVector3Dd f1 = d1 * a11 + d2 * a12;
			clVector3Dd f2 = d1 * a12 + d2 * a22;

			vls.F(i0) -=  f1;
			vls.F(i1) +=  f1 + f2;
			vls.F(i2) -=  f2;
		}
	}
	
	
/////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////LB solver//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
	void solveLB()
	{
		bcLB();
		propagationLB();
		collisionLB();
#ifdef LB_COARSE
		vlb.C.update();
#endif //LB_COARSE
	}

	void bcLB()
	{
		vlb.F.copy();
#ifdef BINARY_FLUID
		vlb.G.copy();
#endif
		vlb.M.update();

#ifdef LB_COARSE
#pragma omp parallel sections
		{
#pragma omp section
		{if(c.lb.bcCoarseX0 == TRUE) bcCoarseLBx0();}
#pragma omp section
		{if(c.lb.bcCoarseX1 == TRUE) bcCoarseLBx1();}
#pragma omp section
		{if(c.lb.bcCoarseY0 == TRUE) bcCoarseLBy0();}
#pragma omp section
		{if(c.lb.bcCoarseY1 == TRUE) bcCoarseLBy1();}
#pragma omp section
		{if(c.lb.bcCoarseZ0 == TRUE) bcCoarseLBz0();}
#pragma omp section
		{if(c.lb.bcCoarseZ1 == TRUE) bcCoarseLBz1();}
		}
#endif //LB_COARSE

		vlb.J.updateinvel();
		bcSetBoundaryLS(); //if there is a problem with bc (periodic) try to move the funciton down

		if(c.lb.bcPeriodicX == TRUE && (c.lb.bcCoarseX0 == FALSE || c.lb.bcCoarseX1 == FALSE)) bcPeriodicLBx();
		if(c.lb.bcPeriodicY == TRUE && (c.lb.bcCoarseY0 == FALSE || c.lb.bcCoarseY1 == FALSE)) bcPeriodicLBy();
		if(c.lb.bcPeriodicZ == TRUE && (c.lb.bcCoarseZ0 == FALSE || c.lb.bcCoarseZ1 == FALSE)) bcPeriodicLBz();

		if(c.lb.bcFreeX0 == TRUE && c.lb.bcCoarseX0 == FALSE) bcFreeFlowLBx0();
		if(c.lb.bcFreeX1 == TRUE && c.lb.bcCoarseX1 == FALSE) bcFreeFlowLBx1();
		if(c.lb.bcFreeY0 == TRUE && c.lb.bcCoarseY0 == FALSE) bcFreeFlowLBy0();
		if(c.lb.bcFreeY1 == TRUE && c.lb.bcCoarseY1 == FALSE) bcFreeFlowLBy1();
		if(c.lb.bcFreeZ0 == TRUE && c.lb.bcCoarseZ0 == FALSE) bcFreeFlowLBz0();
		if(c.lb.bcFreeZ1 == TRUE && c.lb.bcCoarseZ1 == FALSE) bcFreeFlowLBz1();

		if(c.lb.bcWallX0 == TRUE && c.lb.bcCoarseX0 == FALSE) bcWallLBx0();
		if(c.lb.bcWallX1 == TRUE && c.lb.bcCoarseX1 == FALSE) bcWallLBx1();
		if(c.lb.bcWallY0 == TRUE && c.lb.bcCoarseY0 == FALSE) bcWallLBy0();
		if(c.lb.bcWallY1 == TRUE && c.lb.bcCoarseY1 == FALSE) bcWallLBy1();
		if(c.lb.bcWallZ0 == TRUE && c.lb.bcCoarseZ0 == FALSE) bcWallLBz0();
		if(c.lb.bcWallZ1 == TRUE && c.lb.bcCoarseZ1 == FALSE) bcWallLBz1();

		if(c.lb.bcSymmetryX0 == TRUE && c.lb.bcCoarseX0 == FALSE) bcSymmetryLBx0();
		if(c.lb.bcSymmetryX1 == TRUE && c.lb.bcCoarseX1 == FALSE) bcSymmetryLBx1();
		if(c.lb.bcSymmetryY0 == TRUE && c.lb.bcCoarseY0 == FALSE) bcSymmetryLBy0();
		if(c.lb.bcSymmetryY1 == TRUE && c.lb.bcCoarseY1 == FALSE) bcSymmetryLBy1();
		if(c.lb.bcSymmetryZ0 == TRUE && c.lb.bcCoarseZ0 == FALSE) bcSymmetryLBz0();
		if(c.lb.bcSymmetryZ1 == TRUE && c.lb.bcCoarseZ1 == FALSE) bcSymmetryLBz1();

		if(c.lb.bcPressX0 == TRUE && c.lb.bcCoarseX0 == FALSE) bcPressLBx0();		// Keep this after bounce-back BCs
		if(c.lb.bcPressX1 == TRUE && c.lb.bcCoarseX1 == FALSE) bcPressLBx1();
		if(c.lb.bcPressY0 == TRUE && c.lb.bcCoarseY0 == FALSE) bcPressLBy0();
		if(c.lb.bcPressY1 == TRUE && c.lb.bcCoarseY1 == FALSE) bcPressLBy1();
		if(c.lb.bcPressZ0 == TRUE && c.lb.bcCoarseZ0 == FALSE) bcPressLBz0();
		if(c.lb.bcPressZ1 == TRUE && c.lb.bcCoarseZ1 == FALSE) bcPressLBz1();


#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
		vlb.B.save2file();
		vlb.M.save2file();
		vlb.Ro.save2file();
		vlb.J.save2file();
#endif //LB_DEBUG

		vlb.M.test();
		vlb.F.change();
#ifdef BINARY_FLUID
		vlb.G.change();
#endif


#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
		vlb.B.save2file();
#endif //LB_DEBUG
	}

	void bcSetBoundaryLS()
	{
		vls.F.zerosB();
		c.ls.BLFb.zeros();
		bcFindBoundaryNodes();

//		vls.F.smooth();

	}

	void bcSetNewF(clVector3Di ii0, clVector3Dd UU, const char map, const int no, clVector3Dd NN)
	{
		return;			// Uncomment to ignore this section (Wenbin said it is inaccurate, and it also seems to have made code not run on cluster)
		//char map0 = vlb.M(ii0);
		double Ro = c.lb.Ro[map];
		double countRo = 0., sumRo = 0.;
		double N = 0., Nd2 = 0.; //set constants for defult values
		double sumN = 0., sumNd2 = 0.;
/////////==========================================================
		//if(map == 0)
		//	printf("\nmap = 0 in bcSetNewF\n");
/////////==========================================================
#ifdef BINARY_FLUID
		for(int i = 1;  i < c.lb.nL; i++)
		{
			if(NN * c.lb.C[i] < 0.2)
				continue;

			clVector3Di ii = MOD(ii0 + c.lb.C[i], c.lb.nn);
			if(vlb.M(ii) == map)
			{
				countRo++;
//				sumRo += vlb.Ro(ii);
				sumN += vlb.N(ii);
				sumNd2 += vlb.N.d2(ii).abs();
			}
		}
//////============================================================
		if(countRo > 0.)
		{
//			Ro = sumRo / countRo;
			N = sumN / countRo;
			Nd2 = sumNd2 / countRo;
		}

//		Ro = vlb.Ro(ii0);
#endif

#ifdef LB_USE_LBGK
		double Gtot = 0., Ftot = 0.;
		for(int i = 1; i < c.lb.nL; i++)
		{
			clVector3Dd UUc = UU - vlb.J.cor(no, map) / Ro;
			double U = (c.lb.C[i] * UUc) / c.lb.Cs2;
//			double U = (c.lb.C[i] * UU) / c.lb.Cs2;
			double F = Ro * c.lb.F0[i] * (1. + U) / c.lb.Ro[map];
			vlb.F(i)(ii0) = F;
			Ftot += F;

#ifdef BINARY_FLUID
			double GMu = LBB_GAMMA_NUMBER * ((-LBB_A + LBB_B * N * N) * N - LBB_KAPPA * Nd2) / c.lb.Cs2;
			double G = c.lb.F0[i] * (GMu + U * N);
			vlb.G(i)(ii0) = G;
			Gtot += G;
#endif

			clVector3Di ii1t = ii0 + c.lb.C[i], ii1 = MOD(ii1t, c.lb.nn);
			int BMflag = vlb.M.Bflag(ii1);
#ifdef LB_THIN_WALL
			if(!(BMflag & c.lb.LF[c.lb.rev[i]]) || testperiodicindex(ii1t, ii1) == FALSE)
			{
				vlb.F(i).nc(ii0) = F;
#ifdef BINARY_FLUID
				vlb.G(i).nc(ii0) = G;
#endif
			}
#endif //LB_THIN_WALL
		}
#endif //LB_USE_MRTLB
#ifdef LB_USE_MRTLB
		clVector3Dd  J = UU * Ro;

		double Meq[LB_NUMBER_CONNECTIONS];
		for(int i = 0; i < c.lb.nL; i++)
			Meq[i] = vlb.Meq(i, Ro, J);

		for(int i = 0; i < c.lb.nL; i++)
		{
			double F = 0.;
			for(int n = 0; n < c.lb.nL; n++)
				F += c.lb.M1[i][n] * Meq[n];

			vlb.F(i)(ii0) = F /  c.lb.Ro[map];

			clVector3Di ii1t = ii0 + c.lb.C[i], ii1 = MOD(ii1t, c.lb.nn);
			int BMflag = vlb.M.Bflag(ii1);
#ifdef LB_THIN_WALL
			if(!(BMflag & c.lb.LF[c.lb.rev[i]]) || testperiodicindex(ii1t, ii1) == FALSE)
				vlb.F(i).nc(ii0) = F;
#endif //LB_THIN_WALL
		}
#endif //LB_USE_MRTLB

		vlb.F(0)(ii0) = Ro / c.lb.Ro[map] - Ftot;
#ifdef BINARY_FLUID
		vlb.G(0)(ii0) = N - Gtot;
#endif
#ifdef LB_THIN_WALL
		vlb.F(0).nc(ii0) = Ro / c.lb.Ro[map] - Ftot;
#ifdef BINARY_FLUID
		vlb.G(0).nc(ii0) = N - Gtot;
#endif
#endif //LB_THIN_WALL

#ifdef BINARY_FLUID
		vlb.N.a(ii0) = N;
		vlb.N.updatexyz(ii0);
#endif

		vlb.M(ii0) = map;
		vlb.Ro.update(ii0);
		vlb.J.update(ii0);
		vlb.J.cor_add(vlb.J(ii0), no, map);
	}

	void bcSetNodeBC(clVector3Di ii0, int dir, double D, clVector3Dd UU,
									int bl, clVector3Dd vCc, clVector3Dd vC0, 
									clVector3Dd vC1, clVector3Dd vC2)
	{
		double Fn, F = vlb.F(dir)(ii0);
		clVector3Dd CC = c.lb.C[dir];
		char map = vlb.M(ii0);
		double Ro = c.lb.Ro[map];
		int dirp = c.lb.rev[dir];
		int inoutlet = c.ls.Inout(bl);
		int outletflag = 0;
		if (inoutlet > 1 && c.t.Time > 1)
		{
			UU = vlb.J(ii0) / vlb.Ro(ii0);//Ro;
#ifdef OUTLET_P_BC
			outletflag = inoutlet;
#endif	// OUTLET_P_BC
		}
		else if(inoutlet == 1)
		{
			//double osc = ::c.t.oscillator1();
			//UU.x = bcLungUniOsc();
			UU.x = vlb.J.invel;
			UU.y = 0.;
			UU.z = 0.;
#ifdef INLET_P_BC
			outletflag = inoutlet;
#endif // INLET_P_BC
		}

		double U = (CC * UU) * c.lb.F0[dir] / c.lb.Cs2  * vlb.Ro(ii0) / Ro;
#ifdef BINARY_FLUID
		double Gs = 0., Gn, G = vlb.G(dir)(ii0);
		double N = vlb.N(ii0);
		clVector3Dd Nd = vlb.N.d(ii0);
		double UUCC = (CC * UU);
		double UUUU = (UU * UU);
		double F0N = c.lb.F0[dir] * N;
#endif

		/////////////////////////////////////////////////////////
		/*clVector3Dd l01 = vC1 - vC0;
		clVector3Dd l02 = vC2 - vC0;
		clVector3Dd n = l01 ^ l02;
		if ( CC * n < 0)
			n = n / n.abs();
		else
			n = -n / n.abs();
		double temp = Nd * n;
		Nd = Nd - n * temp;*/

		//if(fabs(U) > 0.001)
		//	printf("\nError U bcSetNodeBC!\n");
		clVector3Di iis1t = ii0 - (clVector3Di)CC, iis1 = MOD(iis1t, c.lb.nn);

		if(map != vlb.M(iis1) || testperiodicindex(iis1t, iis1) == FALSE)
		{

//			Fn = F - 2. * U;
			Fn = c.lb.F0[dir] - U;
#ifdef BINARY_FLUID
#ifdef BINARY_SCALAR
			Gn = -G + 2 * F0N * (1 + 4.5 * UUCC*UUCC - 1.5 * UUUU);
#else
			Gn = G - 1. * U * (N + (CC * Nd) / 2.);////?
#endif // BINARY_SCALAR

#endif
		}
		else
		{
			double Fs1 = vlb.F(dir)(iis1);
			clVector3Di iis2t = ii0 - (clVector3Di)CC * 2, iis2 = MOD(iis2t, c.lb.nn);
//			if(D > 0. && testperiodicindex(iis2t, iis2) == TRUE)
//			{
//				double Fs2 = vlb.F(dir)(iis2);
//				double Fs1p = vlb.F(dirp)(ii0);
//				double Fs2p = vlb.F(dirp)(iis1);
//				double kt = 3. + 2. * D, ks = 1. + 6. * D + 4 * D * D, kf = 2 * (1 + D) * (1 + D), kst = 2 * ks / kf;
//				double k1 = kst * kf / ks - 1., k0 = 2. - kst, km1p = 2. - kst * kt / ks, km2p = kst * (kt - 2.) / ks - 1., km1 = k1 - km2p - 1., wq = 4. * kst / ks;
////				Fn = k1 * F + k0 * Fs1 + km1 * Fs2 + km1p * Fs1p + km2p * Fs2p - wq * U;
//				Fn = F + (1. - 2.*D - 2.*D*D) / (1 + D) / (1 + D) * (Fs1 - Fs1p) + D*D / (1 + D) / (1 + D) * (Fs2 - Fs2p) - 2. * U;
// 			}
//			else
//			{
#ifdef LB_HW_BB
			D = 0.5;
#endif // LB_HW_BB
			if(D < 0.5)
			{
				double Fs1 = vlb.F(dir)(iis1);
				//clVector3Di iis2t = ii0 - (clVector3Di)CC * 2, iis2 = MOD(iis2t, c.lb.nn);
//				if(testperiodicindex(iis2t, iis2) == FALSE)
					Fn = 2. * D * F + (1. - 2. * D) * Fs1 - 2. * U;
//				else
//				{
//					double Fs2 = vlb.F(dir)(iis2);
//					Fn = D * (2. * D + 1.) * F + (1. + 2. * D) * (1. - 2. * D) * Fs1 - D * (1. - 2. * D) * Fs2 - 2. * U;
//				}
#ifdef BINARY_FLUID
#ifdef BINARY_SCALAR
					Gn = -G + 2 * F0N * (1 + 4.5 * UUCC*UUCC - 1.5 * UUUU);
#else
				Gs = vlb.G(dir)(iis1);
				Gn = 2. * D * G + (1. - 2. * D) * Gs - 2. * U * (N + (CC * Nd) / 2.);
#endif // BINARY_SCALAR
#endif // BINARY_FLUID
			}
			else
			{
				double Fs1p = vlb.F(dirp)(ii0);
				clVector3Di iis1t = ii0 - (clVector3Di)CC, iis1 = MOD(iis1t, c.lb.nn);
//				if(testperiodicindex(iis1t, iis1) == FALSE)
					Fn = 1. / (2. * D) * F + (2. * D - 1.) / (2. * D) * Fs1p - 1./D * U;
//				else
//				{
//					double Fs2p = vlb.F(dirp)(iis1);
//					Fn = 1. / (D * (2. * D + 1.)) * F + (2. * D - 1.) / D * Fs1p + (1. - 2. * D) / (1. + 2. * D) * Fs2p - 2. / (D * (2. * D + 1.)) * U;
//				}
#ifdef BINARY_FLUID
#ifdef BINARY_SCALAR
					Gn = -G + 2 * F0N * (1 + 4.5 * UUCC*UUCC - 1.5 * UUUU);
#else
				Gs = vlb.G(dirp)(ii0);
				Gn = 1. / (2. * D) * G + (2. * D - 1.) / (2. * D) * Gs - 1./D * U * (N + (CC * Nd) / 2.);
#endif // BINARY_SCALAR
#endif

			}
//			}
		}
		//if(Fn < c.eps)
		//	printf("\nFn error!\n");
		if(outletflag == 0)
		{
			double Fbb = F - 2. * U;
			vlb.Ro.cor(ii0) += Fbb - Fn;

			vlb.F(dirp).n(ii0) = Fn;
#ifdef BINARY_FLUID
			vlb.G(dirp).n(ii0) = Gn;
#endif

			double dF = (Fn + F - 2. * c.lb.F0[dir]) * Ro;
			//double dF = (Fn + F) * Ro;
			//dF -= 2. * c.lb.F0[dir] * (1 + 4.5 * (CC * UU) * (CC * UU) - 1.5 * (UU * UU)) * Ro;
			clVector3Dd P = CC * dF;

			bcSetForce(P, bl, vCc, vC0, vC1, vC2);
		}
		else
		{
			//Fn = 2 * vlb.F(dirp)(ii0) - vlb.F(dirp)(iis1);
			Fn = vlb.F(dirp)(ii0);
			//Fn = vlb.F(dirp)(iis1);
			//vlb.Ro.cor(ii0) += F - Fn;
			vlb.F(dirp).n(ii0) = Fn;
			vlb.M.PBC(ii0) = outletflag;
		}
	}
	BOOL testperiodicindex(clVector3Di i0, clVector3Di i1)
	{
		if(c.lb.bcPeriodicX == FALSE && i0.x != i1.x)
			return FALSE;
		if(c.lb.bcPeriodicY == FALSE && i0.y != i1.y)
			return FALSE;
		if(c.lb.bcPeriodicZ == FALSE && i0.z != i1.z)
			return FALSE;
		return TRUE;
	}
	void bcSetForce(clVector3Dd F, int bl, clVector3Dd vCc, 
									clVector3Dd vC0, clVector3Dd vC1, clVector3Dd vC2)
	{
////////////////////////////dTlb = 1////////////////////////////
		if(F.r2() == 0.) return;
		clVector3Dd R0 = vC0 - vCc,  R1 = vC1 - vCc,  R2 = vC2 - vCc;
		clVector3Dd K0 = R0 ^ F,  K1 = R1 ^ F,  K2 = R2 ^ F;

		double den = K0.x * (R1.x - R2.x) + K0.y * (R1.y - R2.y) + K0.z * (R1.z - R2.z) + K1.x * R2.x + K1.y * R2.y + K1.z * R2.z;
		double c0, c1, c2;

		if(den != 0.)
		{
			double nom0 = K1.x * R2.x + K1.y * R2.y + K1.z * R2.z;
			double nom1 = K0.x * R2.x + K0.y * R2.y + K0.z * R2.z;
			double nom2 = K0.x * R1.x + K0.y * R1.y + K0.z * R1.z;
			c0 = nom0 / den;
			c1 = -nom1 / den;
			c2 = nom2 / den;
		}
		else
		{
			c0 = 1./3.;
			c1 = 1./3.;
			c2 = 1./3.;
			printf("\nnominator equals to zero! (bcSetForce)\n");
		}

		clVector3Dd F0 = F * c0,  F1 = F * c1,  F2 = F * c2;
		//clVector3Dd vMtot = R0 ^ F0 + R1 ^ F1 +  R2 ^ F2;
		
		//Modified
		//bcSetForce(bl, F, F0, F1, F2);
		bcSetForce(bl, vCc, F, F0, F1, F2);
		//Modified

	}

	//Modified
	//void bcSetForce(int bl, clVector3Dd F, clVector3Dd F0, clVector3Dd F1, clVector3Dd F2)
	void bcSetForce(int bl, clVector3Dd vCc, clVector3Dd F, clVector3Dd F0, clVector3Dd F1, clVector3Dd F2)
	{			
		clVector3Dd vCn = clVector3Dd(0.0,0.0,0.0), vRad = clVector3Dd(0.0,0.0,0.0);
		
		if( c.ls.O(c.ls.BL(bl,0)) == (NUM_VESSELS-1))
		{
			clVector3Dd vC0 = vls.C(c.ls.BL(bl,0)), vC1 = vls.C(c.ls.BL(bl,1)), vC2 = vls.C(c.ls.BL(bl,2));
			vCn = (vC1 - vC0) ^ (vC2 - vC0);
			double midy = (c.lb.nY-1.0)/2.0, midz = (c.lb.nZ-1.0)/2.0;
			clVector3Dd midDiff = vCc - clVector3Dd(vCc.x,midy,midz);
			vRad = clVector3Dd(0.0,midDiff.y/midDiff.abs(),midDiff.z/midDiff.abs());
		}

		//if(c.t.Time==100.00)
		//{
			//clVector3Dd pt0 = vls.C(c.ls.BL(bl,0));
			//clVector3Dd pt1 = vls.C(c.ls.BL(bl,1));
			//clVector3Dd pt2 = vls.C(c.ls.BL(bl,2));
			//printf("\n Face=%d, vRad= %f, %f, %f, with vCn= %f, %f, %f, with vCc= %f, %f, %f, pt1= %f, %f, %f, pt2=  %f, %f, %f, pt3= %f, %f, %f", \
			//	bl, vRad.x, vRad.y, vRad.z, vCn.x, vCn.y, vCn.z, vCc.x, vCc.y, vCc.z, pt0.x, pt0.y, pt0.z, pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z);
			//printf("\n Dot product is = %f", (vRad.x*vCn.x + vRad.y*vCn.y + vRad.z*vCn.z));
		//}

		c.ls.BLFb(bl,0) += F0;
		c.ls.BLFb(bl,1) += F1;
		c.ls.BLFb(bl,2) += F2;

		if( (c.ls.O(c.ls.BL(bl,0)) != (NUM_VESSELS-1)) || ((vRad.x*vCn.x + vRad.y*vCn.y + vRad.z*vCn.z) < 0 ) )
		{	
			//printf("\n FSI force added at face num = %d\n",bl); //Outputs which row in lsbl the F.b is added			
			vls.F.b(c.ls.BL(bl,0)) += F0;
			vls.F.b(c.ls.BL(bl,1)) += F1;
			vls.F.b(c.ls.BL(bl,2)) += F2;
		}
		//Modified
	}

	double bcLungUniOsc()
	{
		double vel;
#ifdef OSC_FLOW
		if(MODFAST(c.t.Time + c.t.StartTime,CYCLE_PERIOD) <= INHALE)
			vel = ::c.t.oscillator1();
		else
			vel = ::c.t.oscillator2();
#else
		vel = INLET_VEL_MAG;
#endif // OSC_FLOW
		return vel;
	}

	void bcSetBoundaryNode(clVector3Di ii0, clVector3Di ii1, int dir, double dist, clVector3Dd vVavr, 
									int bl, clVector3Dd vCc, clVector3Dd vC0, 
									clVector3Dd vC1, clVector3Dd vC2, clVector3Dd vCn)
	{
#ifdef FIND_JUNCTIONS // Find nodes where distr. functions cross both the valve and the vessel
		int mo = vlb.M.O(ii0);
		int ns0 = c.ls.BL(bl,1);
		int o = c.ls.O(ns0);

		if(mo == -1){
			vlb.M.O(ii0) = o;
			//vlb.M.V(ii0) = vls.V(ns0);
		}
		else if(mo == 0 && o > 0)
		{
			vlb.M.O(ii0) = -2;
			vlb.M.OS(ii0) = 1;
//#ifdef SOLID_JUNCTIONS
//			vlb.M(ii0) = LB_SOLID;
//#endif // SOLID_JUNCTIONS
			return;
		}
		else if(mo > 0 && o == 0)
		{
			vlb.M.O(ii0) = -2;
			vlb.M.OS(ii0) = 1;
			//vlb.M.V(ii0) = vls.V(ns0);
//#ifdef SOLID_JUNCTIONS
			//vlb.M(ii0) = LB_SOLID;
//#endif // SOLID_JUNCTIONS
			return;
		}
#endif //  FIND_JUNCTIONS

		ii0 = MOD(ii0, c.lb.nn);
		if(vlb.M.s(ii0) == TRUE)
				return;

		int BMflag = vlb.M.Bflag(ii0);
		if(BMflag & c.lb.LF[dir])
			return;
		vlb.M.Bflag(ii0) = BMflag | c.lb.LF[dir];

		//clVector3Dd NL = c.ls.directionBL(c.ls.BNind(c.ls.BL(bl,0))); //fix for capsules
		//double norm = vCn * NL;
		//if(norm > 0.)
		//	vlb.M(ii0) = c.ls.BF(c.ls.BNind(c.ls.BL(bl,0))) + 1;
		//else
		//	vlb.M(ii0) = c.ls.BF(c.ls.BNind(c.ls.BL(bl,0)));
		//if(vlb.M(ii0) == LB_SOLID || bcTestCrossSolid(ii0, vCc, vCn, vVavr) == TRUE)
		if(vlb.M(ii0) == LB_SOLID )
		{
			//clVector3Dd NL = c.ls.directionBL(c.ls.BNind(c.ls.BL(bl,0))) + c.ls.directionBL(c.ls.BNind(c.ls.BL(bl,1))) + c.ls.directionBL(c.ls.BNind(c.ls.BL(bl,2))); //fix for capsules
			//if(vCn * NL > 0)
			//	bcSetNewF(ii0, vVavr, c.ls.BF(c.ls.BNind(c.ls.BL(bl,0))) + 1, c.ls.O(c.ls.BL(bl,0)), vCn.n());
			//else
				bcSetNewF(ii0, vVavr, c.ls.BF(c.ls.BNind(c.ls.BL(bl,0))), c.ls.O(c.ls.BL(bl,0)), vCn.n());
		}

		bcSetNodeBC(ii0, dir, dist, vVavr, bl, vCc, vC0, vC1, vC2);
		
#ifdef LB_DEBUG
		vlb.B.D[dir](ii0) = dist;
		vlb.B.Ux[dir](ii0) = vVavr.x;
		vlb.B.Uy[dir](ii0) = vVavr.y;
		vlb.B.Uz[dir](ii0) = vVavr.z;
#endif //LB_DEBUG
	}
	BOOL bcTestCrossSolid(clVector3Di ii0, clVector3Dd vCc, clVector3Dd vCn, clVector3Dd vVavr)
	{
		//return FALSE;

		//if(vCc.y < 9.)
		//	printf("node %d %d %d %g %g %g\n", ii0.x, ii0.y, ii0.z, vCc.x, vCc.y, vCc.z);
		clVector3Dd N = vCn.n();
		double cutdist = vVavr * N * 2.;

		double dist = ((clVector3Dd)ii0 - vCc) * N;

		if(dist > fabs(cutdist) || dist == 0.)
			return FALSE;

		clVector3Dd Dnew = N * dist;
		clVector3Dd D = vlb.M.D(ii0);
		double Dabs = D.abs();

		if(Dabs == 0 || Dabs > dist)
			vlb.M.D(ii0) = Dnew;

		if(Dabs > 0.)
			return FALSE;

		clVector3Dd Dold = vlb.M.Dold(ii0);
		double Doldabs = Dold.abs();

		if(Doldabs == 0.)
			return FALSE;

		double norm = Dnew.n() * Dold.n();

		if(norm < -0.2)
		{
//			printf("set new distribution %d %d %d %g %g %g\n", ii0.x, ii0.y, ii0.z, vCc.x, vCc.y, vCc.z);
			return TRUE;
		}

		return FALSE;
	}

	BOOL bcFindIntersection(clVector3Di &vCcut0, clVector3Di &vCcut1, clVector3Dd &vCcut, double &dist, 
						clVector3Dd vL0, clVector3Dd vLd, 
						clVector3Dd vC0, clVector3Dd vC1, clVector3Dd vC2, clVector3Dd vCn, clVector3Di flPO)
	{
		double distF;
		if(bcFindIntersectionLinePlane(vCcut, distF, vL0, vLd, vC0, vC1, vC2) == FALSE)
			return FALSE;

		if(testInside(vCcut, vC0, vC1, vC2) == FALSE)
			return FALSE;

		int n_cut_0 = (int)ceil(distF);
		int n_cut = (int)ceil(distF);
		double dist_0 = 1. - (n_cut - (double)distF);
		dist = 1. - (n_cut - (double)distF);

		vCcut1 = (clVector3Di)(vL0 + vLd * n_cut);
		vCcut0 = vCcut1 - (clVector3Di)vLd;

		clVector3Di ii0 = MOD(vCcut0, c.lb.nn);

		if(flPO.x == FALSE && ii0.x != vCcut0.x)
		{
//			printf("node: %d %d %d\n", ii0.x, ii0.y, ii0.z);
			return FALSE; 
		}

		if(flPO.y == FALSE && ii0.y != vCcut0.y)
			return FALSE;

		if(flPO.z == FALSE && ii0.z != vCcut0.z)
			return FALSE;


		if(vlb.M.s(ii0) == TRUE)
		{
			if(dist < c.eps)
			{
				dist = 1.;
				n_cut--;
				vCcut1 = (clVector3Di)(vL0 + vLd * n_cut);
				vCcut0 = vCcut1 - (clVector3Di)vLd;

				clVector3Di ii0 = MOD(vCcut0, c.lb.nn);

				if(flPO.x == FALSE && ii0.x != vCcut0.x)
					return FALSE;

				if(flPO.y == FALSE && ii0.y != vCcut0.y)
					return FALSE;

				if(flPO.z == FALSE && ii0.z != vCcut0.z)
					return FALSE;

				if(vlb.M.s(ii0) == TRUE)
					return FALSE;
			}
			else
				return FALSE;
		}

		return TRUE;
	}

	BOOL bcFindIntersectionLinePlane(clVector3Dd &vC, double &dist, clVector3Dd vL0, clVector3Dd vLd, 
									clVector3Dd vP0, clVector3Dd vP1, clVector3Dd vP2)
	{
		clVector3Dd vP10 = vP1 - vP0, vP20 = vP2 - vP0;
		clVector3Dd vN = vP10 ^ vP20;
		double den = vN * vLd;
		if(den == 0.)
			return FALSE;

		dist = vN * (vP0 - vL0) / den;
		vC = vL0 + vLd * dist;
		return TRUE;
	}

	BOOL testInside(clVector3Dd vd, clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
	{
		clVector3Dd v10 = v1 - v0, v20 = v2 - v0, v21 = v2 - v1, vd0 = vd - v0, vd1 = vd - v1;
		clVector3Dd stot = v10 ^ v20, s0 = vd0 ^ v10, s1 = vd0 ^ v20, s2 = vd1 ^ v21;

		if(fabs(stot.r() - s0.r() - s1.r() - s2.r()) < c.eps)
			return TRUE;
		return FALSE;
	}

	clVector3Di min3(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
	{
		double x = v0.x, y = v0.y, z = v0.z;
		if(v1.x < x) x = v1.x;
		if(v2.x < x) x = v2.x;
		if(v1.y < y) y = v1.y;
		if(v2.y < y) y = v2.y;
		if(v1.z < z) z = v1.z;
		if(v2.z < z) z = v2.z;
		return clVector3Di((int)floor(x), (int)floor(y), (int)floor(z));
	}

	clVector3Di max3(clVector3Dd v0, clVector3Dd v1, clVector3Dd v2)
	{
		double x = v0.x, y = v0.y, z = v0.z;
		if(v1.x > x) x = v1.x;
		if(v2.x > x) x = v2.x;
		if(v1.y > y) y = v1.y;
		if(v2.y > y) y = v2.y;
		if(v1.z > z) z = v1.z;
		if(v2.z > z) z = v2.z;
		return clVector3Di((int)ceil(x), (int)ceil(y), (int)ceil(z));
	}

	// One lattice link cut by a surface triangle, found by the (read-only) search.
	struct BCHit
	{
		clVector3Di c0, c1;
		int dir;
		double dist;
		clVector3Dd vavr, vcc, v0, v1, v2, vn;
	};
	std::vector< std::vector<BCHit> > bcHits;

	void bcFindBoundaryNodes()
	{
#ifdef LB_DEBUG
		vlb.B.zeros();
#endif //LB_DEBUG

#ifdef BC_USE_OPENMP
		// Two phases. The geometric search (~98% of a time step) only reads the
		// structure and the static mask, so it runs in parallel and records each
		// triangle's cut links. The links are then applied serially in the original
		// triangle order: bcSetBoundaryNode claims a link for the FIRST triangle that
		// reaches it (Bflag) and accumulates into shared arrays, so applying in order
		// reproduces the serial result exactly.
		if((int)bcHits.size() != c.ls.nBL)
			bcHits.resize(c.ls.nBL);
		int i;
		#pragma omp parallel for default(shared) private(i) schedule(dynamic, 64)
		for(i = 0; i < c.ls.nBL; i++)
		{
			bcHits[i].clear();
			bcFindNodes(i, &bcHits[i]);
		}
		for(i = 0; i < c.ls.nBL; i++)
			for(size_t k = 0; k < bcHits[i].size(); k++)
			{
				BCHit &h = bcHits[i][k];
				bcSetBoundaryNode(h.c0, h.c1, h.dir, h.dist, h.vavr, i, h.vcc, h.v0, h.v1, h.v2, h.vn);
			}
#else //BC_USE_OPENMP
		for(int i = 0; i < c.ls.nBL; i++)
			bcFindNodes(i);
#endif //BC_USE_OPENMP
	}

	void bcEmitHit(std::vector<BCHit> *out, clVector3Di vCcut0, clVector3Di vCcut1, int dir, double dist,
					clVector3Dd vVavr, int bl, clVector3Dd vCcut, clVector3Dd vC0, clVector3Dd vC1,
					clVector3Dd vC2, clVector3Dd vCn)
	{
		if(out == NULL)
		{
			bcSetBoundaryNode(vCcut0, vCcut1, dir, dist, vVavr, bl, vCcut, vC0, vC1, vC2, vCn);
			return;
		}
		BCHit h;
		h.c0 = vCcut0; h.c1 = vCcut1; h.dir = dir; h.dist = dist; h.vavr = vVavr;
		h.vcc = vCcut; h.v0 = vC0; h.v1 = vC1; h.v2 = vC2; h.vn = vCn;
		out->push_back(h);
	}

	clVector3Dd bcVavr(clVector3Dd vCc, clVector3Dd vC0, clVector3Dd vC1, clVector3Dd vC2, 
						clVector3Dd vV0, clVector3Dd vV1, clVector3Dd vV2)
	{
		clVector3Dd vN = (vC1 - vC0) ^ (vC2 - vC0);
		clVector3Dd vP12;
		double distCP;
		if(bcFindIntersectionLinePlane(vP12, distCP, vC0, vCc - vC0, vC1, vC1+vN, vC2) == FALSE)
		{
			if(vCc == vC0)
				return vV0;
			if(vCc == vC1)
				return vV1;
			if(vCc == vC2)
				return vV2;

			return (vV0 + vV1 + vV2) / 3.;
		}

		clVector3Dd vP1 = vP12 - vC1, vP2 = vP12 - vC2;
		double dP1 = vP1.r(), dP2 = vP2.r();
		clVector3Dd vVp = (vV1 * dP2 +  vV2 * dP1) / (dP1 + dP2);

		clVector3Dd vc0 = vCc - vC0, vcP = vCc - vP12;
		double dc0 = vc0.r(), dcP = vcP.r();
		clVector3Dd vVavr = (vV0 * dcP +  vVp * dc0) / (dc0 + dcP);
		return vVavr;
	}

	void bcFindNodes(int bl, std::vector<BCHit> *out = NULL)
	{
		int n0 = c.ls.BL(bl,0), n1 = c.ls.BL(bl,1),  n2 = c.ls.BL(bl,2);
		if(n0 < 0)
			return;
		clVector3Dd vC0 = vls.C.s(n0, n0), vC1 = vls.C.s(n1, n0), vC2 = vls.C.s(n2, n0);
		clVector3Dd vV0 = vls.V.avr(n0), vV1 = vls.V.avr(n1), vV2 = vls.V.avr(n2);
		clVector3Dd vCn = (vC1 - vC0) ^ (vC2 - vC0);
		clVector3Di vCmin = min3(vC0, vC1, vC2), vCmax = max3(vC0, vC1, vC2);
		clVector3Di flPO = c.ls.PO(c.ls.O(n0));
		if(flPO.x == FALSE)
		{
			if(vCmin.x < 0) vCmin.x = 0;
			if(vCmin.x > c.nX-1) return;
			if(vCmax.x < 0) return;
			if(vCmax.x > c.nX-1) vCmax.x = c.nX-1;
			//if(vC0.x < vCmin.x) vC0.x = vCmin.x;
			//if(vC1.x < vCmin.x) vC1.x = vCmin.x;
			//if(vC2.x < vCmin.x) vC2.x = vCmin.x;
		}
		if(flPO.y == FALSE)
		{
			if(vCmin.y < 0) vCmin.y = 0;
			if(vCmin.y > c.nY-1) return;
			if(vCmax.y < 0) return;
			if(vCmax.y > c.nY-1) vCmax.y = c.nY-1;
		}
		if(flPO.z == FALSE)
		{
			if(vCmin.z < 0) vCmin.z = 0;
			if(vCmin.z > c.nZ-1) return;
			if(vCmax.z < 0) return;
			if(vCmax.z > c.nZ-1) vCmax.z = c.nZ-1;
		}

		for(int dir = 1; dir < c.lb.nL; dir+=2)
			bcFindDirection(dir, bl, vC0, vC1, vC2, vCn, vCmin, vCmax, vV0, vV1, vV2, flPO, out);
	}

	void bcFindDirection(int dir, int bl, clVector3Dd vC0, clVector3Dd vC1, clVector3Dd vC2, 
								clVector3Dd vCn, clVector3Di vCmin, clVector3Di vCmax, 
								clVector3Dd vV0, clVector3Dd vV1, clVector3Dd vV2, clVector3Di flPO,
								std::vector<BCHit> *out = NULL)
	{
		double dist;
		clVector3Dd Cdir = c.lb.C[dir];
		double Cndir = vCn * Cdir;
		if(Cndir == 0.)
			return;
		if(Cndir > 0.)
		{
			dir = c.lb.rev[dir];
			Cdir = c.lb.C[dir];
		}
		clVector3Di vCcut0, vCcut1;
		clVector3Dd vCcut;
		clVector3Dd vVavr;

		if(c.lb.Cx[dir] != 0)
		{
			int i;
			if(c.lb.Cx[dir] == 1) i = vCmin.x;
			if(c.lb.Cx[dir] == -1) i = vCmax.x;

			for(int j = vCmin.y; j <= vCmax.y; j++)
				for(int z = vCmin.z; z <= vCmax.z; z++)
				{
					if(bcFindIntersection(vCcut0, vCcut1, vCcut, dist, clVector3Dd(i, j, z), 
												Cdir, vC0, vC1, vC2, vCn, flPO) == TRUE)
					{
						vVavr = bcVavr(vCcut, vC0, vC1, vC2, vV0, vV1, vV2);
						bcEmitHit(out, vCcut0, vCcut1, dir, dist, vVavr, bl, vCcut, vC0, vC1, vC2, vCn);
					}
				}
		}
		if(c.lb.Cy[dir] != 0)
		{
			int j;
			if(c.lb.Cy[dir] == 1) j = vCmin.y;
			if(c.lb.Cy[dir] == -1) j = vCmax.y;

			for(int i = vCmin.x; i <= vCmax.x; i++)
				for(int z = vCmin.z; z <= vCmax.z; z++)
				{
					if(bcFindIntersection(vCcut0, vCcut1, vCcut, dist, clVector3Dd(i, j, z), 
												Cdir, vC0, vC1, vC2, vCn, flPO) == TRUE)
					{
						vVavr = bcVavr(vCcut, vC0, vC1, vC2, vV0, vV1, vV2);
						bcEmitHit(out, vCcut0, vCcut1, dir, dist, vVavr, bl, vCcut, vC0, vC1, vC2, vCn);
					}
				}
		}
		if(c.lb.Cz[dir] != 0)
		{
			int z;
			if(c.lb.Cz[dir] == 1) z = vCmin.z;
			if(c.lb.Cz[dir] == -1) z = vCmax.z;

			for(int i = vCmin.x; i <= vCmax.x; i++)
				for(int j = vCmin.y; j <= vCmax.y; j++)
				{
					if(bcFindIntersection(vCcut0, vCcut1, vCcut, dist, clVector3Dd(i, j, z), 
												Cdir, vC0, vC1, vC2, vCn, flPO) == TRUE)
					{
						vVavr = bcVavr(vCcut, vC0, vC1, vC2, vV0, vV1, vV2);
						bcEmitHit(out, vCcut0, vCcut1, dir, dist, vVavr, bl, vCcut, vC0, vC1, vC2, vCn);
					}
				}
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////boundary conditions///////////////////////////////////////////////////////////////////////////////

////////////////////////////////solid wall BC/////////////////////////////////////////////////////////////////////////
	

	void bcWallLB(int i, int j, int z, int dir, clVector3Dd UU2Cs2)
	{
		if(vlb.M.Bflag(i,j,z) & c.lb.LF[dir])
			return;

		int dirp = c.lb.rev[dir];
		double CUU2Cs2 = c.lb.C[dir] * UU2Cs2;
		vlb.F(dirp).n(i,j,z) = vlb.F(dir)(i,j,z) - c.lb.F0[dir] * CUU2Cs2;
#ifdef BINARY_FLUID
		double N = vlb.N(i,j,z);
		double CdN = c.lb.C[dir] * vlb.N.d(i,j,z);
		vlb.G(dirp).n(i,j,z) = vlb.G(dir)(i,j,z) - (N + CdN * 0.5) * c.lb.F0[dir] * CUU2Cs2; //check how it works
#endif
		vlb.M.Bflag(i,j,z) = vlb.M.Bflag(i,j,z) | c.lb.LF[dir];
	}
	void bcWallLBx0()
	{
		clVector3Dd U0 = c.lb.UX0 * 2. / c.lb.Cs2;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
//				if(vlb.M(0,j,z) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == -1)
							bcWallLB(0, j, z, k, U0);
			}
	}
	void bcWallLBx1()
	{
		int nXm = c.lb.nX-1;
		clVector3Dd U1 = c.lb.UX1 * 2. / c.lb.Cs2;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
//				if(vlb.M(nXm,j,z) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cx[k] == 1)
							bcWallLB(nXm, j, z, k, U1);
			}
	}

	void bcWallLBy0()
	{
		clVector3Dd U0 = c.lb.UY0 * 2. / c.lb.Cs2;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
//				if(vlb.M(i,0,z) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == -1)
							bcWallLB(i, 0, z, k, U0);
			}
	}
	void bcWallLBy1()
	{
		int nYm = c.lb.nY-1;
		clVector3Dd U1 = c.lb.UY1 * 2. / c.lb.Cs2;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
//				if(vlb.M(i,nYm,z) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cy[k] == 1)
							bcWallLB(i, nYm, z, k, U1);
			}
	}

	void bcWallLBz0()
	{
		clVector3Dd U0 = c.lb.UZ0 * 2. / c.lb.Cs2;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
//				if(vlb.M(i,j,0) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == -1)
							bcWallLB(i, j, 0, k, U0);
			}
	}
	void bcWallLBz1()
	{
		int nZm = c.lb.nZ-1;

		clVector3Dd U1 = c.lb.UZ1 * 2. / c.lb.Cs2;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
//				if(vlb.M(i,j,nZm) > LB_SOLID)
					for(int k = 1; k < c.lb.nL; k++)
						if(c.lb.Cz[k] == 1)
							bcWallLB(i, j, nZm, k, U1);
			}
	}

////////////////////////////////////////// Pressure BC ////////////////////////////////////
	clVector3Dd bcPressVel(int i, int j, int z, int norm, int n, double Ro)
	{
		clVector3Dd U0;
		double sum;
		int nm = (int)(-1 * n);
		for(int ii = 0; ii < 3; ii++)
		{
			sum = 0;
			if(ii == norm)
			{
				for(int ic = 0; ic < 19; ic++)
				{
					if(c.lb.C[ic].e(ii) == nm)
						sum += 2 * vlb.F(ic)(i,j,z);
					else if(c.lb.C[ic].e(ii) == 0)
						sum += vlb.F(ic)(i,j,z);
				}
				U0.e(ii) = n * (1 - sum / Ro);
			}
			else
			{
				for(int ic = 1; ic < 19; ic++)
				{
					if(c.lb.C[ic].e(norm) == 0)
						sum += c.lb.C[ic].e(ii) * vlb.F(ic)(i,j,z);
				}
				U0.e(ii) = 1.5 * sum / Ro;
			}
		}
		return U0;
	}
	void bcPressLBx0()
	{
		//clVector3Dd U0 = c.lb.UX0 * 2. / c.lb.Cs2;
		double Ro = 1. + c.lb.RoX0 * vlb.J.invel;		
		
		int n = 1;
		int norm = 0;
		clVector3Dd U0;
		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				if(vlb.M(0,j,z) == LB_SOLID)
					continue;

				U0 = bcPressVel(0, j, z, norm, n, Ro);
				U0 = clVector3Dd(U0.x,0.0,0.0);
				//U0 = U0 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cx[k] == -1)
					{
						//bcWallLB(0, j, z, k, U0);
						if(vlb.M.Bflag(0,j,z) & c.lb.LF[k])
							return;
						
						int dirp = c.lb.rev[k];
						clVector3Dd surfn = clVector3Dd(-1.0,0.0,0.0);
						clVector3Dd tang = c.lb.C[k] - clVector3Dd((c.lb.C[k]*surfn) * surfn.x, (c.lb.C[k]*surfn) * surfn.y, (c.lb.C[k]*surfn) * surfn.z);					

						double refsum = 0.0;
						for(int m = 0; m < c.lb.nL; m++)
							refsum = 0.5 * vlb.F(m)(0,j,z) * (tang*c.lb.C[k]) * (1.0-fabs(c.lb.C[k]*surfn));
						
						vlb.F(dirp).n(0,j,z) = vlb.F(k)(0,j,z) - Ro/6.0*(c.lb.C[k]*U0) - Ro/3.0*(tang*U0) + refsum;

						vlb.M.Bflag(0,j,z) = vlb.M.Bflag(0,j,z) | c.lb.LF[k];
					}
			}
	}
	void bcPressLBx1()
	{
		double Ro = 1. + c.lb.RoX1 * vlb.J.invel;

		int nXm = c.lb.nX-1;		
		int n = -1;
		int norm = 0;
		clVector3Dd U1;
		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				if(vlb.M(nXm,j,z) == LB_SOLID)
					continue;

				U1 = bcPressVel(nXm, j, z, norm, n, Ro);
				U1 = clVector3Dd(U1.x,0.0,0.0);
				//U1 = U1 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cx[k] == 1)
					{
						//bcWallLB(nXm, j, z, k, U1);
						if(vlb.M.Bflag(nXm,j,z) & c.lb.LF[k])
							return;
						
						int dirp = c.lb.rev[k];
						clVector3Dd surfn = clVector3Dd(1.0,0.0,0.0);
						clVector3Dd tang = c.lb.C[k] - clVector3Dd((c.lb.C[k]*surfn) * surfn.x, (c.lb.C[k]*surfn) * surfn.y, (c.lb.C[k]*surfn) * surfn.z);					

						double refsum = 0.0;
						for(int m = 0; m < c.lb.nL; m++)
							refsum = 0.5 * vlb.F(m)(nXm,j,z) * (tang*c.lb.C[k]) * (1.0-fabs(c.lb.C[k]*surfn));
						
						vlb.F(dirp).n(nXm,j,z) = vlb.F(k)(nXm,j,z) - Ro/6.0*(c.lb.C[k]*U1) - Ro/3.0*(tang*U1) + refsum;

						vlb.M.Bflag(nXm,j,z) = vlb.M.Bflag(nXm,j,z) | c.lb.LF[k];
					}
			}
	}
	void bcPressLBy0()
	{
		double Ro = 1. + c.lb.RoY0 * vlb.J.invel;
		int n = 1;
		int norm = 1;
		clVector3Dd U0;
		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				if(vlb.M(i,0,z) == LB_SOLID)
					continue;

				int n = 1;
				int norm = 1;
				U0 = bcPressVel(i, 0, z, norm, n, Ro);
				U0 = U0 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cy[k] == -1)
						bcWallLB(i, 0, z, k, U0);
			}
	}
	void bcPressLBy1()
	{
		int nYm = c.lb.nY-1;
		double Ro = 1. + c.lb.RoY1 * vlb.J.invel;
		int n = -1;
		int norm = 1;
		clVector3Dd U1;
		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				if(vlb.M(i,nYm,z) == LB_SOLID)
					continue;

				U1 = bcPressVel(i, nYm, z, norm, n, Ro);
				U1 = U1 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cy[k] == 1)
						bcWallLB(i, nYm, z, k, U1);
			}
	}
	void bcPressLBz0()
	{
		double Ro = 1. + c.lb.RoZ0 * vlb.J.invel;
		int n = 1;
		int norm = 2;
		clVector3Dd U0;
		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
				if(vlb.M(i,j,0) == LB_SOLID)
					continue;

				U0 = bcPressVel(i, j, 0, norm, n, Ro);
				U0 = U0 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cz[k] == -1)
						bcWallLB(i, j, 0, k, U0);
			}
	}
	void bcPressLBz1()
	{
		int nZm = c.lb.nZ-1;
		double Ro = 1. + c.lb.RoZ1 * vlb.J.invel;
		int n = -1;
		int norm = 2;
		clVector3Dd U1;
		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
				if(vlb.M(i,j,nZm) == LB_SOLID)
					continue;

				U1 = bcPressVel(i, j, nZm, norm, n, Ro);
				U1 = U1 * Ro * 2. / c.lb.Cs2;
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cz[k] == 1)
						bcWallLB(i, j, nZm, k, U1);
			}
	}

////////////////////////////////symmetry BC/////////////////////////////////////////////////////////////////////////
	void bcSymmetryLB(int i, int j, int z, int ip, int jp, int zp, int dir, int ind)
	{
		if(vlb.M.Bflag(i,j,z) & c.lb.LF[dir])
			return;

		int dirp = c.lb.rev[dir];
		if(dir > 2 * NUMBER_DIMENSIONS)
		{
			int dirn = c.lb.per[dir];
			if(c.lb.C[dir].e(ind) == -c.lb.C[dirn].e(ind))
				dirn = c.lb.rev[dirn];
			vlb.F(dirp).n(i,j,z) = vlb.F(dirn)(ip,jp,zp);
#ifdef BINARY_FLUID
			vlb.G(dirp).n(i,j,z) = vlb.G(dirn)(ip,jp,zp);
#endif
		}
		else
		{
			vlb.F(dirp).n(i,j,z) = vlb.F(dir)(i,j,z);
#ifdef BINARY_FLUID
			vlb.G(dirp).n(i,j,z) = vlb.G(dir)(i,j,z);
#endif
		}
		vlb.M.Bflag(i,j,z) = vlb.M.Bflag(i,j,z) | c.lb.LF[dir];
	}
	void bcSymmetryLBx0()
	{
		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cx[k] == -1)
						bcSymmetryLB(0, j, z, 0, j+c.lb.Cy[k], z+c.lb.Cz[k], k, 0);
	}
	void bcSymmetryLBx1()
	{
		int nXm = c.lb.nX-1;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cx[k] == 1)
						bcSymmetryLB(nXm, j, z, nXm, j+c.lb.Cy[k], z+c.lb.Cz[k], k, 0);
	}

	void bcSymmetryLBy0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cy[k] == -1)
						bcSymmetryLB(i, 0, z, i+c.lb.Cx[k], 0, z+c.lb.Cz[k], k, 1);
	}
	void bcSymmetryLBy1()
	{
		int nYm = c.lb.nY-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cy[k] == 1)
						bcSymmetryLB(i, nYm, z, i+c.lb.Cx[k], nYm, z+c.lb.Cz[k], k, 1);
	}

	void bcSymmetryLBz0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cz[k] == -1)
						bcSymmetryLB(i, j, 0, i+c.lb.Cx[k], j+c.lb.Cy[k], 0, k, 2);
	}
	void bcSymmetryLBz1()
	{
		int nZm = c.lb.nZ-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
				for(int k = 1; k < c.lb.nL; k++)
					if(c.lb.Cz[k] == 1)
						bcSymmetryLB(i, j, nZm, i+c.lb.Cx[k], j+c.lb.Cy[k], nZm, k, 2);
	}
/////////////////////////////////////coupling with coarse grid//////////////////////////////////////////
	void bcCoarseLB(int i, int j, int z, int dir)
	{
		vlb.F(dir).n(i,j,z) = vlb.C.Fc2f(dir, clVector3Di(i,j,z));
	}
	void bcCoarseLBx0()
	{
		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cx[k] == 1)
						bcCoarseLB(0, j, z, k);
	}
	void bcCoarseLBx1()
	{
		int nXm = c.lb.nX-1;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cx[k] == -1)
						bcCoarseLB(nXm, j, z, k);
	}

	void bcCoarseLBy0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cy[k] == 1)
						bcCoarseLB(i, 0, z, k);
	}
	void bcCoarseLBy1()
	{
		int nYm = c.lb.nY-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cy[k] == -1)
						bcCoarseLB(i, nYm, z, k);
	}

	void bcCoarseLBz0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cz[k] == 1)
						bcCoarseLB(i, j, 0, k);
	}
	void bcCoarseLBz1()
	{
		int nZm = c.lb.nZ-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
				for(int k = 0; k < c.lb.nL; k++)
//					if(c.lb.Cz[k] == -1)
						bcCoarseLB(i, j, nZm, k);
	}
/////////////////////////////////////free flow BC//////////////////////////////////////////	void bcFreeFlowLBx0()
	void bcFreeFlowLBx0()
	{
		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cx[k] == 1)
					{
						vlb.F.F[k].npx(0,j,z) = vlb.F.F[k](1,j,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npx(0,j,z) = vlb.G.F[k](1,j,z);
#endif
						vlb.M.Bflag(0,j,z) = vlb.M.Bflag(0,j,z) | c.lb.LF[k];
					}
				}
			}
	}
	void bcFreeFlowLBx1()
	{
		int nXm = c.lb.nX-1;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cx[k] == -1)
					{
						vlb.F.F[k].npx(nXm,j,z) = vlb.F.F[k](nXm-1,j,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npx(nXm,j,z) = vlb.G.F[k](nXm-1,j,z);
#endif
						vlb.M.Bflag(nXm,j,z) = vlb.M.Bflag(nXm,j,z) | c.lb.LF[k];
					}
				}
			}
	}

	void bcFreeFlowLBy0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cy[k] == 1)
					{
						vlb.F.F[k].npy(i,0,z) = vlb.F.F[k](i,1,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npy(i,0,z) = vlb.G.F[k](i,1,z);
#endif
						vlb.M.Bflag(i,0,z) = vlb.M.Bflag(i,0,z) | c.lb.LF[k];
					}
				}
			}
	}
	void bcFreeFlowLBy1()
	{
		int nYm = c.lb.nY-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cy[k] == -1)
					{
						vlb.F.F[k].npy(i,nYm,z) = vlb.F.F[k](i,nYm-1,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npy(i,nYm,z) = vlb.G.F[k](i,nYm-1,z);
#endif
						vlb.M.Bflag(i,nYm,z) = vlb.M.Bflag(i,nYm,z) | c.lb.LF[k];
					}
				}
			}
	}

	void bcFreeFlowLBz0()
	{
		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cz[k] == 1)
					{
						vlb.F.F[k].npz(i,j,0) = vlb.F.F[k](i,j,1);
#ifdef BINARY_FLUID
						vlb.G.F[k].npz(i,j,0) = vlb.G.F[k](i,j,1);
#endif
						vlb.M.Bflag(i,j,0) = vlb.M.Bflag(i,j,0) | c.lb.LF[k];
					}
				}
			}
	}

	void bcFreeFlowLBz1()
	{
		int nZm = c.lb.nZ-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cz[k] == -1)
					{
						vlb.F.F[k].npz(i,j,nZm) = vlb.F.F[k](i,j,nZm-1);
#ifdef BINARY_FLUID
						vlb.G.F[k].npz(i,j,nZm) = vlb.G.F[k](i,j,nZm-1);
#endif
						vlb.M.Bflag(i,j,nZm) = vlb.M.Bflag(i,j,nZm) | c.lb.LF[k];
					}
				}
			}
	}
////////////////////////////////////////////periodic BC//////////////////////////////////////////
	void bcPeriodicLBx()
	{
		int nXm = c.lb.nX-1;

		for(int j = 0; j < c.lb.nY; j++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cx[k] == -1) {
						vlb.F.F[k].npx(nXm,j,z) = vlb.F.F[k](0,j,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npx(nXm,j,z) = vlb.G.F[k](0,j,z);
#endif
					}
					if(c.lb.Cx[k] == 1) {
						vlb.F.F[k].npx(0,j,z) = vlb.F.F[k](nXm,j,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npx(0,j,z) = vlb.G.F[k](nXm,j,z);
#endif
					}
				}
			}
	}

	void bcPeriodicLBy()
	{
		int nYm = c.lb.nY-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int z = 0; z < c.lb.nZ; z++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cy[k] == -1) {
						vlb.F.F[k].npy(i,nYm,z) = vlb.F.F[k](i,0,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npy(i,nYm,z) = vlb.G.F[k](i,0,z);
#endif
					}
					if(c.lb.Cy[k] == 1) {
						vlb.F.F[k].npy(i,0,z) = vlb.F.F[k](i,nYm,z);
#ifdef BINARY_FLUID
						vlb.G.F[k].npy(i,0,z) = vlb.G.F[k](i,nYm,z);
#endif
					}
				}
			}
	}

	void bcPeriodicLBz()
	{
		int nZm = c.lb.nZ-1;

		for(int i = 0; i < c.lb.nX; i++)
			for(int j = 0; j < c.lb.nY; j++)
			{
				for(int k = 1; k < c.lb.nL; k++)
				{
					if(c.lb.Cz[k] == -1) {
						vlb.F.F[k].npz(i,j,nZm) = vlb.F.F[k](i,j,0);
#ifdef BINARY_FLUID
						vlb.G.F[k].npz(i,j,nZm) = vlb.G.F[k](i,j,0);
#endif
					}
					if(c.lb.Cz[k] == 1) {
						vlb.F.F[k].npz(i,j,0) = vlb.F.F[k](i,j,nZm);
#ifdef BINARY_FLUID
						vlb.G.F[k].npz(i,j,0) = vlb.G.F[k](i,j,nZm);
#endif
					}
				}
			}
	}
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////LB propagation step/////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
	void propagationLB()
	{
#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
#endif //LB_DEBUG
		vlb.F.propagation();
#ifdef BINARY_FLUID
		vlb.G.propagation();
#endif

#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
#endif //LB_DEBUG
	}
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////LB collision step///////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
	void collisionLB()
	{
		vlb.Ro.update();
		vlb.J.update();
#ifdef BINARY_FLUID
		vlb.N.update();
#endif

#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
		vlb.Ro.save2file();
		vlb.J.save2file();
#endif //LB_DEBUG
		

#ifdef LB_USE_MRTLB
		collisionMRTLB();
#endif //LB_USE_MRTLB
#ifdef LB_USE_LBGK
		BOOL flStress = FALSE;
		//BOOL flStress = TRUE;
		for(int i = 1; i < c.lb.nF; i++)
		{
			if(c.lb.lamda[i] != -1. || c.lb.lamdab[i] != -1.)
			{
				flStress = TRUE;
				break;
			}
		}
#ifdef BINARY_SCALAR
		flStress = FALSE;
#endif // BINARY_SCALAR
		if(flStress == TRUE)
			collision_with_stress_LBGK();
		else
			collision_without_stress_LBGK();
#endif //LB_USE_LBGK

		// builk free energy density
#ifdef LB_DEBUG
		vlb.F.save2file();
		vlb.G.save2file();
		vlb.Ro.save2file();
		vlb.J.save2file();
		vlb.PP.save2file();
#endif //LB_DEBUG
	}

	void collision_with_stress_LBGK()
	{
		vlb.PP.update();
		//return;
		int i, inout;
#ifdef LB_USE_OPENMP
#pragma omp parallel default(shared) private(i)
#endif //LB_USE_OPENMP
		{
			double lamda, lamdab;
			double PRoCorr, Jmag;
			clVector3Dd Fbody; 
			clTensorS3Dd Tunit;
			Tunit.unit();

#ifdef LB_USE_OPENMP
#pragma omp for
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
			{
				for(int j = 0; j < c.lb.nY; j++)
				{
					for(int z = 0; z < c.lb.nZ; z++)
					{
						char M = vlb.M(i,j,z);
						if(M != LB_SOLID)
						{
							lamda = c.lb.lamda[M];
							lamdab = c.lb.lamdab[M];
							Fbody = c.lb.F[M] * c.lb.FM(i,j,z);
							double Ro = vlb.Ro(i,j,z);
							inout = vlb.M.PBC(i,j,z);
							if(inout == 0)
							{
								PRoCorr = 1.;
							}
							else if(inout == 1)
							{
								PRoCorr = (1. + DELTA_RO * vlb.J.invel) / Ro;
							}
							else
							{
#ifndef OUTLET_INTERP
								PRoCorr = (1. - DELTA_RO * vlb.J.invel) / Ro;
#endif // OUTLET_INTERP
							}
							clVector3Dd  J = vlb.J(i,j,z);
#ifdef INOUTNORMAL
							if(inout != 0)
							{
								Jmag = J.r();
								J = c.ls.InoutNorms(inout) * Jmag;
								if(vlb.J.invel < 0.){
									J = J * (-1.);
								}
							}
#endif // INOUTNORMAL
							clTensorS3Dd PP = vlb.PP(i,j,z);
							clTensorS3Dd RoUU = clTensorS3Dd(J) / Ro;							

							clTensorS3Dd PPneq = PP	- (RoUU + Tunit * 1./3. * Ro);

							//clTensorS3Dd PPneqtr = Tunit * 1./3. * PPneq.tr_sum();
							//clTensorS3Dd PPneq_s = (PPneq - PPneqtr) * (1. + lamda) + PPneqtr * (1. + lamdab);

							clTensorS3Dd FFbody2 = clTensorS3Dd(J/Ro, Fbody);

							clVector3Dd JF = J + Fbody;
#ifdef BINARY_FLUID
							
							double N = vlb.N(i,j,z);
							clVector3Dd NU = J / Ro * N;


							double Fi = Ro * c.lb.Cs2 * log(Ro) - (N*N)*LBB_A/2. + (N*N*N*N)*LBB_B/4.;
							double dFidN = (-LBB_A + LBB_B*N*N) * N;
							double dFidRo = (log(Ro) + 1.) / 3.;
							double P0 = Ro * dFidRo + N * dFidN - Fi;

							double d2N = vlb.N.d2(i,j,z).sum();
							double dN2 = vlb.N.d(i,j,z).r2();
							double kd2N = LBB_KAPPA*d2N;
														
							double A = ( P0 + N*kd2N - 0.5*LBB_KAPPA*dN2) / c.lb.Cs2;

							clTensorS3Dd dN2t(vlb.N.d(i,j,z));

//							clTensorS3Dd RoUUPF = RoUU + PPneq_s + FFbody2 + d2N;
							clTensorS3Dd RoUUPF = RoUU + FFbody2 + dN2t*LBB_KAPPA;

							clTensorS3Dd NUU = RoUU / Ro * N;
							double GMu = LBB_GAMMA_NUMBER * (dFidN - kd2N) / c.lb.Cs2;
							
							double Ftot = 0., Gtot = 0., gTrackTot = 0.;

							for(int k = 1; k < c.lb.nL; k++)
							{
								double F, G;
								const clVector3Dd C = c.lb.C[k];
								clTensorS3Dd CC = clTensorS3Dd(C) - Tunit * 1./3.;

								double t =  RoUUPF * CC * 4.5;
								double Fnew = (A + 3. * (JF * C) + t) / c.lb.Ro[M];
								F = c.lb.F0[k] * Fnew;
								vlb.F(k)(i,j,z) = vlb.F(k)(i,j,z) * (1. - LB_TAU_NUMBER) + F * LB_TAU_NUMBER;
								Ftot += F;

								double Gnew = (GMu + 3. * (NU * C) + NUU * CC * 4.5);
								G = c.lb.F0[k] * Gnew;
								vlb.G(k)(i,j,z) = (vlb.G(k)(i,j,z) * (1. - LB_TAUN_NUMBER) + G * LB_TAUN_NUMBER);
								gTrackTot += vlb.G(k)(i,j,z);
								Gtot += G;

							}
							vlb.F(0)(i,j,z) = Ro / c.lb.Ro[M] - Ftot;
							vlb.G(0)(i,j,z) = N - Gtot;
#endif
#ifndef BINARY_FLUID
							clTensorS3Dd PPneqtr = Tunit * 1./3. * PPneq.tr_sum();
							clTensorS3Dd PPneq_s = (PPneq - PPneqtr) * (1. + lamda) + PPneqtr * (1. + lamdab);
							clTensorS3Dd RoUUPF = RoUU + PPneq_s + FFbody2;
							for(int k = 0; k < c.lb.nL; k++)
							{
								const clVector3Dd C = c.lb.C[k];
								clTensorS3Dd CC = clTensorS3Dd(C) - Tunit * 1./3.;
#ifdef STOKES_FLOW
									double Fnew = c.lb.F0[k] * (Ro + 3. * (JF * C));
#else
									double Fnew = c.lb.F0[k] * PRoCorr * (Ro + 3. * (JF * C) + RoUUPF * CC * 4.5);
#endif
								vlb.F(k)(i,j,z) = Fnew / c.lb.Ro[M];
							}
#endif
						}
					}
				}
			}
		}

	}

	void collision_without_stress_LBGK()
	{
		//return;
		int i, inout;
		const double c13 = 1./3.;
		double PRoCorr, Jmag;
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
			for(i = 0; i < c.lb.nX; i++)
			{
				for(int j = 0; j < c.lb.nY; j++)
				{
					for(int z = 0; z < c.lb.nZ; z++)
					{						
						char M = vlb.M(i,j,z);
						//printf("\n M: %d\n", vlb.M(i,j,z));
						if(M != LB_SOLID)
						{
						#ifdef RAMP
						Fbody = c.lb.F[M] * c.lb.FM(i,j,z) * c.t.ramping();						
						#else
							#ifdef SQUARE
								Fbody = c.lb.F[M] * c.lb.FM(i,j,z) * c.t.square();
							#else
								Fbody = c.lb.F[M] * c.lb.FM(i,j,z);
								//printf("\n Fbody: %f, %f, %f \n", Fbody.x, Fbody.y, Fbody.z);
							#endif
						#endif

							//Fbody = vlb.G(i,j,z);
							double Ro = vlb.Ro(i,j,z);
							double Ro1 = 1./Ro;
							double RoM1 = 1./c.lb.Ro[M];
							inout = vlb.M.PBC(i,j,z);
							if(inout == 0)
							{
								PRoCorr = 1.;
							}
							else if(inout == 1)
							{
								PRoCorr = (1. + DELTA_RO * vlb.J.invel) / Ro;
							}
							else
							{
#ifndef OUTLET_INTERP
								PRoCorr = (1. - DELTA_RO * vlb.J.invel) / Ro;
#endif // OUTLET_INTERP
							}
							clVector3Dd  J = vlb.J(i,j,z); 
							#ifdef INOUTNORMAL
							if(inout != 0)
							{
								Jmag = J.r();
								J = c.ls.InoutNorms(inout) * Jmag;
								if(vlb.J.invel < 0.){
									J = J * (-1.);
								}
							}
#endif // INOUTNORMAL
							clTensorS3Dd RoUU = clTensorS3Dd(J) * Ro1;
							clTensorS3Dd FFbody2 = clTensorS3Dd(J*Ro1, Fbody);
							clVector3Dd JF = J + Fbody;
							clTensorS3Dd RoUUF = RoUU + FFbody2;
#ifdef BINARY_FLUID
							double N = vlb.N(i,j,z);
#endif // BINARY_FLUID
#ifdef LB_STOKES
							RoUUF = clTensorS3Dd(0.,0.,0.,0.,0.,0.);
#endif  //LB_STOKES
							for(int k = 0; k < c.lb.nL; k++)
							{
								const clVector3Dd C = c.lb.C[k];
								clTensorS3Dd CC = clTensorS3Dd(C) - Tunit * c13;
								double Fnew = c.lb.F0[k] * PRoCorr * (Ro + 3. * (JF * C)+ RoUUF * CC * 4.5);
								//vlb.F(k).ff(i,j,z) = Fnew * RoM1;
								vlb.F(k)(i,j,z) = Fnew * RoM1;
#ifdef BINARY_FLUID
								double Gnew = vlb.F(k)(i,j,z) * N / Ro;
								vlb.G(k)(i,j,z) = (vlb.G(k)(i,j,z) * (1. - LB_TAUN_NUMBER) + Gnew * LB_TAUN_NUMBER);
#endif // BINARY_FLUID
							}
						}
					}
				}
			}
		}
	}

//	void collision_without_stress_LBGK()
//	{
//		//return;
//		int i;
//#ifdef LB_USE_OPENMP
//#pragma omp parallel default(shared) private(i)
//#endif //LB_USE_OPENMP
//		{
//			clVector3Dd Fbody; 
//			clTensorS3Dd Tunit;
//			Tunit.unit();
//
//#ifdef LB_USE_OPENMP
//#pragma omp for
//#endif //LB_USE_OPENMP
//			for(i = 0; i < c.lb.nX; i++)
//			{
//				for(int j = 0; j < c.lb.nY; j++)
//				{
//					for(int z = 0; z < c.lb.nZ; z++)
//					{
//						char M = vlb.M(i,j,z);
//						if(M != LB_SOLID)
//						{
//							Fbody = c.lb.F[M] * c.lb.FM(i,j,z);
//							double Ro = vlb.Ro(i,j,z);
//
//							clVector3Dd  J = vlb.J(i,j,z); 
//							clTensorS3Dd RoUU = clTensorS3Dd(J) / Ro;
//							clTensorS3Dd FFbody2 = clTensorS3Dd(J/Ro, Fbody);
//							clVector3Dd JF = J + Fbody;
//							
//
//							double N = vlb.N(i,j,z);
//							clVector3Dd NU = J / Ro * N;
//							clTensorS3Dd NUU = RoUU / Ro * N;
//
//
//							double Fi = Ro * c.lb.Cs2 * log(Ro) - LBB_A/2*N*N + LBB_B/4*N*N*N*N;
//							double dFidRo = (log(Ro) + 1) / 3;
//							double dFidN = (-LBB_A + LBB_B*N*N) * N;
//
//							
//							clTensorS3Dd PP = vlb.PP(i,j,z);
//							clTensorS3Dd d2N = Tunit * vlb.N.d2(i,j,z);
//							double P0 = Ro * dFidRo + N * dFidN - Fi;
//							double A = ( P0 + LBB_KAPPA*N*vlb.N.d2(i,j,z).sum() 
//								- 0.5*LBB_KAPPA*(vlb.N.d(i,j,z).r2()) ) /c.lb.Cs2;
//							clTensorS3Dd At(A,A,A,A,A,A); 
//
//
//							clTensorS3Dd RoUUF = At + RoUU + FFbody2 + clTensorS3Dd(vlb.N.d(i,j,z));//?
//							
//
//							for(int k = 0; k < c.lb.nL; k++)
//							{
//								const clVector3Dd C = c.lb.C[k];
//								clTensorS3Dd CC = clTensorS3Dd(C) - Tunit * 1./3.;
//								double Fnew = c.lb.F0[k] * (Ro + 3. * (JF * C) + RoUUF * CC * 4.5);
//								vlb.F(k)(i,j,z) = Fnew / c.lb.Ro[M];
//								
//								double G = vlb.G(k)(i,j,z);
//								double GMu = G * (dFidN - LBB_KAPPA *vlb.N.d2(i,j,z).sum()) / c.lb.Cs2;
//
//								double Gnew = c.lb.F0[k] * (N + 3. * (NU * C) + NUU * CC * 4.5);
//								vlb.G(k)(i,j,z) = Gnew;
//
//	
//							}
//
//							//vlb.F(0)(i,j,z) = Ro / c.lb.Ro - Ftot;
//
//
//
//						}
//					}
//				}
//			}
//		}
//	}

	void collisionMRTLB()
	{
		int i;
#ifdef LB_USE_OPENMP
#pragma omp parallel default(shared) private(i)
#endif //LB_USE_OPENMP
		{

#ifdef LB_USE_OPENMP
#pragma omp for
#endif //LB_USE_OPENMP
			for(i = 0; i < c.lb.nX; i++)
			{
				for(int j = 0; j < c.lb.nY; j++)
				{
					for(int z = 0; z < c.lb.nZ; z++)
					{
						char M = vlb.M(i,j,z);
						if(M != LB_SOLID)
						{
							double lamda = c.lb.lamda[M];
							clVector3Di ii = clVector3Di(i,j,z);
							clVector3Dd Fbody = c.lb.F[M] * c.lb.FM(ii) / c.lb.Cs2;
							double Ro = vlb.Ro(ii);
							double Ro0 = c.lb.Ro[M];
							clVector3Dd  J = vlb.J(ii);

							double Mf[LB_NUMBER_CONNECTIONS];
							for(int k = 0; k < c.lb.nL; k++)
							{
								Mf[k] = 0.;
								double S = vlb.S(k, lamda);
								if(S == 0) 
									continue;

								for(int n = 0; n < c.lb.nL; n++)
									Mf[k] += c.lb.M[k][n] * vlb.F(n)(ii) * Ro0;

								Mf[k] = (Mf[k] - vlb.Meq(k, Ro, J)) * S;
							}

							for(int k = 0; k < c.lb.nL; k++)
							{
								double Fs = 0.;
								for(int n = 0; n < c.lb.nL; n++)
									Fs += c.lb.M1[k][n] * Mf[n];

								vlb.F(k)(i,j,z) += (Fs + Fbody * c.lb.C[k] * c.lb.F0[k]) / Ro0;
							}
						}
					}
				}
			}
		}
	}
};


#endif // !defined(AFX_CLMETHOD_H__INCLUDED_)
