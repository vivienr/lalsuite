#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#ifndef _LALSIMSPINTAYLOR_H
#define _LALSIMSPINTAYLOR_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_ST4_TEST_ENERGY                  1025
#define LALSIMINSPIRAL_ST4_TEST_OMEGADOT                1026
#define LALSIMINSPIRAL_ST4_TEST_COORDINATE              1027
#define LALSIMINSPIRAL_ST4_TEST_OMEGANAN                1028
#define LALSIMINSPIRAL_ST4_TEST_FREQBOUND               1029
#define LALSIMINSPIRAL_ST4_DERIVATIVE_OMEGANONPOS       1030
#define LALSIMINSPIRAL_PST4_TEST_OMEGAMATCH              1031

/* Number of variables used for precessing waveforms */
#define LAL_NUM_PST4_VARIABLES 12
#define LAL_NUM_ST4_VARIABLES 14
/* absolute and relative tolerance for adaptive Runge-Kutta ODE integrator */
/* 1.e-06 is too large for end of 1.4--1.4 M_sun BNS inspiral */
/* (phase difference at end will be ~10% of GW cycle). */
/* 1.e-12 is used so last data point isn't nan for 6PN tidal, */
/* since larger values probably cause larger step sizes. */
#define LAL_ST4_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_ST4_RELATIVE_TOLERANCE 1.e-12

/**
 * Struct containing all of the non-dynamical coefficients needed
 * to evolve a spinning, precessing binary and produce a waveform.
 */

typedef struct tagLALSimInspiralSpinTaylorT4Coeffs
{
  REAL8 M; 			   // total mass in seconds
  REAL8 eta; 			   // symmetric mass ratio
  REAL8 m1Bym2;                      // ratio m1/m2
  REAL8 m1ByM;                       // ratio m1/M
  REAL8 m2ByM;                       // ratio m2/M
  REAL8 dmByM;                       // ratio dm/M
  REAL8 wdotnewt;                    //leading order coefficient of wdot = \f$\dot{\omega}\f$
  REAL8 wdotcoeff[LAL_MAX_PN_ORDER]; // coeffs. of PN corrections to wdot
  REAL8 wdotlogcoeff; 		   // coefficient of log term in wdot
  REAL8 wdotSO15s1, wdotSO15s2; 	   // non-dynamical 1.5PN SO corrections
  REAL8 wdotSS2,wdotSSO2;
  REAL8 wdotSSselfS1,wdotSSselfS2;      // non-dynamical 2PN self-spin correction
  REAL8 wdotSSselfS1L,wdotSSselfS2L;    // non-dynamical 2PN self-spin correction
  REAL8 wdotQM2S1,wdotQM2S1L;
  REAL8 wdotQM2S2,wdotQM2S2L;
  REAL8 wdotSO25s1,wdotSO25s2; 	   // non-dynamical 2.5PN SO corrections
  REAL8 wdotSO3s1,wdotSO3s2; 	   // non-dynamical 2.5PN SO corrections
  REAL8 Enewt;                       // coeffs. of PN corrections to energy
  REAL8 Ecoeff[LAL_MAX_PN_ORDER];    // coeffs. of PN corrections to energy
  REAL8 ESO15s1, ESO15s2; 	   // non-dynamical 1.5PN SO corrections
  REAL8 ESS2,ESSO2; 		   // non-dynamical 2PN SS correction
  REAL8 ESelfSSO2s1,ESelfSSO2s2; 	   // non-dynamical 2PN self-spin correction
  REAL8 ESelfSS2s1,ESelfSS2s2; 	           // non-dynamical 2PN self-spin correction
  REAL8 ESO25s1, ESO25s2; 	   // non-dynamical 2.5PN SO corrections 
  REAL8 LNhatSO15s1, LNhatSO15s2;    // non-dynamical 1.5PN SO corrections
  REAL8 LNhatSS2; 		   // non-dynamical 2PN SS correction
  REAL8 S1dot15,S2dot15;             // non-dynamical leading SO term in S1,2dot
  REAL8 S1dot25,S2dot25;             // non-dynamical Next to leading SO correction to S1,2dot
  REAL8 wdottidal5pn;		   // leading order tidal correction 
  REAL8 wdottidal6pn;	           // next to leading order tidal correction
  REAL8 Etidal5pn;	           // leading order tidal correction to energy
  REAL8 Etidal6pn;                   // next to leading order tidal correction to energy
  REAL8 EQM2S1; ///< non-dynamical S1^2 2PN quadrupole-monopole correction
  REAL8 EQM2S1L;///< non-dynamical (S1.L)^2 2PN quadrupole-monopole correction
  REAL8 EQM2S2; ///< non-dynamical S2^2 2PN quadrupole-monopole correction
  REAL8 EQM2S2L;///< non-dynamical (S2.L)^2 2PN quadrupole-monopole correction
  REAL8 fStart; 			   // starting GW frequency of integration
  REAL8 fEnd; 			   // ending GW frequency of integration
  REAL8 dt;                          // sampling in seconds
} LALSimInspiralSpinTaylorT4Coeffs;

typedef struct tagLALSimInspiralInclAngle {
  REAL8 cHi;
  REAL8 sHi;
  REAL8 ci;
  REAL8 si;
  REAL8 ci2;
  REAL8 si2;
  REAL8 cHi2;
  REAL8 sHi2;
  REAL8 cHi3;
  REAL8 sHi3;
  REAL8 cHi4;
  REAL8 sHi4;
  REAL8 cHi5;
  REAL8 sHi5;
  REAL8 cHi6;
  REAL8 sHi6;
  REAL8 cHi8;
  REAL8 sHi8;
  REAL8 cDi;
  REAL8 sDi;
} LALSimInspiralInclAngle;

/**
 * Comments here
 **/

int XLALSimSpinInspiralFillL2Modes(COMPLEX16Vector *hL2,
				   REAL8 v,
				   REAL8 eta,
				   REAL8 dm,
				   REAL8 Psi,
				   REAL8 alpha,
				   LALSimInspiralInclAngle *an);

/**
 * Comments here
 */

int XLALSimSpinInspiralFillL3Modes(COMPLEX16Vector *hL3,
				   REAL8 v,
				   REAL8 eta,
				   REAL8 dm,
				   REAL8 Psi,
				   REAL8 alpha,
				   LALSimInspiralInclAngle *an);

/**
 * Comments here
 */

int XLALSimSpinInspiralFillL4Modes(COMPLEX16Vector* hL4,
				   UNUSED REAL8 v,
				   REAL8 eta,
				   UNUSED REAL8 dm,
				   REAL8 Psi,
				   REAL8 alpha,
				   LALSimInspiralInclAngle *an);

/**
 * Comments here
 **/

int XLALGenerateWaveDerivative (REAL8Vector *dwave,
				REAL8Vector *wave,
				REAL8 dt);

/**
 * Comments here
 **/

int XLALSimInspiralComputeInclAngle(REAL8 ciota, LALSimInspiralInclAngle *angle);

/**
 * Comments here
 **/

int XLALSimIMRPSpinInspiralSetAxis(REAL8 m1,
				   REAL8 m2,
				   REAL8 iota,
				   REAL8 *yinit,
				   LALSimInspiralFrameAxis axisChoice);

#endif /* _LALSIMIMRPHENSPIN_H */
