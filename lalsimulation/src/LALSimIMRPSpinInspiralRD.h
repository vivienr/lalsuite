#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>

#ifndef _LALSIMIMRPHENSPIN_H
#define _LALSIMIMRPHENSPIN_H

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define LAL_NUM_PST4_VARIABLES 12
#define LAL_PST4_ABSOLUTE_TOLERANCE 1.e-12
#define LAL_PST4_RELATIVE_TOLERANCE 1.e-12

/* use error codes above 1024 to avoid conflicts with GSL */
#define LALSIMINSPIRAL_PHENSPIN_TEST_ENERGY	 1025
#define LALSIMINSPIRAL_PHENSPIN_TEST_OMEGADOT    1026
#define LALSIMINSPIRAL_PHENSPIN_TEST_FREQBOUND   1027
#define LALSIMINSPIRAL_PHENSPIN_TEST_OMEGANAN	 1028
#define LALSIMINSPIRAL_PHENSPIN_TEST_OMEGAMATCH  1029
#define LALSIMINSPIRAL_PHENSPIN_TEST_OMEGANONPOS 1031

/**
 * Comments here
 */

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

int XLALSimIMRPhenSpinGenerateQNMFreq(COMPLEX16Vector *modefreqs,
				      UINT4 l,
				      INT4  m,
				      REAL8 finalMass,
				      REAL8 finalSpin,
				      REAL8 totalMass);

#endif /* _LALSIMIMRPHENSPIN_H */
