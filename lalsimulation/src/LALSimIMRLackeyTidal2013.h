#ifndef _LALSIM_IMR_LACKEY_TIDAL_2013_H
#define _LALSIM_IMR_LACKEY_TIDAL_2013_H

void XLALLackeyTidal2013tidalPNAmplitudeCoefficient(
  double *C,
  const double eta,
  const double chi_BH,
  const double Lambda
);

double XLALLackeyTidal2013tidalCorrectionAmplitude(
  const double mf,
  const double C,
  const double eta,
  const double Lambda
);

// precompute a0, a1 and G which do not depend on frequency
void XLALLackeyTidal2013tidalPNPhaseCoefficients(
  double *a0,
  double *a1,
  double *G,
  const double eta,
  const double chi_BH,
  const double Lambda
);

// Implements Eq. 34 of Lackey et al
double XLALLackeyTidal2013tidalCorrectionPhase(
  const double mf,
  const double a0,
  const double a1,
  const double G,
  const double eta,
  const double Lambda
);

#endif /* _LALSIM_IMR_LACKEY_TIDAL_2013_H */
