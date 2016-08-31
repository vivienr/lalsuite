int XLALSimRingdownFD(
                                      COMPLEX16FrequencySeries **hlmmodetilde_out_plus,  /**< complex waveform for lm mode */
                                      COMPLEX16FrequencySeries **hlmmodetilde_out_cross,  /**< complex waveform for lm mode */
                                      REAL8 deltaF,                     /**< frequency resolution */
                                      REAL8 phi0,                   /**< initial phase of ringdown (rad) */
                                      const REAL8 fStart,                    /**< lowest GW frequency (Hz) of waveform generation */
                                      const REAL8 fEnd,                      /**< highest GW frequency (Hz) of waveform generation */
                                      REAL8 mass,                       /**< black hole mass (kg) */
                                      REAL8 eta,                   /**< symmetric mass ratio progenitors */
                                      REAL8 a,  /**< black hole dimensionless spin parameter */
                                      REAL8 chiEff,      /**< effective spin parameter for initial spins */
                                      REAL8 t0,  /**< ringdown starting time (offset from merger time) */
                                      REAL8 distance,           /**< distance to source (m) */
                                      REAL8 inclination,                /**< inclination of source's spin axis (rad) */
                                      LALSimInspiralTestGRParam *TGRParams   /** testing GR params - shifts dfreq and dtau */                             
);
REAL8 XLALSimRingdownFDAmplitudes(INT4 l, INT4 m, INT4 n, REAL8 eta, REAL8 chiEff);
COMPLEX16 XLALSimRingdownFitComplexOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a);
REAL8 XLALQNMOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a, REAL8 mtot);
REAL8 XLALQNMTau(UINT4 l, INT4 m, UINT4 n, REAL8 a, REAL8 mtot);

