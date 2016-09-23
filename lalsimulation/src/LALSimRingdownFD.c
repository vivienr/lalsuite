/* Frequency domain waveform for ringdown signal (TODO: more info)*/

#include <stdlib.h>
#include <lal/TimeDelay.h>
#include <gsl/gsl_vector.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/SphericalHarmonics.h>
#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimRingdownFD.h>

/*STATIC DECLARATIONS*/
static COMPLEX16 XLALSimSphericalHarmonicPlus(UINT4 l, INT4 m, REAL8 iota);
static COMPLEX16 XLALSimSphericalHarmonicCross(UINT4 l, INT4 m, REAL8 iota);
static void XLALShiftParams(LALSimInspiralTestGRParam *TGRParams, REAL8 *omega, REAL8 *tau); 

/*DEFINITIONS*/
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
                                      REAL8 distance,           /**< distance to source (m) */
                                      REAL8 inclination,                /**< inclination of source's spin axis (rad) */
                                      LALSimInspiralTestGRParam *TGRParams   /** testing GR params - shifts dfreq and dtau */                          
                                      )
{

  printf("\nCalled XLALSimRingdownFD with parameters:\n");
  printf("BH_mass = %f\n",mass/LAL_MSUN_SI);
  printf("BH_spin = %f\n",a);
  printf("chieff = %f\n",chiEff);
  printf("eta = %f\n",eta);
  printf("phi0 = %f\n",phi0);
  printf("costheta_jn = %f\n",cos(inclination));
  printf("logdistance = %f\n",log(distance/(LAL_PC_SI * 1.0e6)));
  printf("fStart = %f\n",fStart);
  printf("fEnd = %f\n",fEnd);
  printf("deltaF = %f\n",deltaF);

  /* switch to units where c=G=1, expressed in s */
  mass = mass*LAL_MTSUN_SI/LAL_MSUN_SI;
  REAL8 dist_sec = distance/LAL_C_SI;

  /*allocate memory*/
  UINT4 iEnd = (UINT4) ceil(fEnd / deltaF);
  LIGOTimeGPS tC = {0, 0};
  REAL8 shift = LAL_TWOPI*(tC.gpsSeconds + 1e-9 * tC.gpsNanoSeconds);
  XLALGPSAdd(&tC, -1 / deltaF);  /*coalesce at t=0*/
  COMPLEX16FrequencySeries *hlmmodetilde_plus = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, iEnd);
  COMPLEX16FrequencySeries *hlmmodetilde_cross = XLALCreateCOMPLEX16FrequencySeries("htilde: FD waveform", &tC, 0.0, deltaF, &lalStrainUnit, iEnd);

  /*Compute mode dependent quantities*/
  UINT4 l;
  INT4 m;
  UINT4 n=0;
  UINT4 i = 0;
  UINT4 iStart = (UINT4) ceil(fStart / deltaF);
  for (i = iStart; i < iEnd; i++) {
    REAL8 a=0;
    REAL8 b=0;
    for(l=0;l<5;l++){
     for(m=0;m<5;m++){
      REAL8 A = XLALSimRingdownFDAmplitudes(l, m, n, eta, chiEff);
      if(A!=0){
        COMPLEX16 Yplus = XLALSimSphericalHarmonicPlus(l, m, inclination);
        COMPLEX16 Ycross = XLALSimSphericalHarmonicCross(l, m, inclination);
        REAL8 omega = XLALQNMOmega(l, m, n, a, mass);
        REAL8 tau = XLALQNMTau( l, m, n, a, mass);
        //TODO: not checked
        // maybe make an other static function and loop over all possible shifts, name dfreq21, dfrac22, and so on
        if (TGRParams!=NULL){XLALShiftParams(TGRParams,&omega,&tau);}

        REAL8 f = i * deltaF;
        a += mass/dist_sec*A*Yplus*(cos(f*shift) - I*sin(f*shift))*(tau*((-1-I*2*f*LAL_PI*tau)*cos(m*phi0)-omega*tau*sin(m*phi0)))/(-1-I*4*f*LAL_PI*tau-omega*omega*tau*tau+4*f*f*LAL_PI*LAL_PI*tau*tau);
        a += -mass/dist_sec*A*Ycross*(cos(f*shift) - I*sin(f*shift))*tau*(omega*tau*cos(m*phi0)+(-1-I*2*f*LAL_PI*tau)*sin(m*phi0))/(1+tau*(omega*omega*tau-4*f*LAL_PI*(-I+f*LAL_PI*tau)));
      }
     }
    }
    hlmmodetilde_plus->data->data[i] = a;
    hlmmodetilde_cross->data->data[i] = b;
  }
  *hlmmodetilde_out_plus = hlmmodetilde_plus;
  *hlmmodetilde_out_cross = hlmmodetilde_cross;
  return 0;
}

static COMPLEX16 XLALSimSphericalHarmonicPlus(UINT4 l, INT4 m, REAL8 iota){
  COMPLEX16 Yplus = 0.0;
  Yplus = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) + (1.0 - 2.0*(l % 2))*XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, -m);
  return Yplus;
}

static COMPLEX16 XLALSimSphericalHarmonicCross(UINT4 l, INT4 m, REAL8 iota){
  COMPLEX16 Ycross = 0.0;
  Ycross = XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, m) - (1.0 - 2.0*(l % 2))*XLALSpinWeightedSphericalHarmonic(iota, 0.0, -2, l, -m);
  return Ycross;
}


/**
 * Calculates the amplitudes for the QNM l, m, n=0 for
 * a given symmetric mass-ratio eta of the initial binary.
 * Based on I. Kamaretsos, M. Hannam, and B.S. Sathyaprakash, Phys. Rev. Lett. 109, 141102 (2012
 **/
REAL8 XLALSimRingdownFDAmplitudes(INT4 l, INT4 m, INT4 n, REAL8 eta, REAL8 chiEff){
  //TODO: include spin and overtones
  
  REAL8 A = 0.8639*eta;
  if (l==2 && m==2){ A *= 1.0; }
  else if (l==2 && abs(m)==1){ A *= 0.43*(sqrt(1.0 - 4.*eta) - chiEff); }
  //else if (l==2 && abs(m)==1){ A *= 0.52*pow(1.0 - 4.*eta, 0.71); }        
  else if (l==3 && abs(m)==3){ A *= 0.44*pow(1.0 - 4.*eta, 0.45); }
  else if (l==3 && abs(m)==2){ A *= 3.69*(eta - 0.2)*(eta - 0.2) + 0.053; }
  else if (l==4 && abs(m)==4){ A *= 5.41*((eta - 0.22)*(eta - 0.22) + 0.04); }
  else A = 0.0*n;
  return A;
}

/**
 * Computes the complex dimensionless eigen-frequency for the
 * quasinormal mode (l,m), given a dimensionless spin a \in [-1,1].
 * This is based on the interpolation formula from Berti et al 2006
 * [arXiv: gr-qc/0512160].
 * The dimensionful angular frequency is Re(omega)/M.
 * The dimensionful damping time is M/Im(omega).
 *
 * \note The dimensionless spin assumes values between -1 and 1.
 * The complex omega is M*(w_lmn+i/tau_lmn)
 */
COMPLEX16 XLALSimRingdownFitComplexOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a){
	COMPLEX16 omega=0.0;
	if (fabs(a) > 1.){
		fprintf(stderr, "ERROR: Dimensionless spin magnitude larger than 1! Aborting...\n");
		exit(-1);
	}
	switch (l){
	case 2:
		switch (abs(m)){
		case 2:
			omega = (1.5251 - 1.1568*pow(1.e0-a,0.1292));
			omega += I*omega/(2.*(0.7000 + 1.4187*pow(1.e0-a,-0.4990)));
			break;
		case 1:
			omega = (0.6000 - 0.2339*pow(1.e0-a,0.4175));
			omega += I*omega/(2.*(-0.3000 + 2.3561*pow(1.e0-a,-0.2277)));
			break;
		default:
			break;
		}
		break;
	case 3:
		switch (abs(m)){
		case 3:
			omega = (1.8956 - 1.3043*pow(1.e0-a,0.1818));
			omega += I*omega/(2.*(0.9000 + 2.3430*pow(1.e0-a,-0.4810)));
			break;
		case 2:
			omega = (1.1481 - 0.5552*pow(1.e0-a,0.3002));
			omega += I*omega/(2.*(0.8313 + 2.3773*pow(1.e0-a,-0.3655)));
			break;
		default:
			break;
		}
		break;
	case 4:
		switch (abs(m)){
		case 4:
			omega = (2.3000 - 1.5056*pow(1.e0-a,0.2244));
			omega += I*omega/(2.*(1.1929 + 3.1191*pow(1.e0-a,-0.4825)));
			break;
		default:
			break;
		}
		break;
	default:
		break;
	}
	if (omega == 0.0){
		fprintf(stderr, "ERROR: Invalid mode l = %u, m = %d, n = %u\n", l, m, n);
		exit(-1);
	}
	return omega;
}

/** 
 * Frequency in rad/sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMOmega(UINT4 l, INT4 m, UINT4 n, REAL8 a, REAL8 mtot){
        COMPLEX16 omega = XLALSimRingdownFitComplexOmega(l,m,n,a);
	return creal(omega)/mtot;
}

/** 
 * Damping time in sec given dimensionless complex frequency and M
 **/
REAL8 XLALQNMTau(UINT4 l, INT4 m, UINT4 n, REAL8 a, REAL8 mtot){
        COMPLEX16 omega = XLALSimRingdownFitComplexOmega(l,m,n,a);
	REAL8 tau = 0;
	tau = mtot/cimag(omega);
	return tau;
}

void XLALShiftParams(LALSimInspiralTestGRParam *TGRParams, REAL8 *omega, REAL8 *tau){
        REAL8 dtau22, domega22, dtau33, domega33;
        if ( !XLALSimInspiralTestGRParamExists(TGRParams, "domega22") ){
          domega22=XLALSimInspiralGetTestGRParam(TGRParams, "domega22");
          *omega = *omega*(1+domega22);
        }
        if ( !XLALSimInspiralTestGRParamExists(TGRParams, "domega33") ){
          domega33=XLALSimInspiralGetTestGRParam(TGRParams, "domega33");
          *omega = *omega*(1+domega33);
        }
        if ( !XLALSimInspiralTestGRParamExists(TGRParams, "dtau22") ){
          dtau22=XLALSimInspiralGetTestGRParam(TGRParams, "dtau22");
          *tau=*tau*(1+dtau22);
        }    
        if ( !XLALSimInspiralTestGRParamExists(TGRParams, "dtau33") ){
          dtau33=XLALSimInspiralGetTestGRParam(TGRParams, "dtau33");
          *tau=*tau*(1+dtau33);
        }
}




