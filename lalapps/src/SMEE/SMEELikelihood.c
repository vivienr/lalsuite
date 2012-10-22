/* 
 *  LALInferenceLikelihood.c:  Bayesian Followup likelihood functions
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
 *
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <SMEELikelihood.h>
#include <lal/LALInference.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/Sequence.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>


#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


/** Internal functions which are not exported */
static UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series);
static UINT4 LIGOTimeGPSToNearestIndex(const LIGOTimeGPS *tm, const REAL8TimeSeries *series) {
  REAL8 dT = XLALGPSDiff(tm, &(series->epoch));

  return (UINT4) (round(dT/series->deltaT));
}

/* The time-domain weight corresponding to the (two-sided) noise power
   spectrum S(f) is defined to be:

   s(tau) == \int_{-infty}^\infty df \frac{\exp(2 \pi i f tau)}{S(f)}

*/

static UINT4 nextPowerOfTwo(const UINT4 n);
static UINT4 nextPowerOfTwo(const UINT4 n) {
  UINT4 np2 = 1;
  
  for (np2 = 1; np2 < n; np2 *= 2) ; /* Keep doubling until >= n */

  return np2;
}

static void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
static void padREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  if (padded->length < data->length) {
    fprintf(stderr, "padREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, padded->length*sizeof(padded->data[0]));
  memcpy(padded->data, data->data, data->length*sizeof(data->data[0]));
}

static void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data);
static void padWrappedREAL8Sequence(REAL8Sequence *padded, const REAL8Sequence *data) {
  UINT4 i;
  UINT4 np = padded->length;
  UINT4 nd = data->length;

  if (np < nd) {
    fprintf(stderr, "padWrappedREAL8Sequence: padded sequence too short (in %s, line %d)", __FILE__, __LINE__);
    exit(1);
  }

  memset(padded->data, 0, np*sizeof(padded->data[0]));

  padded->data[0] = data->data[0];
  for (i = 1; i <= (nd-1)/2; i++) {
    padded->data[i] = data->data[i]; /* Positive times/frequencies. */
    padded->data[np-i] = data->data[nd-i]; /* Wrapped, negative times/frequencies. */
  }
  if (nd % 2 == 0) { /* If even, take care of singleton positive frequency. */
    padded->data[nd/2] = data->data[nd/2];
  }
}


REAL8 LALInferenceSMEEFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, 
                              LALInferenceTemplateFunction *template)
/***************************************************************/
/* (log-) likelihood function.                                 */
/* Returns the non-normalised logarithmic likelihood.          */
/* Slightly slower but cleaner than							   */
/* UndecomposedFreqDomainLogLikelihood().          `		   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/
{
 
  REAL8 loglikeli, totalChiSquared=0.0;
  LALInferenceIFOData *ifoPtr=data;
  COMPLEX16Vector *freqModelResponse=NULL;
 
  /* loop over data (different interferometers): */
  while (ifoPtr != NULL) {
    
    ifoPtr->loglikelihood = 0.0;
//fprintf(stdout, "working \n");
	if(freqModelResponse==NULL)
		freqModelResponse= XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		freqModelResponse= XLALResizeCOMPLEX16Vector(freqModelResponse, ifoPtr->freqData->data->length);
	//fprintf(stdout, "working \n");
	/*compute the response*/
	LALInferenceComputeSMEEFreqDomainResponse(currentParams, ifoPtr, template, freqModelResponse);
	
	/*if(residual==NULL)
		residual=XLALCreateCOMPLEX16Vector(ifoPtr->freqData->data->length);
	else
		residual=XLALResizeCOMPLEX16Vector(residual, ifoPtr->freqData->data->length);
	
	COMPLEX16VectorSubtract(residual, ifoPtr->freqData->data, freqModelResponse);
	totalChiSquared+=ComputeFrequencyDomainOverlap(ifoPtr, residual, residual); 
	*/
        REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data)
          -2.0*LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, freqModelResponse)
          +LALInferenceComputeFrequencyDomainOverlap(ifoPtr, freqModelResponse, freqModelResponse);
	totalChiSquared+=temp;
        ifoPtr->loglikelihood -= 0.5*temp;
    
    ifoPtr = ifoPtr->next;
  }
  
  loglikeli = -0.5   *totalChiSquared; // note (again): the log-likelihood is unnormalised!
  XLALDestroyCOMPLEX16Vector(freqModelResponse);
 
  return(loglikeli);
}

void LALInferenceComputeSMEEFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr, 
                              LALInferenceTemplateFunction * template UNUSED, COMPLEX16Vector *freqWaveform)
/***************************************************************/
/* Frequency-domain single-IFO response computation.           */
/* Computes response for a given template.                     */
/* Will re-compute template only if necessary                  */
/* (i.e., if previous, as stored in data->freqModelhCross,     */
/* was based on different parameters or template function).    */
/* Carries out timeshifting for a given detector               */
/* and projection onto this detector.                          */
/* Result stored in freqResponse, assumed to be correctly      */
/* initialized												   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, >0)                      */
/*   - "time"            (REAL8, GPS sec.)                     */
/***************************************************************/							  
{
  //static int timeDomainWarning = 0;

	double ra, dec, psi, distMpc, gmst;
	//dataPtr->freqWaveform=&get_model;
	gsl_matrix *PCvectorsreal = NULL;
	gsl_matrix *PCvectorsimag = NULL;
	//get_model(currentParams, dataPtr, freqWaveform);
	INT4 ph;
	REAL8 beta1, beta2, beta3, beta4, beta5, beta6, beta7;
	double GPSdouble;
	double timeTmp;
	LIGOTimeGPS GPSlal;
	double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
	double timeshift;  /* time shift (not necessarily same as above)                   */
	double deltaT, deltaF, f, re, im;
	double twopit;
	//int A, B;
	//int different;
	LALInferenceVariables intrinsicParams;
	LALStatus status;
	memset(&status,0,sizeof(status));
	
	double Fplus, Fcross;
	double FplusScaled, FcrossScaled;
	//REAL8Vector *plainTemplate = NULL;
	//REAL8Vector *PCs = NULL;
	REAL8 plainTemplatereal, plainTemplateimag;
	UINT4 i;
	//REAL8 mc;
	/* Fill in derived parameters if necessary */
	//if(LALInferenceCheckVariable(currentParams,"logdistance")){
	//	distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
	//	LALInferenceAddVariable(currentParams,"distance",&distMpc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	//}

	/*if(LALInferenceCheckVariable(currentParams,"logmc")){
		mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
		LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
	}*/
		
	
	/* determine source's sky location & orientation parameters: */
	ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
	dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
	psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
	GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
	//distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */
	//ra=LAL_PI;
	//dec=0.0;
	//psi=LAL_PI/2;
	distMpc=1;
	//GPSdouble=944697615.0;
	/* figure out GMST: */
	//XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5));
	XLALGPSSetREAL8(&GPSlal, GPSdouble);
	//UandA.units    = MST_RAD;
	//UandA.accuracy = LALLEAPSEC_LOOSE;
	//LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA);
	gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
	intrinsicParams.head      = NULL;
	intrinsicParams.dimension = 0;
	LALInferenceCopyVariables(currentParams, &intrinsicParams);
	LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
	LALInferenceRemoveVariable(&intrinsicParams, "declination");
	LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
	LALInferenceRemoveVariable(&intrinsicParams, "time");
	//LALInferenceRemoveVariable(&intrinsicParams, "distance");	
	
// TODO: add pointer to template function here.
	// (otherwise same parameters but different template will lead to no re-computation!!)
      
	/* The parameters the response function can handle by itself     */
    /* (and which shouldn't affect the template function) are        */
    /* sky location (ra, dec), polarisation and signal arrival time. */
    /* Note that the template function shifts the waveform to so that*/
	/* t_c corresponds to the "time" parameter in                    */
	/* IFOdata->modelParams (set, e.g., from the trigger value).     */
    
    /* Compare parameter values with parameter values corresponding  */
    /* to currently stored template; ignore "time" variable:         */
   //if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
     // timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
     // LALInferenceRemoveVariable(dataPtr->modelParams, "time");
   // }
    /*else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);*/
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */
  
   // if (different) { /* template needs to be re-computed: */
     /*LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
		  if (!timeDomainWarning) {
			  timeDomainWarning = 1;
			  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
		  }
		  LALInferenceExecuteFT(dataPtr);*/
		  /* note that the dataPtr->modelParams "time" element may have changed here!! */
		  /* (during "template()" computation)                                      
		}
    }*/
    
    //else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      timeTmp=GPSdouble;
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    //}

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                            --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross, dataPtr->detector->response,
			     ra, dec, psi, gmst);
		 
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.) */

    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    
    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }
	
	deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    
	beta1        = *(REAL8*) LALInferenceGetVariable(currentParams, "beta1"); /* radian      */
	beta2      = *(REAL8*) LALInferenceGetVariable(currentParams, "beta2");    /* radian      */
	beta3        = *(REAL8*) LALInferenceGetVariable(currentParams, "beta3"); /* radian      */
	beta4      = *(REAL8*) LALInferenceGetVariable(currentParams, "beta4");    /* radian      */
	beta5      = *(REAL8*) LALInferenceGetVariable(currentParams, "beta5");    /* radian      */
	beta6        = *(REAL8*) LALInferenceGetVariable(currentParams, "beta6"); /* radian      */
	beta7      = *(REAL8*) LALInferenceGetVariable(currentParams, "beta7");    /* radian      */
	
	PCvectorsimag = *(gsl_matrix**) LALInferenceGetVariable(dataPtr->modelParams, "PCvectorsimag");
	PCvectorsreal = *(gsl_matrix**) LALInferenceGetVariable(dataPtr->modelParams, "PCvectorsreal");
	ph = *(INT4**) LALInferenceGetVariable(dataPtr->modelParams, "phase");
	//fprintf(stdout, "working \n");
	if(freqWaveform->length!=dataPtr->freqModelhPlus->data->length){
		printf("fW%d data%d\n", freqWaveform->length, dataPtr->freqModelhPlus->data->length);
		printf("Error!  Frequency data vector must be same length as original data!\n");
		exit(1);
	}
#ifdef DEBUG
FILE* file=fopen("TempSignal.dat", "w");	
#endif
	for(i=0; i<freqWaveform->length; i++){
		/* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
		plainTemplatereal = FplusScaled * (beta1 * gsl_matrix_get(PCvectorsreal, i, 0)  + beta2 * gsl_matrix_get(PCvectorsreal, i, 1)
		+ beta3 * gsl_matrix_get(PCvectorsreal, i, 2)  + beta4 * gsl_matrix_get(PCvectorsreal, i, 3)
		+ beta5 * gsl_matrix_get(PCvectorsreal, i, 4)  + beta6 * gsl_matrix_get(PCvectorsreal, i, 5)
		+ beta7 * gsl_matrix_get(PCvectorsreal, i, 6));
		/*plainTemplatereal = (beta1 * gsl_matrix_get(PCvectorsreal, i, 0)  + beta2 * gsl_matrix_get(PCvectorsreal, i, 1)
		+ beta3 * gsl_matrix_get(PCvectorsreal, i, 2)  + beta4 * gsl_matrix_get(PCvectorsreal, i, 3)
		+ beta5 * gsl_matrix_get(PCvectorsreal, i, 4)  + beta6 * gsl_matrix_get(PCvectorsreal, i, 5)
		+ beta7 * gsl_matrix_get(PCvectorsreal, i, 6));*/
		plainTemplateimag = FplusScaled * (beta1 * gsl_matrix_get(PCvectorsimag, i, 0)  + beta2 * gsl_matrix_get(PCvectorsimag, i, 1)
		+ beta3 * gsl_matrix_get(PCvectorsimag, i, 2)  + beta4 * gsl_matrix_get(PCvectorsimag, i, 3)
		+ beta5 * gsl_matrix_get(PCvectorsimag, i, 4)  + beta6 * gsl_matrix_get(PCvectorsimag, i, 5)
		+ beta7 * gsl_matrix_get(PCvectorsimag, i, 6));
		/*plainTemplateimag = (beta1 * gsl_matrix_get(PCvectorsimag, i, 0)  + beta2 * gsl_matrix_get(PCvectorsimag, i, 1)
		+ beta3 * gsl_matrix_get(PCvectorsimag, i, 2)  + beta4 * gsl_matrix_get(PCvectorsimag, i, 3)
		+ beta5 * gsl_matrix_get(PCvectorsimag, i, 4)  + beta6 * gsl_matrix_get(PCvectorsimag, i, 5)
		+ beta7 * gsl_matrix_get(PCvectorsimag, i, 6));*/

		/* do time-shifting...             */
		/* (also un-do 1/deltaT scaling): */
		f = ((double) i) * deltaF;
		/* real & imag parts of  exp(-2*pi*i*f*deltaT): */
		re = cos(twopit * f);
		//re = cos(0 * f);
		im = - sin(twopit * f);
		//im = - sin(0 * f);
		if(ph==0){
		//freqWaveform->data[i].re= (plainTemplatereal*re - plainTemplateimag*im);
		//freqWaveform->data[i].im= (plainTemplatereal*im + plainTemplateimag*re);
		freqWaveform->data[i].re= (plainTemplatereal)*sqrt(re*re +im*im);
		freqWaveform->data[i].im= (plainTemplateimag*re); }
		else if(ph==1){
		  freqWaveform->data[i].re= (plainTemplatereal*re);
		  freqWaveform->data[i].im= (plainTemplateimag*re); }
		  
		
#ifdef DEBUG
		fprintf(file, "%lg %lg \t %lg\n", f, freqWaveform->data[i].re, freqWaveform->data[i].im);
#endif
	}
#ifdef DEBUG
fclose(file);
#endif
	
	LALInferenceAddVariable(dataPtr->modelParams,"PCvectorsreal",&PCvectorsreal, LALINFERENCE_gslMatrix_t,
					LALINFERENCE_PARAM_FIXED);
	LALInferenceAddVariable(dataPtr->modelParams,"PCvectorsimag",&PCvectorsimag, LALINFERENCE_gslMatrix_t,
					LALINFERENCE_PARAM_FIXED);
	LALInferenceDestroyVariables(&intrinsicParams);
	
}


REAL8 LALInferenceSMEENullLogLikelihood(LALInferenceIFOData *data)
/*Identical to FreqDomainNullLogLikelihood                        */
{      
	REAL8 loglikeli, totalChiSquared=0.0;
	LALInferenceIFOData *ifoPtr=data;
	FILE *noise=fopen("noisefile2.dat", "w");
	int i=0;
	for (i=1; i<=6143; ++i){  	
	 fprintf(noise, "%lg \n",data->oneSidedNoisePowerSpectrum->data->data[i] );
	//fprintf(stderr, "%lg \n",data->oneSidedNoisePowerSpectrum->data->data[i] ); 
	 
		 
  }
fclose(noise);
	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
	   
          REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
	 
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	
	loglikeli = -0.5  * totalChiSquared; // note (again): the log-likelihood is unnormalised!
	
	return(loglikeli);
}

REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * dataPtr,
                                    COMPLEX16Vector * freqData1, 
                                    COMPLEX16Vector * freqData2)
{
  
  if (dataPtr==NULL || freqData1 ==NULL || freqData2==NULL){
  	XLAL_ERROR_REAL8(XLAL_EFAULT); 
  	}
  	
  int lower, upper, i;
  double deltaT, deltaF;
  double norm;
  double overlap=0.0;
  
  /* determine frequency range & loop over frequency bins: */
  deltaT = dataPtr->timeData->deltaT;
  
  norm = (deltaT /(((double)dataPtr->timeData->data->length)));
  //norm = (deltaT /(12276));
  norm = norm;
  deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
 //fprintf(stdout,"the number= %lg \n",norm);
  lower = ceil(dataPtr->fLow / deltaF);
  upper = floor(dataPtr->fHigh / deltaF);
 
  for (i=lower; i<=upper; ++i){  	
    
    overlap  += ((4.0*norm*(freqData1->data[i].re*freqData2->data[i].re+freqData1->data[i].im*freqData2->data[i].im)) 
                 / (dataPtr->oneSidedNoisePowerSpectrum->data->data[i]));
		
	//fprintf(stderr, "%lg \n",dataPtr->oneSidedNoisePowerSpectrum->data->data[i] ); 
	// overlap  += ((4.0*deltaF*(fabs(freqData1->data[i].re*freqData2->data[i].re))) 
          //       / fabs(dataPtr->oneSidedNoisePowerSpectrum->data->data[i]));
		 
  }

  return overlap;
  
}

REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data)
/*Identical to FreqDomainNullLogLikelihood                        */
{      
	REAL8 loglikeli, totalChiSquared=0.0;
	LALInferenceIFOData *ifoPtr=data;
	
	/* loop over data (different interferometers): */
	while (ifoPtr != NULL) {
          ifoPtr->nullloglikelihood = 0.0;
	   
          REAL8 temp = LALInferenceComputeFrequencyDomainOverlap(ifoPtr, ifoPtr->freqData->data, ifoPtr->freqData->data);
	 
          totalChiSquared+=temp;
          ifoPtr->nullloglikelihood -= 0.5*temp;
		ifoPtr = ifoPtr->next;
	}
	
	loglikeli = -0.5 *totalChiSquared; // note (again): the log-likelihood is unnormalised!
	
	return(loglikeli);
}

/***************************************************************/
/* Student-t (log-) likelihood function                        */
/* as described in Roever/Meyer/Christensen (2011):            */
/*   "Modelling coloured residual noise                        */
/*   in gravitational-wave signal processing."                 */
/*   Classical and Quantum Gravity, 28(1):015010.              */
/*   http://dx.doi.org/10.1088/0264-9381/28/1/015010           */
/*   http://arxiv.org/abs/0804.3853                            */
/* Returns the non-normalised logarithmic likelihood.          */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`currentParams') parameters are:                  */
/*   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       */
/*   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  */
/*   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        */
/*   - "distance"        (REAL8, Mpc, > 0)                     */
/*   - "time"            (REAL8, GPS sec.)                     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* This function is essentially the same as the                */
/* "UndecomposedFreqDomainLogLikelihood()" function.           */
/* The additional parameter to be supplied is the (REAL8)      */
/* degrees-of-freedom parameter (nu) for each Ifo.             */
/* The additional "df" argument gives the corresponding        */
/* d.f. parameter for each element of the "*data" list.        */
/* The names of "df" must match the "->name" slot of           */
/* the elements of "data".                                     */
/*                                                             */
/* (TODO: allow for d.f. parameter to vary with frequency,     */
/*        i.e., to be a set of vectors corresponding to        */
/*        frequencies)                                         */
/***************************************************************/

REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data, 
                                      LALInferenceTemplateFunction *template)
{
  static int timeDomainWarning = 0;
  double Fplus, Fcross;
  double FplusScaled, FcrossScaled;
  double diffRe, diffIm, diffSquared;
  double dataReal, dataImag;
  REAL8 loglikeli;
  REAL8 plainTemplateReal, plainTemplateImag;
  REAL8 templateReal, templateImag;
  int i, lower, upper;
  LALInferenceIFOData *dataPtr;
  double ra, dec, psi, distMpc, gmst,mc;
  double GPSdouble;
  LIGOTimeGPS GPSlal;
  double chisquared;
  double timedelay;  /* time delay b/w iterferometer & geocenter w.r.t. sky location */
  double timeshift;  /* time shift (not necessarily same as above)                   */
  double deltaT, FourDeltaToverN, deltaF, twopit, f, re, im, singleFreqBinTerm;
  double degreesOfFreedom, nu;
  double timeTmp;
  int different;
  LALStatus status;
  memset(&status,0,sizeof(status));
  LALInferenceVariables intrinsicParams;
  
  /* Fill in derived parameters if necessary */
  if(LALInferenceCheckVariable(currentParams,"logdistance")){
    distMpc=exp(*(REAL8 *) LALInferenceGetVariable(currentParams,"logdistance"));
    LALInferenceAddVariable(currentParams,"distance",&distMpc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }

  if(LALInferenceCheckVariable(currentParams,"logmc")){
    mc=exp(*(REAL8 *)LALInferenceGetVariable(currentParams,"logmc"));
    LALInferenceAddVariable(currentParams,"chirpmass",&mc,LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_OUTPUT);
  }
  /* determine source's sky location & orientation parameters: */
  ra        = *(REAL8*) LALInferenceGetVariable(currentParams, "rightascension"); /* radian      */
  dec       = *(REAL8*) LALInferenceGetVariable(currentParams, "declination");    /* radian      */
  psi       = *(REAL8*) LALInferenceGetVariable(currentParams, "polarisation");   /* radian      */
  GPSdouble = *(REAL8*) LALInferenceGetVariable(currentParams, "time");           /* GPS seconds */
  distMpc   = *(REAL8*) LALInferenceGetVariable(currentParams, "distance");       /* Mpc         */

  /* figure out GMST: */
  /* XLALINT8NSToGPS(&GPSlal, floor(1e9 * GPSdouble + 0.5)); */
  XLALGPSSetREAL8(&GPSlal, GPSdouble);
  /* UandA.units    = MST_RAD;                       */
  /* UandA.accuracy = LALLEAPSEC_LOOSE;              */
  /* LALGPStoGMST1(&status, &gmst, &GPSlal, &UandA); */
  gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
  intrinsicParams.head      = NULL;
  intrinsicParams.dimension = 0;
  LALInferenceCopyVariables(currentParams, &intrinsicParams);
  LALInferenceRemoveVariable(&intrinsicParams, "rightascension");
  LALInferenceRemoveVariable(&intrinsicParams, "declination");
  LALInferenceRemoveVariable(&intrinsicParams, "polarisation");
  LALInferenceRemoveVariable(&intrinsicParams, "time");
  LALInferenceRemoveVariable(&intrinsicParams, "distance");
  /*  TODO: add pointer to template function here.                                         */
  /*  (otherwise same parameters but different template will lead to no re-computation!!)  */

  chisquared = 0.0;
  /* loop over data (different interferometers): */
  dataPtr = data;

  while (dataPtr != NULL) {
    /* The parameters the Likelihood function can handle by itself    */
    /* (and which shouldn't affect the template function) are         */
    /* sky location (ra, dec), polarisation and signal arrival time.  */
    /* Note that the template function shifts the waveform to so that */
    /* t_c corresponds to the "time" parameter in                     */
    /* IFOdata->modelParams (set, e.g., from the trigger value).      */
    
    /* Reset log-likelihood */
    dataPtr->loglikelihood = 0.0;

    /* Compare parameter values with parameter values corresponding */
    /* to currently stored template; ignore "time" variable:        */
    if (LALInferenceCheckVariable(dataPtr->modelParams, "time")) {
      timeTmp = *(REAL8 *) LALInferenceGetVariable(dataPtr->modelParams, "time");
      LALInferenceRemoveVariable(dataPtr->modelParams, "time");
    }
    else timeTmp = GPSdouble;
    different = LALInferenceCompareVariables(dataPtr->modelParams, &intrinsicParams);
    /* "different" now may also mean that "dataPtr->modelParams" */
    /* wasn't allocated yet (as in the very 1st iteration).      */

    if (different) { /* template needs to be re-computed: */
      LALInferenceCopyVariables(&intrinsicParams, dataPtr->modelParams);
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
      template(dataPtr);
      if (dataPtr->modelDomain == LALINFERENCE_DOMAIN_TIME) {
	if (!timeDomainWarning) {
	  timeDomainWarning = 1;
	  fprintf(stderr, "WARNING: using time domain template with frequency domain likelihood (in %s, line %d)\n", __FILE__, __LINE__);
	}
        LALInferenceExecuteFT(dataPtr);
        /* note that the dataPtr->modelParams "time" element may have changed here!! */
        /* (during "template()" computation)  */
      }
    }
    else { /* no re-computation necessary. Return back "time" value, do nothing else: */
      LALInferenceAddVariable(dataPtr->modelParams, "time", &timeTmp, LALINFERENCE_REAL8_t,LALINFERENCE_PARAM_LINEAR);
    }

    /*-- Template is now in dataPtr->freqModelhPlus and dataPtr->freqModelhCross. --*/
    /*-- (Either freshly computed or inherited.)                                  --*/

    /* determine beam pattern response (F_plus and F_cross) for given Ifo: */
    XLALComputeDetAMResponse(&Fplus, &Fcross,
                             dataPtr->detector->response,
			     ra, dec, psi, gmst);
    /* signal arrival time (relative to geocenter); */
    timedelay = XLALTimeDelayFromEarthCenter(dataPtr->detector->location,
                                             ra, dec, &GPSlal);
    /* (negative timedelay means signal arrives earlier at Ifo than at geocenter, etc.)    */
    /* amount by which to time-shift template (not necessarily same as above "timedelay"): */
    timeshift =  (GPSdouble - (*(REAL8*) LALInferenceGetVariable(dataPtr->modelParams, "time"))) + timedelay;
    twopit    = LAL_TWOPI * timeshift;

    /* include distance (overall amplitude) effect in Fplus/Fcross: */
    FplusScaled  = Fplus  / distMpc;
    FcrossScaled = Fcross / distMpc;

    if (LALInferenceCheckVariable(currentParams, "crazyInjectionHLSign") &&
        *((INT4 *)LALInferenceGetVariable(currentParams, "crazyInjectionHLSign"))) {
      if (strstr(dataPtr->name, "H") || strstr(dataPtr->name, "L")) {
        FplusScaled *= -1.0;
        FcrossScaled *= -1.0;
      }
    }

    dataPtr->fPlus = FplusScaled;
    dataPtr->fCross = FcrossScaled;
    dataPtr->timeshift = timeshift;

    /* extract the element from the "df" vector that carries the current Ifo's name: */
    CHAR df_variable_name[64];
    sprintf(df_variable_name,"df_%s",dataPtr->name);
    if(LALInferenceCheckVariable(currentParams,df_variable_name)){
      degreesOfFreedom = *(REAL8*) LALInferenceGetVariable(currentParams,df_variable_name);
    }
    else {
      degreesOfFreedom = dataPtr->STDOF;
    }
    if (!(degreesOfFreedom>0)) {
      XLALPrintError(" ERROR in StudentTLogLikelihood(): degrees-of-freedom parameter must be positive.\n");
      XLAL_ERROR_REAL8(XLAL_EDOM);
    }

    /* determine frequency range & loop over frequency bins: */
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    lower = (UINT4)ceil(dataPtr->fLow / deltaF);
    upper = (UINT4)floor(dataPtr->fHigh / deltaF);
    FourDeltaToverN = 4.0 * deltaT / ((double) dataPtr->timeData->data->length);
    for (i=lower; i<=upper; ++i){
      /* degrees-of-freedom parameter (nu_j) for this particular frequency bin: */
      nu = degreesOfFreedom;
      /* (for now constant across frequencies)                                  */
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * dataPtr->freqModelhPlus->data->data[i].re  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].re;
      plainTemplateImag = FplusScaled * dataPtr->freqModelhPlus->data->data[i].im  
                          +  FcrossScaled * dataPtr->freqModelhCross->data->data[i].im;

      /* do time-shifting...            */
      /* (also un-do 1/deltaT scaling): */
      f = ((double) i) * deltaF;
      /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
      re = cos(twopit * f);
      im = - sin(twopit * f);
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      dataReal     = dataPtr->freqData->data->data[i].re / deltaT;
      dataImag     = dataPtr->freqData->data->data[i].im / deltaT;
      /* compute squared difference & 'chi-squared': */
      diffRe       = dataReal - templateReal;         /* Difference in real parts...                     */
      diffIm       = dataImag - templateImag;         /* ...and imaginary parts, and...                  */
      diffSquared  = diffRe*diffRe + diffIm*diffIm ;  /* ...squared difference of the 2 complex figures. */
      singleFreqBinTerm = ((nu+2.0)/2.0) * log(1.0 + (FourDeltaToverN * diffSquared) / (nu * dataPtr->oneSidedNoisePowerSpectrum->data->data[i]));
      chisquared  += singleFreqBinTerm;   /* (This is a sum-of-squares, or chi^2, term in the Gaussian case, not so much in the Student-t case...)  */
      dataPtr->loglikelihood -= singleFreqBinTerm;
    }
    dataPtr = dataPtr->next;
  }
  loglikeli = -1.0 * chisquared; /* note (again): the log-likelihood is unnormalised! */
  LALInferenceDestroyVariables(&intrinsicParams);  
  return(loglikeli);
}

