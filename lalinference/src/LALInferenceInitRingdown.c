#include <stdio.h>
#include <lal/Date.h>
#include <lal/GenerateInspiral.h>
#include <lal/LALInference.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/StringInput.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/TimeSeries.h>
#include <lal/LALInferencePrior.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceProposal.h>
#include <lal/LALInferenceLikelihood.h>
#include <lal/LALInferenceReadData.h>
#include <lal/LALInferenceInit.h>
#include <lal/LALInferenceCalibrationErrors.h>

static void LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model);
static int checkParamInList(const char *list, const char *param);

static int checkParamInList(const char *list, const char *param)
{
  /* Check for param in comma-seperated list */
  char *post=NULL,*pos=NULL;
  if (list==NULL) return 0;
  if (param==NULL) return 0;
  
  if(!(pos=strstr(list,param))) return 0;
  
  /* The string is a substring. Check that it is a token */
  /* Check the character before and after */
  if(pos!=list)
  if(*(pos-1)!=',')
  return 0;
  
  post=&(pos[strlen(param)]);
  if(*post!='\0')
  if(*post!=',')
  return 0;
  return 1;
}

/*
 * Initialize threads in memory, using LALInferenceInitRingdownModel() to init models.
 */
void LALInferenceInitRingdownThreads(LALInferenceRunState *run_state, INT4 nthreads) {
  if (run_state == NULL){
    LALInferenceInitRingdownModel(run_state);
    return;
  }

  ProcessParamsTable *commandLine=run_state->commandLine;

  LALInferenceThreadState *thread;
  INT4 t, nifo;
  INT4 randomseed;
  LALInferenceIFOData *data = run_state->data;
  run_state->nthreads = nthreads;
  run_state->threads = LALInferenceInitThreads(nthreads);

  for (t = 0; t < nthreads; t++) {
    thread = run_state->threads[t];

    /* Link back to run-state */
    thread->parent = run_state;

    /* Set up Ringdown model and parameter array */
    thread->model = LALInferenceInitRingdownModel(run_state);
    thread->model->roq_flag = 0;

    /* Allocate IFO likelihood holders */
    nifo = 0;
    while (data != NULL) {
        data = data->next;
        nifo++;
    }
    thread->currentIFOLikelihoods = XLALCalloc(nifo, sizeof(REAL8));

    /* Setup ROQ */
    if (LALInferenceGetProcParamVal(commandLine, "--roqtime_steps")){

        LALInferenceSetupROQmodel(thread->model, commandLine);
        fprintf(stderr, "done LALInferenceSetupROQmodel\n");

    }else{
      thread->model->roq_flag=0;
    }

    LALInferenceCopyVariables(thread->model->params, thread->currentParams);
    LALInferenceCopyVariables(run_state->proposalArgs, thread->proposalArgs);

    /* Link thread-state prior-args to the parent run-state's */
    thread->priorArgs = run_state->priorArgs;

    /* Use clocktime if seed isn't provided */
    thread->GSLrandom = gsl_rng_alloc(gsl_rng_mt19937);
    randomseed = gsl_rng_get(run_state->GSLrandom);
    gsl_rng_set(thread->GSLrandom, randomseed);
  }

  return;
}


/* Setup the variables to control template generation for the RingdownFD model */
/* Includes specification of prior ranges. Returns address of new LALInferenceVariables */

LALInferenceModel *LALInferenceInitRingdownModel(LALInferenceRunState *state) {
  char help[]="Good luck on your own ... harhar.\n";

  /* Print command line arguments if state was not allocated */
  if(state==NULL)
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  /* Print command line arguments if help requested */
  if(LALInferenceGetProcParamVal(state->commandLine,"--help"))
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  LALStatus status;
  memset(&status,0,sizeof(status));
  //int errnum;
  SimInspiralTable *injTable=NULL;
  LALInferenceVariables *priorArgs=state->priorArgs;
  LALInferenceVariables *proposalArgs=state->proposalArgs;
  ProcessParamsTable *commandLine=state->commandLine;
  ProcessParamsTable *ppt=NULL;
  //ProcessParamsTable *ppt_order=NULL;
  LALPNOrder PhaseOrder=-1;
  LALPNOrder AmpOrder=-1;
  Approximant approx=NumApproximants;
  REAL8 f_ref = 100.0;
  LALInferenceIFOData *dataPtr;
  UINT4 event=0;
  UINT4 i=0;

  /* Default priors */
  REAL8 Dmin=20.0;
  REAL8 Dmax=80.0;
  REAL8 Dinitial = (Dmax + Dmin)/2.0;
  REAL8 BH_massMin = 5.0;
  REAL8 BH_massMax = 70.0;
  REAL8 chieffMin = 0.0;
  REAL8 chieffMax = 1.0;
  REAL8 BH_spinMin = 0.0;
  REAL8 BH_spinMax = 1.0;  
  REAL8 etaMin=0.0312;
  REAL8 etaMax=0.25;
  REAL8 psiMin=0.0,psiMax=LAL_PI;
  REAL8 decMin=-LAL_PI/2.0,decMax=LAL_PI/2.0;
  REAL8 raMin=0.0,raMax=LAL_TWOPI;
  REAL8 phiMin=0.0,phiMax=LAL_TWOPI;
  REAL8 costhetaJNmin=-1.0 , costhetaJNmax=1.0;
  REAL8 dt=0.1;  /* Half the width of time prior */
  gsl_rng *GSLrandom=state->GSLrandom;
  REAL8 endtime=0.0, timeParam=0.0;
  REAL8 timeMin=endtime-dt,timeMax=endtime+dt;
  REAL8 zero=0.0; /* just a number that will be overwritten anyway*/

  /* Over-ride prior bounds if analytic test */
  if (LALInferenceGetProcParamVal(commandLine, "--correlatedGaussianLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence(state));
  }
  else if (LALInferenceGetProcParamVal(commandLine, "--bimodalGaussianLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence_bimod(state));
  }
  else if (LALInferenceGetProcParamVal(commandLine, "--rosenbrockLikelihood"))
  {
    return(LALInferenceInitModelReviewEvidence(state)); /* CHECKME: Use the default prior for unimodal */
  }

  LALInferenceModel *model = XLALMalloc(sizeof(LALInferenceModel));
  model->params = XLALCalloc(1, sizeof(LALInferenceVariables));
  memset(model->params, 0, sizeof(LALInferenceVariables));

  UINT4 signal_flag=1;
  ppt = LALInferenceGetProcParamVal(commandLine, "--noiseonly");
  if(ppt)signal_flag=0;
  LALInferenceAddVariable(model->params, "signalModelFlag", &signal_flag,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);

  /* Read injection XML file for parameters if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--inj");
  if(ppt){
    SimInspiralTableFromLIGOLw(&injTable,ppt->value,0,0);
    if(!injTable){
      fprintf(stderr,"Unable to open injection file %s\n",ppt->value);
      exit(1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--event");
    if(ppt){
      event= atoi(ppt->value);
      fprintf(stderr,"Reading event %d from file\n",event);
      i=0;
      while(i<event) {i++; injTable=injTable->next;} /* select event */
    }
  }

  /* See if there are any parameters pinned to injection values */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--pinparams"))){
    char *pinned_params=ppt->value;
    LALInferenceVariables tempParams;
    memset(&tempParams,0,sizeof(tempParams));
    char **strings=NULL;
    UINT4 N;
    LALInferenceParseCharacterOptionString(pinned_params,&strings,&N);
    LALInferenceInjectionToVariables(injTable,&tempParams);
    LALInferenceVariableItem *node=NULL;
    while(N>0){
      N--;
      char *name=strings[N];
      fprintf(stdout,"Pinning parameter %s\n",node->name);
      node=LALInferenceGetItem(&tempParams,name);
      if(node) LALInferenceAddVariable(model->params,node->name,node->value,node->type,node->vary);
      else {fprintf(stderr,"Error: Cannot pin parameter %s. No such parameter found in injection!\n",node->name);}
    }
  }

  approx=RingdownFD;

  // Set the model domain appropriately 
  if (XLALSimInspiralImplementedFDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  } else if (XLALSimInspiralImplementedTDApproximants(approx)) {
    model->domain = LAL_SIM_DOMAIN_TIME;
  } else {
    fprintf(stderr,"ERROR. Unknown approximant number %i. Unable to choose time or frequency domain model.",approx);
    exit(1);
  }

  ppt=LALInferenceGetProcParamVal(commandLine, "--fref");
  if (ppt) f_ref = atof(ppt->value);

  ppt=LALInferenceGetProcParamVal(commandLine,"--modeldomain");
  if(ppt){
    if ( ! strcmp( "time", ppt->value ) )
    {
      model->domain = LAL_SIM_DOMAIN_TIME;
    }
    else if ( ! strcmp( "frequency", ppt->value ) )
    {
      model->domain = LAL_SIM_DOMAIN_FREQUENCY;
    }
    else
    {
      fprintf( stderr, "invalid argument to --modeldomain:\n"
              "unknown domain specified: "
              "domain must be one of: time, frequency\n");
      exit( 1 );
    }
  }

  /* This sets the component masses and total mass priors, if given in command line.
   * The prior for other parameters are now read in in RegisterUniformVariable, if given by the user. */
  //LALInferenceInitMassVariables(state);
  /* now we need to update the chirp mass and q limits accordingly */
  /*REAL8 m2_min=LALInferenceGetREAL8Variable(state->priorArgs,"mass2_min");
  REAL8 m1_max=LALInferenceGetREAL8Variable(state->priorArgs,"mass1_max");
  REAL8 mtot_min = *(REAL8 *)LALInferenceGetVariable(state->priorArgs,"MTotMin");
  REAL8 mtot_max = *(REAL8 *)LALInferenceGetVariable(state->priorArgs,"MTotMax");
  qMin = m2_min / m1_max;
  mcMin =mtot_min*pow(qMin/pow(1.+qMin,2.),3./5.);
  mcMax =mtot_max*pow(0.25,3./5.);*/

  /************ Initial Value Related Argument START *************/
  /* Read time parameter from injection file */
  if(injTable)
  {
    endtime=XLALGPSGetREAL8(&(injTable->geocent_end_time));
    fprintf(stdout,"Using end time from injection file: %lf\n", endtime);
  }
  /* Over-ride end time if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--trigtime");
  if(ppt){
    endtime=atof(ppt->value);
  }
  /* Over-ride time prior window if specified */
  ppt=LALInferenceGetProcParamVal(commandLine,"--dt");
  if(ppt)
    dt=atof(ppt->value);
  timeMin=endtime-dt; timeMax=endtime+dt;
  timeParam = timeMin + (timeMax-timeMin)*gsl_rng_uniform(GSLrandom);

  /* Initial Value Related END */
  printf("Adding LAL_APPROXIMANT\n");
  LALInferenceAddVariable(model->params, "LAL_APPROXIMANT", &approx,        LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "LAL_PNORDER",     &PhaseOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "LAL_AMPORDER",     &AmpOrder,        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  LALInferenceAddVariable(model->params, "f_ref", &f_ref, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

  /* flow handling */
  REAL8 fLow = state->data->fLow;
  ppt=LALInferenceGetProcParamVal(commandLine,"--vary-flow");
  if(ppt){
    REAL8 fLow_min = fLow;
    REAL8 fLow_max = 200.0;
    if(LALInferenceCheckVariable(model->params,"f_ref"))
      f_ref = *(REAL8*)LALInferenceGetVariable(model->params, "f_ref");
      if (f_ref > 0.0 && fLow_max > f_ref) {
        fprintf(stdout,"WARNING: flow can't go higher than the reference frequency.  Setting flow-max to %f\n",f_ref);
        fLow_max = f_ref;
      }
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "flow", fLow, fLow_min, fLow_max, LALINFERENCE_PARAM_LINEAR);
  } else {
    LALInferenceAddVariable(model->params, "flow", &fLow,  LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);
  }

  /* Set up the variable parameters */

  /********************* TBL: Adding noise-fitting parameters  *********************/
  UINT4 nscale_block; //number of noise parameters per IFO (1 per frequency block)
  UINT4 nscale_bin;   //number of Fourier bins in each noise block
  REAL8 nscale_dflog; //logarithmic spacing for noise parameters
  REAL8 nscale_min;   //minimum value for psd scale parameter
  REAL8 nscale_max;   //maximum value for psd scale parameters
  UINT4 nscale_dim;   //total dimension of noise model (params X detectors)
  UINT4 nscale_flag;  //flag to tell likelihood function if psd fitting is in use

  REAL8Vector *nscale_prior = NULL; //std. dev. of prior distribution
  REAL8Vector *nscale_sigma = NULL; //std. dev. of prior distribution

  //assume no noise fitting
  nscale_flag=0;

  //set Nblock to default unless specified at command line
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdNblock");
  if(ppt) nscale_block = atoi(ppt->value);
  else nscale_block = 8;

  //First, figure out sizes of dataset to set up noise blocks
  UINT4 nifo; //number of data channels
  UINT4 imin; //minimum Fourier bin for integration in IFO
  UINT4 imax; //maximum Fourier bin for integration in IFO
  UINT4 f_min = 1; //minimum Fourier bin for integration over network
  UINT4 f_max = 1; //maximum Fourier bin for integration over network
  REAL8 df = 1.0; //frequency resolution

  /* Set model sampling rates to be consistent with data */
  model->deltaT = state->data->timeData->deltaT;
  model->deltaF = state->data->freqData->deltaF;

  /* Get number of interferometers */
  nifo=0;
  dataPtr = state->data;
  while (dataPtr != NULL)
  {
    df      = 1.0 / (((double)dataPtr->timeData->data->length) * model->deltaT);
    imin    = (UINT4)ceil( dataPtr->fLow  / df);
    imax    = (UINT4)floor(dataPtr->fHigh / df);

    if(nifo==0)
    {
      f_min=imin;
      f_max=imax;
    }
    else
    {
      if(imin<f_min)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_min=imin;
      }
      if(imax>f_max)
      {
        fprintf(stderr,"Warning: Different IFO's have different minimum frequencies -- bad for noise fitting\n");
        f_max=imax;
      }
    }

    dataPtr = dataPtr->next;
    nifo++;
  }
  imin = f_min;
  imax = f_max;

  UINT4 j = 0;

  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFit");
  if (ppt) {
      printf("WARNING: --psdFit has been deprecated in favor of --psd-fit\n");
  } else {
      ppt = LALInferenceGetProcParamVal(commandLine, "--psd-fit");
  }
  if(ppt)//MARK: Here is where noise PSD parameters are being added to the model
  {

    printf("Setting up PSD fitting for %i ifos...\n",nifo);

    dataPtr = state->data;

    gsl_matrix *bands_min = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix *bands_max = gsl_matrix_alloc(nifo,nscale_block);

    i=0;
    while (dataPtr != NULL)
    {
      printf("ifo=%i  %s\n",i,dataPtr->name);fflush(stdout);

        nscale_bin   = (imax+1-imin)/nscale_block;
        nscale_dflog = log( (double)(imax+1)/(double)imin )/(double)nscale_block;

        int freq_min, freq_max;

        for (j = 0; j < nscale_block; j++)
        {

            freq_min = (int) exp(log((double)imin ) + nscale_dflog*j);
            freq_max = (int) exp(log((double)imin ) + nscale_dflog*(j+1));

            gsl_matrix_set(bands_min,i,j,freq_min);
            gsl_matrix_set(bands_max,i,j,freq_max);
        }


      dataPtr = dataPtr->next;
      i++;

    }

    printf("Running PSD fitting with bands (Hz)...\n");
    dataPtr = state->data;
    i=0;
    while (dataPtr != NULL)
    {
      printf("%s:",dataPtr->name);
      for (j = 0; j < nscale_block; j++)
      {
        printf(" %f-%f ",gsl_matrix_get(bands_min,i,j)*df,gsl_matrix_get(bands_max,i,j)*df);
      }
      printf("\n");
      dataPtr = dataPtr->next;
      i++;
    }

    nscale_bin   = (f_max+1-f_min)/nscale_block;
    nscale_dflog = log( (double)(f_max+1)/(double)f_min )/(double)nscale_block;

    nscale_min   = 1.0e-1;
    nscale_max   = 1.0e+1;
    nscale_dim   = nscale_block*nifo;
    nscale_flag  = 1;

    // Set noise parameter arrays.
    nscale_prior = XLALCreateREAL8Vector(nscale_block);
    nscale_sigma = XLALCreateREAL8Vector(nscale_block);
    for(i=0; i<nscale_block; i++)
    {
      nscale_prior->data[i] = 1.0/sqrt( gsl_matrix_get(bands_max,0,i)-gsl_matrix_get(bands_min,0,i) );
      nscale_sigma->data[i] = nscale_prior->data[i]/sqrt((double)(nifo*nscale_block));
    }

    gsl_matrix *nscale = gsl_matrix_alloc(nifo,nscale_block);
    gsl_matrix_set_all(nscale, 1.0);

    LALInferenceAddVariable(model->params, "psdscale", &nscale, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_LINEAR);
    LALInferenceAddVariable(model->params, "logdeltaf", &nscale_dflog, LALINFERENCE_REAL8_t, LALINFERENCE_PARAM_FIXED);

    LALInferenceAddVariable(model->params, "psdBandsMin", &bands_min, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(model->params, "psdBandsMax", &bands_max, LALINFERENCE_gslMatrix_t, LALINFERENCE_PARAM_FIXED);

    //Set up noise priors
    LALInferenceAddVariable(priorArgs,      "psddim",   &nscale_dim,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddMinMaxPrior(priorArgs,   "psdscale", &nscale_min,  &nscale_max,   LALINFERENCE_REAL8_t);
    LALInferenceAddMinMaxPrior(priorArgs,   "psdrange", &nscale_min,  &nscale_max,   LALINFERENCE_REAL8_t);
    LALInferenceAddVariable(priorArgs,      "psdsigma", &nscale_prior, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);

    //Store meta data for noise model in proposal
    LALInferenceAddVariable(proposalArgs, "psdblock", &nscale_block, LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(proposalArgs, "psdbin",   &nscale_bin,   LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);
    LALInferenceAddVariable(proposalArgs, "psdsigma", &nscale_sigma, LALINFERENCE_REAL8Vector_t, LALINFERENCE_PARAM_FIXED);


  }//End of noise model initialization
  LALInferenceAddVariable(model->params, "psdScaleFlag", &nscale_flag, LALINFERENCE_UINT4_t, LALINFERENCE_PARAM_FIXED);

  UINT4 psdGaussianPrior=1;
  ppt = LALInferenceGetProcParamVal(commandLine, "--psdFlatPrior");
  if(ppt)psdGaussianPrior=0;
  LALInferenceAddVariable(priorArgs, "psdGaussianPrior", &psdGaussianPrior,  LALINFERENCE_INT4_t,  LALINFERENCE_PARAM_FIXED);

  ppt = LALInferenceGetProcParamVal(commandLine, "--glitchFit");
  if (ppt)
      printf("WARNING: --glitchFit has been deprecated in favor of --glitch-fit\n");
  else
      ppt = LALInferenceGetProcParamVal(commandLine, "--glitch-fit");
  if (ppt)
      LALInferenceInitGlitchVariables(state, model->params);

  /* Handle, if present, requests for calibration parameters. */
  LALInferenceInitCalibrationVariables(state, model->params);

  //Only add waveform parameters to model if needed
  if(signal_flag)
  {  
    /* The idea here is the following:
     * We call RegisterUniformVariable with startval=0 and meanigful min and max values.
     * That function will then take care of setting startval to a random value between min and max, or read a value from command line (with --parname VALUE).
     * The user can fix the param to a given value with --fix-parname --parname VALUE
     * */
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--BH_mass-max"))) BH_massMax=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--BH_mass-min"))) BH_massMin=atof(ppt->value);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "BH_mass", zero, BH_massMin, BH_massMax, LALINFERENCE_PARAM_LINEAR);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--BH_spin-max"))) BH_spinMax=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--BH_spin-min"))) BH_spinMin=atof(ppt->value);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "BH_spin", zero, BH_spinMin, BH_spinMax, LALINFERENCE_PARAM_LINEAR);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--eta-max"))) etaMax=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--eta-min"))) etaMin=atof(ppt->value);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "eta", zero, etaMin, etaMax, LALINFERENCE_PARAM_LINEAR);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--chieff-max"))) chieffMax=atof(ppt->value);
    if((ppt=LALInferenceGetProcParamVal(commandLine,"--chieff-min"))) chieffMin=atof(ppt->value);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "chieff", zero, chieffMin, chieffMax, LALINFERENCE_PARAM_LINEAR);
    if(!LALInferenceGetProcParamVal(commandLine,"--margphi") && !LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
      LALInferenceRegisterUniformVariableREAL8(state, model->params, "phi0", zero, phiMin, phiMax, LALINFERENCE_PARAM_CIRCULAR);
    }

  /* Check for distance prior for use if the user samples in logdistance */
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--distance-max"))) Dmax=atof(ppt->value);
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--distance-min"))) Dmin=atof(ppt->value);
  LALInferenceParamVaryType distanceVary = LALINFERENCE_PARAM_LINEAR;
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--fix-distance")))
  {
    Dinitial=atof(ppt->value);
    distanceVary = LALINFERENCE_PARAM_FIXED;
  }

  LALInferenceRegisterUniformVariableREAL8(state, model->params, "logdistance", log(Dinitial), log(Dmin), log(Dmax), distanceVary);
  LALInferenceRegisterUniformVariableREAL8(state, model->params, "polarisation", zero, psiMin, psiMax, LALINFERENCE_PARAM_LINEAR);
  LALInferenceRegisterUniformVariableREAL8(state, model->params, "costheta_jn", zero, costhetaJNmin, costhetaJNmax,LALINFERENCE_PARAM_LINEAR);

  /* Option to use the detector-aligned frame */
  if(!LALInferenceGetProcParamVal(commandLine,"--no-detector-frame") && nifo >1)
  {
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"t0",timeParam,timeMin,timeMax,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"cosalpha",0,-1,1,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state,model->params,"azimuth",0.0,0.0,LAL_TWOPI,LALINFERENCE_PARAM_CIRCULAR);
        /* add the time parameter then remove it so that the prior is set up properly */
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", timeParam, timeMin, timeMax,LALINFERENCE_PARAM_LINEAR);
        LALInferenceRemoveVariable(model->params,"time");
        INT4 one=1;
        LALInferenceAddVariable(model->params,"SKY_FRAME",&one,LALINFERENCE_INT4_t,LALINFERENCE_PARAM_FIXED);
  }
  else
  {
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "rightascension", zero, raMin, raMax,LALINFERENCE_PARAM_CIRCULAR);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "declination", zero, decMin, decMax, LALINFERENCE_PARAM_LINEAR);
        LALInferenceRegisterUniformVariableREAL8(state, model->params, "time", timeParam, timeMin, timeMax,LALINFERENCE_PARAM_LINEAR);
   }

  /* If we are marginalising over the time, remove that variable from the model (having set the prior above) */
  /* Also set the prior in model->params, since Likelihood can't access the state! (ugly hack) */
  if(LALInferenceGetProcParamVal(commandLine,"--margtime") || LALInferenceGetProcParamVal(commandLine, "--margtimephi")){
	  LALInferenceVariableItem *p=LALInferenceGetItem(state->priorArgs,"time_min");
	  LALInferenceAddVariable(model->params,"time_min",p->value,p->type,p->vary);
	  p=LALInferenceGetItem(state->priorArgs,"time_max");
	  LALInferenceAddVariable(model->params,"time_max",p->value,p->type,p->vary);
	  if (LALInferenceCheckVariable(model->params,"time")) LALInferenceRemoveVariable(model->params,"time");
      if (LALInferenceCheckVariable(model->params,"t0")) LALInferenceRemoveVariable(model->params,"t0");
	  if (LALInferenceGetProcParamVal(commandLine, "--margtimephi")) {
		  UINT4 margphi = 1;
		  LALInferenceAddVariable(model->params, "margtimephi", &margphi, LALINFERENCE_UINT4_t,LALINFERENCE_PARAM_FIXED);
	  }
  }

    /* If requested by the user populate the testing GR or PPE model parameters */
  if (LALInferenceGetProcParamVal(commandLine,"--grtest-parameters") || LALInferenceGetProcParamVal(commandLine,"--ppe-parameters"))
  {
    printf("InitNonGRParams\n");
    LALInferenceInitNonGRParams(state, model);
  }
  /* PPE parameters */

  ppt=LALInferenceGetProcParamVal(commandLine, "--TaylorF2ppE");
  if(approx==TaylorF2 && ppt){

    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppealpha",zero, -1000.0 , 1000.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppebeta", zero, -1000.0 , 1000.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppeuppera", zero, -3.0, 3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppeupperb", zero, -3.0, 3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppelowera", zero, -3.0, 2.0/3.0 , LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "ppelowerb", zero, -4.5, 1.0, LALINFERENCE_PARAM_LINEAR);

  }

 /* if(LALInferenceGetProcParamVal(commandLine,"--tidalT")&&LALInferenceGetProcParamVal(commandLine,"--tidal")){
    XLALPrintError("Error: cannot use both --tidalT and --tidal.\n");
    XLAL_ERROR_NULL(XLAL_EINVAL);
  } else if(LALInferenceGetProcParamVal(commandLine,"--tidalT")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambdaT", zero, lambdaTMin, lambdaTMax, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "dLambdaT", zero, dLambdaTMin, dLambdaTMax, LALINFERENCE_PARAM_LINEAR);

  } else if(LALInferenceGetProcParamVal(commandLine,"--tidal")){
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambda1", zero, lambda1Min, lambda1Max, LALINFERENCE_PARAM_LINEAR);
    LALInferenceRegisterUniformVariableREAL8(state, model->params, "lambda2", zero, lambda2Min, lambda2Max, LALINFERENCE_PARAM_LINEAR);

  }*/

  LALSimInspiralSpinOrder spinO = LAL_SIM_INSPIRAL_SPIN_ORDER_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine, "--spinOrder");
  if(ppt) {
    spinO = atoi(ppt->value);
    LALInferenceAddVariable(model->params, "spinO", &spinO,
        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }
  LALSimInspiralTidalOrder tideO = LAL_SIM_INSPIRAL_TIDAL_ORDER_ALL;
  ppt=LALInferenceGetProcParamVal(commandLine, "--tidalOrder");
  if(ppt) {
    tideO = atoi(ppt->value);
    LALInferenceAddVariable(model->params, "tideO", &tideO,
        LALINFERENCE_INT4_t, LALINFERENCE_PARAM_FIXED);
  }

  printf("checkpoint\n");

  /* LALInference uses a proprietary spin parameterization, the conversion from which
   * assumes the LALSimulations default frame */
  LALSimInspiralFrameAxis frameAxis = LAL_SIM_INSPIRAL_FRAME_AXIS_DEFAULT;

  model->waveFlags = XLALSimInspiralCreateWaveformFlags();
  XLALSimInspiralSetSpinOrder(model->waveFlags,  spinO);
  XLALSimInspiralSetTidalOrder(model->waveFlags, tideO);
  XLALSimInspiralSetFrameAxis(model->waveFlags,frameAxis);
  if((ppt=LALInferenceGetProcParamVal(commandLine,"--numreldata"))) {
    XLALSimInspiralSetNumrelData(model->waveFlags, ppt->value);
    fprintf(stdout,"Template will use %s.\n",ppt->value);
  }



  fprintf(stdout,"\n\n---\t\t ---\n");
  //LALInferenceInitSpinVariables(state, model);
  //LALInferenceCheckApproximantNeeds(state,approx);

  if (injTable)
     //print_flags_orders_warning(injTable,commandLine);

     /* Print info about orders and waveflags used for templates */

     fprintf(stdout,"Templates will run using Approximant RingdownFD\n");
     fprintf(stdout,"---\t\t ---\n\n");
  }//end of signal only flag
  else
  {
    /* Print info about orders and waveflags used for templates */
    fprintf(stdout,"\n\n------\n");
    fprintf(stdout,"Noise only run\n");
    fprintf(stdout,"------\n\n");
  }

  /* Initialize waveform buffers */
  model->timehPlus  = XLALCreateREAL8TimeSeries("timehPlus",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->timehCross = XLALCreateREAL8TimeSeries("timehCross",
                                                &(state->data->timeData->epoch),
                                                0.0,
                                                model->deltaT,
                                                &lalDimensionlessUnit,
                                                state->data->timeData->data->length);
  model->freqhPlus = XLALCreateCOMPLEX16FrequencySeries("freqhPlus",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);
  model->freqhCross = XLALCreateCOMPLEX16FrequencySeries("freqhCross",
                                                &(state->data->freqData->epoch),
                                                0.0,
                                                model->deltaF,
                                                &lalDimensionlessUnit,
                                                state->data->freqData->data->length);

  model->freqhs = XLALCalloc(nifo, sizeof(COMPLEX16FrequencySeries *));
  for (i=0; i<nifo; i++)
      model->freqhs[i] = XLALCreateCOMPLEX16FrequencySeries("freqh",
                                                            &(state->data->freqData->epoch),
                                                            0.0,
                                                            model->deltaF,
                                                            &lalDimensionlessUnit,
                                                            state->data->freqData->data->length);

  /* Create arrays for holding single-IFO likelihoods, etc. */
  model->ifo_loglikelihoods = XLALCalloc(nifo, sizeof(REAL8));
  model->ifo_SNRs = XLALCalloc(nifo, sizeof(REAL8));

  /* Choose proper template */
  model->templt = LALInferenceInitRingdownTemplate(state);

  /* Use same window and FFT plans on model as data */
  model->window = state->data->window;
  model->timeToFreqFFTPlan = state->data->timeToFreqFFTPlan;
  model->freqToTimeFFTPlan = state->data->freqToTimeFFTPlan;

  /* Initialize waveform cache */
  model->waveformCache = XLALCreateSimInspiralWaveformCache();

  return(model);
}

/* Setup the template generation */
/* Defaults to using LALSimulation */
LALInferenceTemplateFunction LALInferenceInitRingdownTemplate(LALInferenceRunState *runState)
{
  char help[]="(--template LALInferenceTemplateRingdownFD, some stuff to do?\n";

  /* Print command line arguments if state was not allocated */
  if(runState==NULL)
    {
      fprintf(stdout,"%s",help);
      return(NULL);
    }

  ProcessParamsTable *ppt=NULL;
  ProcessParamsTable *commandLine=runState->commandLine;
  /* Print command line arguments if help requested */
  //Help is taken care of in LALInferenceInitCBCVariables
  /*ppt=LALInferenceGetProcParamVal(commandLine,"--help");
  if(ppt)
  {
      fprintf(stdout,"%s",help);
      return ;
  }*/
  /* This is the LAL template generator for inspiral signals */
  LALInferenceTemplateFunction templt = &LALInferenceTemplateRingdownFD;
  ppt=LALInferenceGetProcParamVal(commandLine,"--template");
  if(ppt) {
    if(!strcmp("LALSim",ppt->value))
      templt=&LALInferenceTemplateXLALSimInspiralChooseWaveform;
        else if(!strcmp("null",ppt->value))
        templt=&LALInferenceTemplateNullFreqdomain;
        else if(!strcmp("multiband",ppt->value)){
        templt=&LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated;
        fprintf(stdout,"Template function called is \"LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated\"\n");
    }
    else {
        XLALPrintError("Error: unknown template %s\n",ppt->value);
        XLALPrintError("%s", help);
        XLAL_ERROR_NULL(XLAL_EINVAL);
    }
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--LALSimulation")){
    fprintf(stderr,"Warning: --LALSimulation is deprecated, the LALSimulation package is now the default. To use LALInspiral specify:\n\
                    --template LALGenerateInspiral (for time-domain templates)\n\
                    --template LAL (for frequency-domain templates)\n");
  }
  else if(LALInferenceGetProcParamVal(commandLine,"--roqtime_steps")){
  templt=&LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence;
        fprintf(stderr, "template is \"LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence\"\n");
  }
  else {
    fprintf(stdout,"Template function called is \"LALInferenceTemplateRingdownFD\"\n");
  }
  return templt;
}

/*******************************************************************
 * LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model)
 * Function to initialise either the TaylorF2Test of SpinTaylorT4Test waveform models
 * or the PPE waveform model
 *******************************************************************/
static void LALInferenceInitNonGRParams(LALInferenceRunState *state, LALInferenceModel *model)
{    
    ProcessParamsTable *commandLine = state->commandLine;
    ProcessParamsTable *ppt=NULL;
    /* check that the user does not request both a TaylorF2Test and a PPE waveform model */
    if (LALInferenceGetProcParamVal(commandLine,"--grtest-parameters") && LALInferenceGetProcParamVal(commandLine,"--ppe-parameters"))
    {
        fprintf(stderr,"--grtest-parameters and --ppe-parameters are not simultaneously supported. Please choose one. Aborting\n");
        exit(-1);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--grtest-parameters");
    if (ppt)
    {
        REAL8 dchi_max=1.;
        REAL8 dchi_min=-1.;
        REAL8 dxi_max=1.;
        REAL8 dxi_min=-1.;
        REAL8 dalpha_max=1.;
        REAL8 dalpha_min=-1.;
        REAL8 dbeta_max=1.;
        REAL8 dbeta_min=-1.;
        REAL8 dsigma_max=1.;
        REAL8 dsigma_min=-1.;
        REAL8 domega22_min = -1.;
        REAL8 domega33_min = -1.;
        REAL8 domega22_max = 1.;
        REAL8 domega33_max = 1.;
        REAL8 dtau22_min = -1.;
        REAL8 dtau33_min = -1.;
        REAL8 dtau22_max = 1.;
        REAL8 dtau33_max = 1.;
        REAL8 tmpVal=0.0;
	/* Relative shifts for inspiral phase PN coefficients (absolute value for dchi1) */
        if (checkParamInList(ppt->value,"dchi0")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi0", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi1", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi2", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi3", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi4", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi5", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi5l")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi5l", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi6")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi6", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi6l")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi6l", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dchi7")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dchi7", tmpVal, dchi_min, dchi_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for pre-merger phase coefficients (PhenomC/P) */
        if (checkParamInList(ppt->value,"dxi1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi1", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi2", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi3", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi4", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi5", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dxi6")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dxi6", tmpVal, dxi_min, dxi_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for merger-ringdown phase coefficients  (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dalpha1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha1", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha2", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha3", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha4", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dalpha5")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dalpha5", tmpVal, dalpha_min, dalpha_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for phenomenological inspiral phase coefficients (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dsigma1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma1", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma2", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma3", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dsigma4")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dsigma4", tmpVal, dsigma_min, dsigma_max, LALINFERENCE_PARAM_LINEAR);
	/* Relative shifts for intermediate phase coefficients (PhenomD/Pv2) */
        if (checkParamInList(ppt->value,"dbeta1")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta1", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dbeta2")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta2", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dbeta3")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dbeta3", tmpVal, dbeta_min, dbeta_max, LALINFERENCE_PARAM_LINEAR);
        /* Relative shifts for ringdown */
        if (checkParamInList(ppt->value,"domega22")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "domega22", tmpVal, domega22_min, domega22_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"domega33")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "domega33", tmpVal, domega33_min, domega33_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dtau22")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dtau22", tmpVal, dtau22_min, dtau22_max, LALINFERENCE_PARAM_LINEAR);
        if (checkParamInList(ppt->value,"dtau33")) LALInferenceRegisterUniformVariableREAL8(state, model->params, "dtau33", tmpVal, dtau33_min, dtau33_max, LALINFERENCE_PARAM_LINEAR);
    }
    ppt=LALInferenceGetProcParamVal(commandLine,"--ppe-parameters");
    if (ppt)
    {
        /* amplitude parameters */
        REAL8 appe_min = -5.0,appe_max=5.0;
        REAL8 alphappe_min = -1000.0,alphappe_max=1000.0;
        REAL8 bppe_min = -5.0,bppe_max=5.0;
        REAL8 betappe_min = -1000.0,betappe_max=1000.0;
        char aPPEparam[64]="";
        char alphaPPEparam[64]="";
        /* phase parameters */
        char bPPEparam[64]="";
        char betaPPEparam[64]="";
        int counters[4]={0};
        do
        {
            sprintf(aPPEparam, "%s%d","aPPE",++counters[0]);
            if (checkParamInList(ppt->value,aPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, aPPEparam, 0.0, appe_min, appe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(alphaPPEparam, "%s%d","alphaPPE",++counters[1]);
            if (checkParamInList(ppt->value,alphaPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, alphaPPEparam, 0.0, alphappe_min, alphappe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(bPPEparam, "%s%d","bPPE",++counters[2]);
            if (checkParamInList(ppt->value,bPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, bPPEparam, 0.0, bppe_min, bppe_max, LALINFERENCE_PARAM_LINEAR);
            sprintf(betaPPEparam, "%s%d","betaPPE",++counters[3]);
            if (checkParamInList(ppt->value,betaPPEparam)) LALInferenceRegisterUniformVariableREAL8(state, model->params, betaPPEparam, 0.0, betappe_min, betappe_max, LALINFERENCE_PARAM_LINEAR);
            
        } while((checkParamInList(ppt->value,aPPEparam))||(checkParamInList(ppt->value,alphaPPEparam))||(checkParamInList(ppt->value,bPPEparam))||(checkParamInList(ppt->value,betaPPEparam)));
        if ((counters[0]!=counters[1])||(counters[2]!=counters[3])) {fprintf(stderr,"Unequal number of PPE parameters detected! Check your command line!\n"); exit(-1);}
    }
    
}


