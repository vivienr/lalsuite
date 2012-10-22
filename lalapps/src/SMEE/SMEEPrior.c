/*
 *  LALInferencePrior.c:  Nested Sampling using LALInference
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

#include <SMEEPrior.h>
#include <lal/LALInferencePrior.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Private helper function prototypes */
//static double qInnerIntegrand(double M2, void *viData);
//static double etaInnerIntegrand(double M2, void *viData);
//static double outerIntegrand(double M1, void *voData);


REAL8 LALInferenceSMEEPrior(LALInferenceRunState *runState, LALInferenceVariables *params)
{
  if (params==NULL || runState == NULL || runState->priorArgs == NULL)
     XLAL_ERROR_REAL8(XLAL_EFAULT, "Null arguments received.");
  
  REAL8 logPrior=0.0;
  if( !in_range( runState->priorArgs, params ) ) return -INFINITY;
  LALInferenceVariableItem *item=params->head;
  LALInferenceVariables *priorParams=runState->priorArgs;
  REAL8 min=0., max=0.;
  //REAL8 A_min=0.0, A_max=0.0;
 // REAL8 logmc=0.0;
  //REAL8 B_min=0.0, B_max=0.0;
  
  for(;item;item=item->next)
  {
    if(LALInferenceCheckMinMaxPrior(priorParams, item->name)){
      LALInferenceGetMinMaxPrior(priorParams, item->name, &min, &max);
      
      logPrior-=log(max-min);
    }
  }
  return(logPrior);
}


UINT4 in_range( LALInferenceVariables *priors, LALInferenceVariables *params ){
 LALInferenceVariableItem *item = params->head;
 REAL8 min, max;

 /* loop over variables */
 for(; item; item = item->next ){
   if( item->vary == LALINFERENCE_PARAM_FIXED ||
     item->vary == LALINFERENCE_PARAM_OUTPUT ){ continue; }

   if( LALInferenceCheckMinMaxPrior( priors, item->name ) ){
     LALInferenceGetMinMaxPrior( priors, item->name, &min, &max );

     /* For cyclic boundaries, mod out by range. */
     if( item->vary == LALINFERENCE_PARAM_CIRCULAR ) {
       REAL8 val = *(REAL8 *)item->value;
       REAL8 delta = max - min;

       if (val > max) {
         REAL8 offset = val - min;

         *(REAL8 *)item->value = min + fmod(offset, delta);
       }
       else {
         REAL8 offset = max - val;

         *(REAL8 *)item->value = max - fmod(offset, delta);
       }

       continue;
     }

     if( (*(REAL8 *) item->value) < min || (*(REAL8 *)item->value) > max )
       return 0;
   }
 }

 return 1;
}

#undef MIN*/



