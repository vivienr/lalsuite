/* Copyright (C) 2012 Walter Del Pozzo, Evan Ochsner and Salvatore Vitale
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
 
#include  <lal/LALSimInspiralTestGRParams.h>

/**
 * @addtogroup LALSimInspiralTestGRParams_c
 * @brief Routines to manipulate non-GR parameter key-value pairs.
 * @{
 */

/**
 * Function that creates the head node of the test GR parameters linked list.
 * It is initialized with a single parameter with given name and value
 */
LALSimInspiralTestGRParam *XLALSimInspiralCreateTestGRParam(
        const char *name, /**< Name of first parameter in new linked list */
        double value 	 /**< Value of first parameter in new linked list */
        )
{
        LALSimInspiralTestGRParam *parameter = (LALSimInspiralTestGRParam *)XLALMalloc(sizeof(LALSimInspiralTestGRParam));
        if (parameter) 
        {
            parameter->data =  (LALSimInspiralTestGRParamData *)XLALMalloc(sizeof(LALSimInspiralTestGRParamData));
            memcpy(parameter->data->name, name, 32);
            parameter->data->value = value;
        }
        parameter->next=NULL;
        return parameter;
}

/**
 * Function that adds a prameter to the test GR parameters linked list. If the
 * parameter already exists, it throws an error.
 */
int XLALSimInspiralAddTestGRParam(
        LALSimInspiralTestGRParam **parameter, /**< Pointer to the head node of the linked list of parameters */
        const char *name, 		/**< Parameter name */
        double value 			/**< Parameter value */
        )
{
    LALSimInspiralTestGRParam *temp;
    temp = *parameter;
    if (*parameter==NULL) 
    {
        temp = XLALSimInspiralCreateTestGRParam(name,value); 
        //temp->next=NULL;
        *parameter=temp;
    }
    else 
    {

        if (!XLALSimInspiralTestGRParamExists(*parameter, name))
        {
            temp = *parameter;
             while(temp->next!=NULL) {temp=temp->next;}
            LALSimInspiralTestGRParam *newParam = XLALSimInspiralCreateTestGRParam(name,value);        
            temp->next = newParam;
        }
        else 
        {
            XLALPrintError("XLAL Error - %s: parameter '%s' exists already! Not added to the structure\n",
                    __func__, name);
            XLAL_ERROR(XLAL_EINVAL);
        }
    }
    return XLAL_SUCCESS;
}

/**
 * Function that sets the value of the desired parameter in the test GR
 * parameters linked list to 'value'.  Throws an error if the parameter
 * is not found
 */
int XLALSimInspiralSetTestGRParam(
        LALSimInspiralTestGRParam *parameter, /**< Linked list to be modified */
        const char *name, 		/**< Name of parameter to be modified */
        const double value 		/**< New value for parameter */
        )
{
    if (XLALSimInspiralTestGRParamExists(parameter, name)) 
    {
        while(parameter)
        {
            if(!strcmp(parameter->data->name, name)) parameter->data->value = value;
            parameter=parameter->next;
        }
        return XLAL_SUCCESS;
    }
    else
    {
        XLALPrintError("XLAL Error - %s: parameter '%s' unknown!\n",
                __func__, name);
        XLAL_ERROR(XLAL_EINVAL);
    }
}

/**
 * Function that returns the value of the desired parameters in the
 * test GR parameters linked list.  Aborts if the parameter is not found
 */
double XLALSimInspiralGetTestGRParam(
        const LALSimInspiralTestGRParam *parameter, /**< Linked list to retrieve from */
        const char *name 	   /**< Name of parameter to be retrieved */
        )
{
    if (XLALSimInspiralTestGRParamExists(parameter, name)) 
        {
            while(parameter) 
            {
                if(!strcmp(parameter->data->name, name)) return parameter->data->value;
                parameter=parameter->next;
            }
        }
    else 
    {
        XLALPrintError("XLAL Error - %s: parameter '%s' unknown!\n",
                __func__, name);
        XLAL_ERROR(XLAL_EINVAL);
    }
    return 0.0; // Should not actually get here!
}

/**
 * Function that checks whether the requested parameter exists within the
 * test GR parameters linked list.  Returns true (1) or false (0) accordingly
 */
bool XLALSimInspiralTestGRParamExists(
        const LALSimInspiralTestGRParam *parameter, 	/**< Linked list to check */
        const char *name 		/**< Parameter name to check for */
        )
{
  if(!parameter) return false;
  while(parameter) {if(!strcmp(parameter->data->name, name)) return true; else parameter=parameter->next;}
  return false;
}

/** Function that prints the whole test GR params linked list */
int XLALSimInspiralPrintTestGRParam(
        FILE *fp, 			/**< FILE pointer to write to */
        LALSimInspiralTestGRParam *parameter 	/**< Linked list to print */
        )
{
    if (parameter!=NULL)
    {
        while(parameter) 
        {
            fprintf(fp,"%s %10.5f\n",parameter->data->name,parameter->data->value);
            parameter=parameter->next;
        }
        return XLAL_SUCCESS;
    }
    else
    {
        XLALPrintError("XLAL Error - %s: parameter not allocated!\n", __func__);
        XLAL_ERROR(XLAL_EINVAL);
    }
}

/** Function that destroys the whole test GR params linked list */
void XLALSimInspiralDestroyTestGRParam(
        LALSimInspiralTestGRParam *parameter 	/**< Linked list to destroy */
        )
{
   LALSimInspiralTestGRParam *tmp;
   while(parameter){
	tmp=parameter->next;
	XLALFree(parameter->data);
	XLALFree(parameter);
	parameter=tmp;
	}
}

int XLALSimInspiralWaveformParamsNonGRAreDefault(LALDict *params)
{
  return (XLALSimInspiralWaveformParamsNonGRPhi1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRPhi2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRPhi3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRPhi4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi0IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi5IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi5LIsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi6IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi6LIsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDChi7IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi5IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDXi6IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDSigma1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDSigma2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDSigma3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDSigma4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDAlpha1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDAlpha2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDAlpha3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDAlpha4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDAlpha5IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDBeta1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDBeta2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRDBeta3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPEIsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE0IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE1IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE2IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE3IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE4IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE5IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE6IsDefault(params)
	  && XLALSimInspiralWaveformParamsNonGRAlphaPPE7IsDefault(params));
}

int XLALSimInspiralWaveformParamsInsertNonGRParamsGeneric(LALDict *LALpars, const char *name, const REAL8 value)
{
  if (strcmp(name,"Phi1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRPhi1(LALpars,value);
  else if (strcmp(name,"Phi2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRPhi2(LALpars,value);
  else if (strcmp(name,"Phi3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRPhi3(LALpars,value);
  else if (strcmp(name,"Phi4")==0)
    XLALSimInspiralWaveformParamsInsertNonGRPhi4(LALpars,value);
  else if (strcmp(name,"DChi0")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi0(LALpars,value);
  else if (strcmp(name,"DChi1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi1(LALpars,value);
  else if (strcmp(name,"DChi2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi2(LALpars,value);
  else if (strcmp(name,"DChi3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi3(LALpars,value);
  else if (strcmp(name,"DChi4")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi4(LALpars,value);
  else if (strcmp(name,"DChi5")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi5(LALpars,value);
  else if (strcmp(name,"DChi5L")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi5L(LALpars,value);
  else if (strcmp(name,"DChi6")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi6(LALpars,value);
  else if (strcmp(name,"DChi6L")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi6L(LALpars,value);
  else if (strcmp(name,"DChi7")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDChi7(LALpars,value);
  else if (strcmp(name,"DXi1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi1(LALpars,value);
  else if (strcmp(name,"DXi2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi2(LALpars,value);
  else if (strcmp(name,"DXi3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi3(LALpars,value);
  else if (strcmp(name,"DXi4")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi4(LALpars,value);
  else if (strcmp(name,"DXi5")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi5(LALpars,value);
  else if (strcmp(name,"DXi6")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDXi6(LALpars,value);
  else if (strcmp(name,"DSigma1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDSigma1(LALpars,value);
  else if (strcmp(name,"DSigma2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDSigma2(LALpars,value);
  else if (strcmp(name,"DSigma3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDSigma3(LALpars,value);
  else if (strcmp(name,"DSigma4")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDSigma4(LALpars,value);
  else if (strcmp(name,"DAlpha1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDAlpha1(LALpars,value);
  else if (strcmp(name,"DAlpha2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDAlpha2(LALpars,value);
  else if (strcmp(name,"DAlpha3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDAlpha3(LALpars,value);
  else if (strcmp(name,"DAlpha4")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDAlpha4(LALpars,value);
  else if (strcmp(name,"DAlpha5")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDAlpha5(LALpars,value);
  else if (strcmp(name,"DBeta1")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDBeta1(LALpars,value);
  else if (strcmp(name,"DBeta2")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDBeta2(LALpars,value);
  else if (strcmp(name,"DBeta3")==0)
    XLALSimInspiralWaveformParamsInsertNonGRDBeta3(LALpars,value);
  else if (strcmp(name,"aPPE")==0 || strcmp(name,"alphaPPE")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE(LALpars,value);
  else if (strcmp(name,"bPPE")==0 ||  strcmp(name,"betaPPE")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE(LALpars,value);
  else if (strcmp(name,"aPPE0")==0 || strcmp(name,"alphaPPE0")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE0(LALpars,value);
  else if (strcmp(name,"bPPE0")==0 ||  strcmp(name,"betaPPE0")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE0(LALpars,value);
  else if (strcmp(name,"alPPE1")==0 || strcmp(name,"alphaPPE1")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE1(LALpars,value);
  else if (strcmp(name,"bPPE1")==0 ||  strcmp(name,"betaPPE1")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE1(LALpars,value);
  else if (strcmp(name,"alPPE2")==0 || strcmp(name,"alphaPPE2")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE2(LALpars,value);
  else if (strcmp(name,"bPPE2")==0 ||  strcmp(name,"betaPPE2")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE2(LALpars,value);
  else if (strcmp(name,"alPPE3")==0 || strcmp(name,"alphaPPE3")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE3(LALpars,value);
  else if (strcmp(name,"bPPE3")==0 ||  strcmp(name,"betaPPE3")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE3(LALpars,value);
  else if (strcmp(name,"alPPE4")==0 || strcmp(name,"alphaPPE4")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE4(LALpars,value);
  else if (strcmp(name,"bPPE4")==0 ||  strcmp(name,"betaPPE4")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE4(LALpars,value);
  else if (strcmp(name,"alPPE5")==0 || strcmp(name,"alphaPPE5")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE5(LALpars,value);
  else if (strcmp(name,"bPPE5")==0 ||  strcmp(name,"betaPPE5")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE5(LALpars,value);
  else if (strcmp(name,"alPPE6")==0 || strcmp(name,"alphaPPE6")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE6(LALpars,value);
  else if (strcmp(name,"bPPE6")==0 ||  strcmp(name,"betaPPE6")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE6(LALpars,value);
  else if (strcmp(name,"alPPE7")==0 || strcmp(name,"alphaPPE7")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRAlphaPPE7(LALpars,value);
  else if (strcmp(name,"bPPE7")==0 ||  strcmp(name,"betaPPE7")==0 )
    XLALSimInspiralWaveformParamsInsertNonGRBetaPPE7(LALpars,value);
  else {
    XLALPrintError("XLAL Error - %s: Invalid nonGR param name %s\n",  __func__, name);
    XLAL_ERROR(XLAL_EINVAL);
  }

  return XLAL_SUCCESS;
}

/** @} */
