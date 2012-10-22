/*
 *
 *  LALInferenceLikelihood.c:   Likelihood functions for LALInference codes        
 *  LALInferenceLikelihood.h:   header file
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

/**
 * \file SMEELikelihood.h
 * \brief Header file for likelihood functions used by LALInference codes
 *
 * LALInferenceLikelihood contains all the necessary routines to compute the likelihood 
 * from a template (computed with LALInferenceTemplate) and the data (initialised with LALInferenceReadData).
 *
 * Likelihood functions follow the basic naming convention: LALInference<type_of>LogLikelihood()
 *
 * Takes as input:
 * - a pointer to a LALInferenceVariable structure containing the parameters to compute the likelihood for,
 * - a pointer to a LALInferenceIFOData structure containing the linked list of interferometer data,
 * - a pointer to the LALInferenceTemplateFunction template function to be used.
 * 
 * Outputs as a REAL8 the natural logarithm value of the likelihood, as defined by:
 *  
 * \f[
 * Likelihood(\vec{x}|\vec{\lambda},M)=\exp(-\tfrac{1}{2}<\vec{x}-\vec{h_M}(\vec{\lambda})|\vec{x}-\vec{h_M}(\vec{\lambda})>)
 * \f] 
 *
 * where: \f$<x|y>=4Re\left ( \int \frac{\tilde{x}\,\tilde{y}^*}{S_f}\, df \right )\f$
 * 
 *
 * Note that the likelihood is reported unnormalised.
 * 
 */



#ifndef SMEELikelihood_h
#define SMEELikelihood_h

#include <lal/LALInference.h>



/***********************************************************//**
 * (log-) likelihood function.                                 
 * Returns the non-normalised logarithmic likelihood.          
 * Slightly slower but cleaner than							   
 * UndecomposedFreqDomainLogLikelihood().          `		   
 *
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, >0)                      
 *   - "time"            (REAL8, GPS sec.)                     
 ***************************************************************/


			      
REAL8 LALInferenceSMEEFreqDomainLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData * data,
                              LALInferenceTemplateFunction *template);			      

/***********************************************************//**
 * Frequency-domain single-IFO response computation.           
 * Computes response for a given template.                    
 * Will re-compute template only if necessary                  
 * (i.e., if previous, as stored in data->freqModelhCross,     
 * was based on different parameters or template function).    
 * Carries out timeshifting for a given detector               
 * and projection onto this detector.                          
 * Result stored in freqResponse, assumed to be correctly      
 * initialized												   
 *
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, >0)                      
 *   - "time"            (REAL8, GPS sec.)                     
 ***************************************************************/				

			      
void LALInferenceComputeSMEEFreqDomainResponse(LALInferenceVariables *currentParams, LALInferenceIFOData * dataPtr,
                              LALInferenceTemplateFunction *template, COMPLEX16Vector *freqWaveform);	
			      
void get_model( LALInferenceVariables *currentParams, LALInferenceIFOData *dataPtr, COMPLEX16Vector *freqWaveform);			      


/**
 * Computes the <x|y> overlap in the Fourrier domain.
 */
REAL8 LALInferenceComputeFrequencyDomainOverlap(LALInferenceIFOData * data,
        COMPLEX16Vector * freqData1, COMPLEX16Vector * freqData2);

/**
 * Identical to LALInferenceFreqDomainNullLogLikelihood, but returns the likelihood of a null template.
 * Used for normalising.
 */
REAL8 LALInferenceNullLogLikelihood(LALInferenceIFOData *data);

REAL8 LALInferenceSMEENullLogLikelihood(LALInferenceIFOData *data);

/***********************************************************//**
 * Student-t (log-) likelihood function                        
 * as described in Roever/Meyer/Christensen (2011):            
 *   "Modelling coloured residual noise                        
 *   in gravitational-wave signal processing."                 
 *   Classical and Quantum Gravity, 28(1):015010.              
 *   http://dx.doi.org/10.1088/0264-9381/28/1/015010           
 *   http://arxiv.org/abs/0804.3853                            
 * Returns the non-normalised logarithmic likelihood.          
 * 
 * Required (`currentParams') parameters are:                  
 *   - "rightascension"  (REAL8, radian, 0 <= RA <= 2pi)       
 *   - "declination"     (REAL8, radian, -pi/2 <= dec <=pi/2)  
 *   - "polarisation"    (REAL8, radian, 0 <= psi <= ?)        
 *   - "distance"        (REAL8, Mpc, > 0)                     
 *   - "time"            (REAL8, GPS sec.)                     
 * 
 * This function is essentially the same as the                
 * "UndecomposedFreqDomainLogLikelihood()" function.           
 * The additional parameter to be supplied is the (REAL8)      
 * degrees-of-freedom parameter (nu) for each Ifo.             
 * The additional "df" argument gives the corresponding        
 * d.f. parameter for each element of the "*data" list.        
 * The names of "df" must match the "->name" slot of           
 * the elements of "data".                                     
 *                                                             
 * (TODO: allow for d.f. parameter to vary with frequency,     
 *        i.e., to be a set of vectors corresponding to        
 *        frequencies)                                         
 ***************************************************************/
REAL8 LALInferenceFreqDomainStudentTLogLikelihood(LALInferenceVariables *currentParams, LALInferenceIFOData *data,
                                      LALInferenceTemplateFunction *template);

#endif
