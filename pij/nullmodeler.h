/* Copyright 2016 Lingfei Wang
 * 
 * This file is part of Findr.
 * 
 * Findr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Findr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with Findr.  If not, see <http://www.gnu.org/licenses/>.
 */
/* This file contains the modelers of LLR histograms of null hypothesis.
 * Modelers are used by sample-model method to produce null histogram.
 * (See pij_nullhist_sample_model.)
 * 
 * Each null modeler is consisted of seven functions and a parameter struct.
 * For function formats, see struct pij_nullmodeler.
 * The nullhist function will first invoke 'init_container' method to obtain
 * the unique model container object. Then for each thread, one modeler object
 * is constructed as a child of the container, via 'init' method. Each thread
 * holds a separate modeler object and use 'input' method to feed in data.
 * Upon finishing, every thread merges modeler into the model container, with
 * 'merge'. Merge takes place for one thread at a time, guaranteed by nullhist.
 * The thread will 'close' modeler before exit. When all child threads completes,
 * the master thread will apply 'output' method on the model container to obtain
 * null histogram, and then use 'close_container' method to free up container.
 *
 * Two modelers are implemented:
 * 1:	Exponential modeler assumes LLR satisfies an exponential distribution.
 * 2:	Naive modeler assumes the true LLR histogram is in exact agreement with
 * 		the sampled histogram.
 */

#ifndef _HEADER_LIB_PIJ_NULLMODELER_H_
#define _HEADER_LIB_PIJ_NULLMODELER_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif


struct pij_nullmodeler
{
	/* Container initialization function. Should be invoked before
	 * any other operation.
	 * Function should perform parameter storage and memory allocations here.
	 * Somer modelers may use only part of the parameters.
	 * Param 1: Input data
	 * Param 2:	Modeler parameter
	 * Param 3:	Histogram as output.
	 * Return:	Modeler container parameter struct pointer, or 0 if failed.
	 */
	void* (*init_container)(const void*,const void*,const gsl_histogram*);
	/* Output modeler container to histogram.
	 * Histogram should sum to 1.
	 * Param 1:	Modeler container parameter struct.
	 * Param 2:	Output histogram. Bin ranges must be identical with the one
	 * 			referenced in init.
	 */
	int (*output)(const void*,gsl_histogram*);
	/* Closes modeler container and frees memory.
	 * Param 1:	Modeler container parameter struct.
	 */
	void (*close_container)(void*);
	/* Initialization function. Should be invoked after container but
	 * before any other operation.
	 * Function should perform parameter storage and memory allocations here.
	 * Somer modelers may use only part of the parameters.
	 * Param 1:	Output model container for reference.
	 * Param 2: Input data
	 * Param 3:	n1. Sample size to be provided each time.
	 * Param 4: n2. Number of times samples will be provided.
	 * 			So n1*n2 is the total sample count.
	 * Return:	Modeler parameter struct pointer, or 0 if failed.
	 */
	void* (*init)(const void*,const void*,size_t,size_t);
	/* Input samples into modeler. Each sample should count 1.
	 * Param 1:	Modeler parameter struct.
	 * Param 2:	Existing samples.
	 */
	void (*input)(void*,const VECTORF*);
	/* Perform final calculations and merge modeler into modeler container.
	 * Function does not need to guarantee threadsafe.
	 * Param 1:	Modeler struct.
	 * Param 2:	Modeler container.
	 */
	void (*merge)(const void*,void*);
	/* Close modeler and frees memory.
	 * Param 1:	Modeler struct.
	 */
	void (*close)(void*);

};

/*************************************************************************
 * Exponential modeler believes null data yields exponential distribution.
 * This requries all data to be non-negative.
 * The modeler seeks the best exponential coefficient via maximization of
 * likelihood.
 ************************************************************************/

//Modeler state and also model container state
struct pij_nullmodeler_exp_state
{
	//Number of samples
	size_t	n;
	//Sum of samples
	FTYPE	v;
};

extern const struct pij_nullmodeler pij_nullmodeler_exp;


/*************************************************************************
 * Naive modeler assumes the pdf of falling into each bin is exactly proportional
 * to the number of instances sampled to be in that bin.
 * Therefore a simple accumulation of samples in existing bins would suffice.
 ************************************************************************/
struct pij_nullmodeler_naive_state
{
//Histogram of sampled null data.
	gsl_histogram* h;
};

extern const struct pij_nullmodeler pij_nullmodeler_naive;





#ifdef __cplusplus
}
#endif
#endif
