/* Copyright 2016, 2017 Lingfei Wang
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
/* This file contains the LLR sampler definition of null hypotheses.
 * Samplers are used in combination with modelers for the estimation
 * of histograms for the distribution of LLR of null hypothesis, as
 * applied in sample-model method. (See pij_nullhist_sample_model.)
 * 
 * Each null sampler is consisted of four functions and a parameter struct.
 * For function formats, see struct pij_nullsampler.
 * 
 * For multithreading, each thread holds a separate and independent sampler.
 */

#ifndef _HEADER_LIB_PIJ_NULLSAMPLER_H_
#define _HEADER_LIB_PIJ_NULLSAMPLER_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#include "../base/random.h"
#ifdef __cplusplus
extern "C"
{
#endif


struct pij_nullsampler
{
	/* Initialization function. Should be invoked before any other operation.
	 * Function should perform parameter storage and memory allocations here.
	 * Somer samplers may use only part of the parameters.
	 * Param 1:	Input data
	 * Param 2:	Sampler parameters
	 * Param 3:	Initial random seed.
	 * Return:	Sampler parameter struct pointer, or 0 if failed.
	 */
	void* (*init)(const void*,const void*,unsigned long);
	/* Obtain dimensionality of samples.
	 * Param 1: Sampler parameter struct.
	 * Param 2:	Return the return size of vector each time sampling function
	 * 			is invoked.
	 * Param 3:	Return the number of cycles sampling function should be invoked.
	 */
	void (*dimensions)(const void*,size_t*,size_t*);
	/* Produce samples by random with the stated amount by function dimensions.
	 * Param 1:	Sampler parameter struct.
	 * Param 2:	Existing buffer to write samples into.
	 */
	void (*sample)(void*,VECTORF*);
	/* Close sampler and frees memory.
	 * Param 1:	Sampler parameter struct.
	 */
	void (*close)(void*);
};

#ifdef __cplusplus
}
#endif
#endif
