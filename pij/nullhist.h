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
/* This file contains the methods to construct null histogram.
 * Each method contains two parts: an interface function to produce
 * the null histogram and a struct for the parameters needed by the
 * function.
 * Currently two methods are implemented.
 * 1:	Sample-model method utilizes sampler functions
 *		to randomly sample null data and modeler functions to model the
 *		null histogram based on sampled null data. This is a generic
 *		method and can be tailored for specific needs.
 * 2:	Analytical method specifically for pij model, where the distribution
 *		of LLR of null hypothesis is exactly calculable. This is
 *		the actually used method for faster and more precise performance.
 * 
 * Regardless of method, the interface function should conform with
 * the following definition:
 * (*int)(const void* p,gsl_histogram* h);
 * p:	Parameter of the method.
 * h:	Existing histogram with range specified. The function should
 * 		fill this histogram in the bins.
 * Return:	0 on success.
 */

#ifndef _HEADER_LIB_PIJ_NULLHIST_H_
#define _HEADER_LIB_PIJ_NULLHIST_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#include "nullsampler.h"
#include "nullmodeler.h"
#ifdef __cplusplus
extern "C"
{
#endif

/****************************************************************
 * The sample-model method for generation of null histogram
 * Uses a sampler function to randomly sample null data and a modeler
 * to model the histogram of null data. The samplers are defined
 * in pij/nullsampler.h. The modelers are defined in
 * pij/nullmodeler.h. The sample-model method is generic
 * but for the specific pij model, analytical method (see below)
 * perform more precise and faster.
 ***************************************************************/

struct pij_nullhist_sample_model_param
{
	//Input data
	const void *d;
	/* Data struct size (for d).
	 * Must be small and use pointers, for memory saving and easier partitioning.
	 */
	size_t		dsize;
	/* Partition function to partition mission/data into smaller chunks.
	 * Look threading_threading_get_startend for help
	 * src:		Larger data to be partitioned, having same type as d
	 * dest:	Smaller data after partitioning, same type as d
	 * id:		Id of the data chuck (out of total) to obtain from partition
	 * total:	Total number of chucks for partitioning.
	 * Return:	'Size' of dest after partition. Must be !=0 for nonzero data size
	 * 			which needs further process, or 0 if no calculation is required
	 * 			due to zero size in the partition.
	 */
	size_t (*partition)(const void* src,void* dest,size_t id,size_t total);
	//Sampler function struct
	const struct pij_nullsampler* sampler;
	//Modeler function struct
	const struct pij_nullmodeler* modeler;
	//Parameter of sampler function
	const void* ps;
	//Parameter of modeler function
	const void* pm;
};

/* Null histogram generation of sample-model method on single thread.
 * g,t,t2,sampler,modeler,ps:	See pij_nullhist_sample_model_param.
 * pmc:		Model container object. See definition in nullmodel.h.
 * Return:	0 if success.
 */
int pij_nullhist_sample_model_single(const void* data,const struct pij_nullsampler* sampler,const struct pij_nullmodeler* modeler,const void* ps,void* pmc);

// Interface function of sample-model method, multithreaded.
int	pij_nullhist_sample_model(const void* p,gsl_histogram* h);


/****************************************************************
 * The analytical methods to calculate null histograms for pij
 * is faster and more precise. This method calculates pdf values at
 * evenly spreaded points within each bin for better accuracy.
 * The three steps of pij share the same analytical formula of
 * null histograms.
 ***************************************************************/


struct pij_nullhist_analytical_pdf_param
{
	/* Number of split within each bin. 2^nsplit central points are
	 * taken and the mean is calculated as the average pdf within the bin.
	 */
	size_t	nsplit;
	/* Null pdf function to draw null histogram from
	 * Param 1:	Coordinates to calculate pdfs
	 * Param 2: Output buffer for pdfs as return
	 * Param 3:	Parameter accepted by function
	 * Return:	0 on success.
	 */
	int (*func)(const VECTORD*,VECTORD*,const void*);
	//Parameter accepted by pdf function (Param 3)
	const void* param;
};

/* Generic interface function of analytical method.
 * Same as interface function except with one extra parameter in the end.
 * func:	to specifiy which function to use to calculate null pdf.
 */
int	pij_nullhist_analytical_pdf(const void* param,gsl_histogram* h);




















#ifdef __cplusplus
}
#endif
#endif
