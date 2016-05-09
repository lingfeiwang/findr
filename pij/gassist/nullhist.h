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
/* This file provides the analytical method to construct histograms
 * for the null pdf of LLRs of genotype assisted pij inference.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_NULLHIST_H_
#define _HEADER_LIB_PIJ_GASSIST_NULLHIST_H_
#include "../../base/config.h"
#include "../../base/gsl/histogram.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

struct pij_gassist_nullhist_analytical_pdf_param
{
	//(ng,ns) Genotype data
	const MATRIXG*	g;
	//Maximum number of possible values each genotype can take.
	//=Number of alleles + 1
	size_t	nv;
	/* Number of split within each bin. 2^nsplit central points are
	 * taken and the mean is calculated as the average pdf within the bin.
	 */
	size_t	nsplit;
};


//Specific interface function for different stages
int	pij_gassist_nullhist_analytical1_pdf(const void* param,gsl_histogram* h);
//Shared identical analytical null pdf for conservative step 2
#define pij_gassist_nullhist_analytical2c_pdf	pij_gassist_nullhist_analytical1_pdf
int	pij_gassist_nullhist_analytical2b_pdf(const void* param,gsl_histogram* h);
int	pij_gassist_nullhist_analytical3_pdf(const void* param,gsl_histogram* h);





















#ifdef __cplusplus
}
#endif
#endif
