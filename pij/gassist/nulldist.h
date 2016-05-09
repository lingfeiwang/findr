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
/* This part contains analytical calculation of the distribution of log likelihood ratio from null hypothesis.
 * Each function is applicable to one or more stages, which are stated in the function name as pij_nulldistX_..., where X is the applicable stage.
 * For each stage, different methods to calculate histogram can coexist. The method is declared in the function name as suffix:
 * _cdf:	Calculate histogram as the difference of cdf.
 * 			This is applicable when distribution is single-variable integrable.
 * _pdf:	Calculate histogram as the pdf mean of points evenly split within the bin. This is applicable when distribution is single-variable non-integrable.
 * _sim:	Construct histogram by sampling. This is applicable when distribution is multi-variable non-integrable.
 */

#ifndef _HEADER_LIB_PIJ_GASSIST_NULLDIST_H_
#define _HEADER_LIB_PIJ_GASSIST_NULLDIST_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/******************************************************************************************************
 * Mixture distributions of LLR for any step at specific locations.
 * The distribution is a mixture one because different genotype locations can have different number
 * of values observed, yielding to different analytical distributions of LLR.
 *****************************************************************************************************/

/* Calculate log total mixture pdf of LLR for any step at specific locations with buffer provided.
 * The distribution is a mixture one because different genotype locations can have different number
 * of values observed, yielding to different analytical distributions of LLR.
 * Specify function for different steps.
 * g:	(ng,ns) Genotype data
 * nv:	Number of values each genotype can take
 * loc:	(nd) Locations of LLR of current step, where log pdf will be calculated
 * ans:	(nd) Output for log pdf calculated.
 * vb2:	(nd) Buffer
 * vb3:	(ng) Buffer
 * vb4:	(nv+1) Buffer
 * mb1:	(nv-1,nd) Saves pdf before summing up w.r.t nv.
 * nd:	loc->size
 * func:Function to calculate pdf buffed (e.g. pij_nulldist1_calcpdf_buffed or pij_nulldist3_calcpdf_buffed)
 */
void pij_gassist_nulldist_mixed_pdf_buffed(const MATRIXG* g,size_t nv,const VECTORD* loc,VECTORD* ans,VECTORD* vb2,VECTORG* vb3,VECTORD* vb4,VECTORUC* vb5,MATRIXD* mb1,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*));

/* Calculate log pdf of LLR for any step at specific locations.
 * g:	(ng,ns) Genotype data
 * nv:	Number of values each genotype can take
 * loc:	(nd) Locations of LLR of current step, where log pdf will be calculated
 * ans:	(nd) Output for log pdf calculated.
 * nd:	loc->size
 * func:Function to calculate pdf buffed (e.g. pij_nulldist1_calcpdf_buffed or pij_nulldist3_calcpdf_buffed)
 * Return:	0 on success.
 */
int pij_gassist_nulldist_mixed_pdf(const MATRIXG* g,size_t nv,const VECTORD* loc,VECTORD* ans,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*));

/* Calculates density histogram of null distribution of log likelihood ratio.
 * Is a pij_allrtopij_nullhist_method (see pij_allrtopij.h).
 * Method is to use central pdf value as the density.
 * g:		(ng,ns) Genotype data
 * nv:		Number of values each genotype can take
 * range:	(nbin+1) Bin boundary values
 * nbin:	Number of bins
 * hist:	(nbin) Output array for density histogram.
 * param:	Redundant, must be 0.
 * func:	Function to calculate pdf buffed (e.g. pij_nulldist1_calcpdf_buffed or pij_nulldist3_calcpdf_buffed)
 * Return:	0 on success.
 */
int pij_gassist_nulldist_nullhist_mixed_pdf0(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*));

/* Calculates density histogram of null distribution of log likelihood ratio.
 * Method is to use central pdf value as the density.
 * g:	(ng,ns) Genotype data
 * nv:	Number of values each genotype can take
 * range:	(nbin+1) Bin boundary values
 * nbin:	Number of bins
 * hist:	(nbin) Output array for density histogram.
 * param:	Redundant, must be 0.
 * func:	Function to calculate pdf buffed (e.g. pij_nulldist1_calcpdf_buffed or pij_nulldist3_calcpdf_buffed)
 * Return:	0 on success.
 */
int pij_gassist_nulldist_nullhist_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param,void (*func)(size_t,size_t,const VECTORD*,MATRIXD*,VECTORD*));


/*************************************************************
 * Specific functions for each step
 *************************************************************/

/* Calculate log pdf of LLR for each step at specific locations with buffer provided.
 * Uses pij_nulldist_calcpdf_buffed.
 * ns:	Number of samples
 * nv:	Number of values each genotype can take. Calculates for value counts=2,...,nv.
 * loc:	(nd) Locations of LLR of each step, where log pdf will be calculated
 * ans:	(nv-1,nd) Output matrix. ans[i,j]=p(loc[j],ns,i+2).
 * vb2:	(nd) Buffer.
 */
void pij_gassist_nulldist1_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2);
#define	pij_gassist_nulldist2_conserv_calcpdf_buffed pij_gassist_nulldist1_calcpdf_buffed
void pij_gassist_nulldist2b_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2);
// #define	pij_nulldist2_calcpdf_buffed pij_nulldist2_bold_calcpdf_buffed
void pij_gassist_nulldist3_calcpdf_buffed(size_t ns,size_t nv,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2);


struct pij_gassist_nulldist_mixed_pdf_data
{
	const MATRIXG* g;
	size_t nv;
};

/* Calculate log pdf of LLR for each step at specific locations using pij_nulldist_mixed_pdf.
 * Steps 1 and 2_conserv are identical.
 * g:	(ng,ns) Genotype data
 * nv:	Number of values each genotype can take
 * loc:	(nd) Locations of LLR of current step, where log pdf will be calculated
 * ans:	(nd) Output for log pdf calculated.
 * Return:	0 on success.
 */
int pij_gassist_nulldist1_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param);
#define pij_gassist_nulldist2c_mixed_pdf	pij_gassist_nulldist1_mixed_pdf
int pij_gassist_nulldist2b_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param);
// #define pij_nulldist2_pdf pij_nulldist2_bold_pdf
int pij_gassist_nulldist3_mixed_pdf(const VECTORD* loc,VECTORD* ans,const void* param);

// #define pij_nulldist2_nullhist_pdf0	pij_nulldist1_nullhist_pdf0
/* Calculates density histogram of null distribution of log likelihood ratio.
 * Steps 1 and 2_conserv are identical.
 * Uses pij_nulldist_nullhist_pdf.
 * Method is to use central pdf value as the density.
 * g:	(ng,ns) Genotype data
 * nv:	Number of values each genotype can take
 * range:	(nbin+1) Bin boundary values
 * nbin:	Number of bins
 * hist:	(nbin) Output array for density histogram.
 * param:	Redundant, must be 0.
 * Return:	0 on success.
 */
int pij_gassist_nulldist_nullhist1_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param);
#define pij_gassist_nulldist2c_nullhist_mixed_pdf	pij_gassist_nulldist1_nullhist_mixed_pdf
int pij_gassist_nulldist_nullhist2b_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param);
// #define pij_nulldist_nullhist2_pdf pij_nulldist_nullhist2_bold_pdf
int pij_gassist_nulldist_nullhist3_mixed_pdf(const MATRIXG* g,size_t nv,const double* restrict range,size_t nbin,double* restrict hist,void* param);


















#ifdef __cplusplus
}
#endif
#endif
