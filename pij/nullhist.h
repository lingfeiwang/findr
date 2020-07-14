/* Copyright 2016-2018, 2020 Lingfei Wang
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
/* This part produces the histogram of parametric null distributions.
 */

#ifndef _HEADER_LIB_PIJ_NULLHIST_H_
#define _HEADER_LIB_PIJ_NULLHIST_H_
#include "../base/config.h"
#include "../base/types.h"
#include "../base/gsl/histogram.h"
#ifdef __cplusplus
extern "C"
{
#endif


/* Construct one null histogram for a specific genotype value count.
 * The function calculates the null density histogram for random variable:
 * x=-0.5*log(1-z1/(z1+z2)), where z1~chi2(n1),z2~chi2(n2),
 * Histogram bin count and width are automatically determined
 * from real data count (nd).
 * For bin range settings, see histogram_unequalbins_fromnullcdf.
 * For null density histogram from pdf, see pij_nulldist_hist_pdf.
 * dmax:	Specifies the histogram bound as [0,dmax).
 * nd:		Count of real data to form real histograms. This is used to
 * 			automatically decide number of bins and widths.
 * n1,
 * n2:		Parameters of null distribution.
 * Return:	Constructed null distribution histograms with preset
 * 			bin ranges and values as density.
 */
gsl_histogram* pij_nullhist_single(double dmax,size_t nd,size_t n1,size_t n2);

/* Construct multiple null histograms for different genotype value counts.
 * The function calculates the null density histogram for random variable:
 * x=-0.5*log(1-z1/(z1+z2)), where z1~chi2(i*n1c+n1d),z2~chi2(-i*n2c+n2d),
 * i=0,...,nv-2. Histogram bin count and width are automatically determined
 * from real data count (nd).
 * For bin range settings, see histogram_unequalbins_fromnullcdf.
 * For null density histogram from pdf, see pij_nulldist_hist_pdf.
 * dmax:	Specifies the histogram bound as [0,dmax).
 * nv:		Maximum number of values each genotype can type. Must be nv>=2.
 * 			This limits the possible values of kv in distribution, and
 * 			also output histogram count.
 * nd:		Count of real data to form real histograms. This is used to
 * 			automatically decide number of bins and widths.
 * n1c,
 * n1d,
 * n2c
 * n2d:		Parameters of null distribution.
 * Return:	[nv-1]. Constructed null distribution histograms with preset
 * 			bin ranges and values as density. Genotypes with i values have
 * 			histogram stored in Return[i-2].
 */
gsl_histogram** pij_nullhist(double dmax,size_t nv,size_t nd,long n1c,size_t n1d,long n2c,size_t n2d);












#ifdef __cplusplus
}
#endif
#endif
