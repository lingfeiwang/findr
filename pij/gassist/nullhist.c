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
#include "../../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include "../../base/gsl/blas.h"
#include "../../base/logger.h"
#include "../../base/macros.h"
#include "../../base/histogram.h"
#include "../nullhist.h"
#include "nullhist.h"
#include "nulldist.h"


int	pij_gassist_nullhist_analytical1_pdf(const void* param,gsl_histogram* h)
{
	const struct pij_gassist_nullhist_analytical_pdf_param	*p=param;
	const struct pij_gassist_nulldist_mixed_pdf_data		p2={p->g,p->nv};
	const struct pij_nullhist_analytical_pdf_param			p3={p->nsplit,pij_gassist_nulldist1_mixed_pdf,&p2};
	
	return pij_nullhist_analytical_pdf(&p3,h);
}

int	pij_gassist_nullhist_analytical2b_pdf(const void* param,gsl_histogram* h)
{
	const struct pij_gassist_nullhist_analytical_pdf_param	*p=param;
	const struct pij_gassist_nulldist_mixed_pdf_data		p2={p->g,p->nv};
	const struct pij_nullhist_analytical_pdf_param			p3={p->nsplit,pij_gassist_nulldist2b_mixed_pdf,&p2};
	
	return pij_nullhist_analytical_pdf(&p3,h);
}

int	pij_gassist_nullhist_analytical3_pdf(const void* param,gsl_histogram* h)
{
	const struct pij_gassist_nullhist_analytical_pdf_param	*p=param;
	const struct pij_gassist_nulldist_mixed_pdf_data		p2={p->g,p->nv};
	const struct pij_nullhist_analytical_pdf_param			p3={p->nsplit,pij_gassist_nulldist3_mixed_pdf,&p2};
	
	return pij_nullhist_analytical_pdf(&p3,h);
}






















