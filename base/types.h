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
// This file contains the numerical type definitions, especially vectors and matrices

#ifndef _HEADER_LIB_TYPES_H_
#define _HEADER_LIB_TYPES_H_
#include "config.h"
#include <float.h>
#include "gsl/vector.h"
#include "gsl/matrix.h"
#include "gsl/blas.h"


#if FTYPEBITS == 32
	// Type definition
	#define FTYPE	float
	// Type suffix definition, for gsl vector and matrix functions
	#define	FTYPE_SUF	_float
	// BLAS function macro
	#define BLASF(X)	BLASFO(X)
	// Minimal value
	#define FTYPE_MIN	FLT_MIN
#elif FTYPEBITS == 64
	#define FTYPE	double
	#define	FTYPE_SUF	
	#define BLASF(X)	BLASFD(X)
	#define FTYPE_MIN	DBL_MIN
#else
	#error Unknown float type bit count.
#endif
#if GTYPEBITS == 8
	#define	GTYPE		unsigned char
	#define	GTYPE_SUF	_uchar
#else
	#error Unknown genotype type bit count.
#endif
#define BLASFO(X)	gsl_blas_s ## X
#define BLASFD(X)	gsl_blas_d ## X

#define CONCATENATE2_(X,Y)		X ## Y
#define CONCATENATE2(X,Y)		CONCATENATE2_(X,Y)
#define CONCATENATE3_(X,Y,Z)	X ## Y ## Z
#define CONCATENATE3(X,Y,Z)		CONCATENATE3_(X,Y,Z)
#define CONCATENATE4_(X,Y,Z,W)	X ## Y ## Z ## W
#define CONCATENATE4(X,Y,Z,W)	CONCATENATE4_(X,Y,Z,W)

// vector type macro
#define VECTORO		gsl_vector_float
#define VECTORD		gsl_vector
#define VECTORC		gsl_vector_char
#define VECTORUC	gsl_vector_uchar
#define VECTORI		gsl_vector_int
#define VECTORL		gsl_vector_long
#define VECTORUL	gsl_vector_ulong
#define VECTORF		CONCATENATE2(gsl_vector,FTYPE_SUF)
#define	VECTORG		CONCATENATE2(gsl_vector,GTYPE_SUF)
// vector function type macro
#define VECTOROF(X)		gsl_vector_float_ ## X
#define VECTORDF(X)		gsl_vector_ ## X
#define VECTORCF(X)		gsl_vector_char_ ## X
#define VECTORUCF(X)	gsl_vector_uchar_ ## X
#define VECTORIF(X)		gsl_vector_int_ ## X
#define VECTORLF(X)		gsl_vector_long_ ## X
#define VECTORULF(X)	gsl_vector_ulong_ ## X
#define VECTORFF(X)		CONCATENATE2(VECTORF,_ ## X)
#define VECTORGF(X)		CONCATENATE2(VECTORG,_ ## X)
// matrix type macro
#define MATRIXO		gsl_matrix_float
#define MATRIXD		gsl_matrix
#define MATRIXC		gsl_matrix_char
#define MATRIXUC	gsl_matrix_uchar
#define MATRIXI		gsl_matrix_int
#define MATRIXL		gsl_matrix_long
#define MATRIXUL	gsl_matrix_ulong
#define MATRIXF		CONCATENATE2(gsl_matrix,FTYPE_SUF)
#define MATRIXG		CONCATENATE2(gsl_matrix,GTYPE_SUF)
// matrix function type macro
#define MATRIXOF(X)		gsl_matrix_float_ ## X
#define MATRIXDF(X)		gsl_matrix_ ## X
#define MATRIXCF(X)		gsl_matrix_char_ ## X
#define MATRIXUCF(X)	gsl_matrix_uchar_ ## X
#define MATRIXIF(X)		gsl_matrix_int_ ## X
#define MATRIXLF(X)		gsl_matrix_long_ ## X
#define MATRIXULF(X)	CONCATENATE2(gsl_matrix_ulong_ ## X
#define MATRIXFF(X)		CONCATENATE2(MATRIXF,_ ## X)
#define MATRIXGF(X)		CONCATENATE2(MATRIXG,_ ## X)




















































#endif
