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
#include "config.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "gsl/errno.h"
#include "random.h"
#include "logger.h"
#include "lib.h"

#define MACROSTR(X)	#X
#define STR(X)	MACROSTR(X)
#define VERSION1_S	STR(VERSION1)
#define VERSION2_S	STR(VERSION2)
#define VERSION3_S	STR(VERSION3)
#define LIBVERSION	VERSION1_S "." VERSION2_S "." VERSION3_S
#define LIBNAME	STR(LIB_NAME)

#ifndef LIBINFO 
#define	LIBINFONAME(X)	X
#else
#define	LIBINFONAME(X)	LIBINFO##X
#endif


void LIBINFONAME(lib_init)(unsigned char loglv,unsigned long rs0,size_t nthread)
{
	unsigned long	rs;
	size_t	nth;
	LOGLV(loglv);
	random_init();
	rs=rs0?rs0:(unsigned long)time(NULL);
	random_seed(rs);
	if(nthread)
		omp_set_num_threads((int)nthread);
	omp_set_nested(0);
	nth=(size_t)omp_get_max_threads();
	gsl_set_error_handler_off();
	LOG(7,"Library started with log level %u, initial random seed %lu, and max thread count "PRINTFSIZET".",loglv,rs,nth)
}

const char* LIBINFONAME(lib_name)()
{
	return LIBNAME;
}

const char* LIBINFONAME(lib_version)()
{
	return LIBVERSION;
}














































