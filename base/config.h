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
//This is the general configuration header

#ifndef _HEADER_LIB_CONFIG_H_
#define _HEADER_LIB_CONFIG_H_
#include "config_auto.h"
#ifdef LIBEXTENSION_R
//#include <R_ext/Print.h>
//#define logprintf REprintf
//#define	logvprintf REvprintf
#define logprintf(...) fprintf(stderr,__VA_ARGS__)
#define logvprintf(...) vfprintf(stderr,__VA_ARGS__)
#ifdef	NDEBUG
#undef	NDEBUG
#endif
#ifdef	GSL_RANGE_CHECK_OFF
#undef	GSL_RANGE_CHECK_OFF
#endif
#ifdef	HAVE_INLINE
#undef	HAVE_INLINE
#endif
#define	NDEBUG	1
#define GSL_RANGE_CHECK_OFF	1
#define HAVE_INLINE	1
#else
#define logprintf(...) fprintf(stderr,__VA_ARGS__)
#define logvprintf(...) vfprintf(stderr,__VA_ARGS__)
#endif
#endif
