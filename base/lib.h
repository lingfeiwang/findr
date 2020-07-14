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
//This file contains library related functions
#ifndef _HEADER_LIB_LIB_H_
#define _HEADER_LIB_LIB_H_
#ifdef __cplusplus
extern "C"
{
#endif


/* The library needs to be initialized before any other function is called,
 * to perform correctly with desired log level and random seed.
 * loglv:	Logging level, see logger.h.
 * rs:		Initial random seed. If rs=0, use current time as random seed.
 * nthread:	Maximum number of threads, If nthread=0, use default setting.
 */
void lib_init(unsigned char loglv,unsigned long rs,size_t nthread);

/* Returns library name
 */
const char* lib_name();
/* Returns library version in a.b.c format, or a, b, or c, for subfunctions ending with 1, 2, or 3 respectively.
 */
const char* lib_version();
size_t lib_version1();
size_t lib_version2();
size_t lib_version3();

#ifdef __cplusplus
}
#endif
#endif
