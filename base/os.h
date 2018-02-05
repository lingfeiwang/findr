/* Copyright 2016-2018 Lingfei Wang
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
// This file contains OS specific routines
#ifndef _HEADER_LIB_OS_H_
#define _HEADER_LIB_OS_H_

#ifdef _NEWLINE_
#undef _NEWLINE_
#endif

// OS dependent new line
#if defined(unix) || defined(__unix__) || defined(__unix) || defined(__APPLE__) || defined(__MACH__) || defined(__linux__)
#define _NEWLINE_	"\n"
#define	PRINTFSIZET	"%zu"
#endif
#if defined(_WIN32) || defined(_WIN64)
#define _NEWLINE_	"\r\n"
#define	PRINTFSIZET	"%Iu"
#endif
#ifndef _NEWLINE_
#error Unsupported OS
#endif

#endif
