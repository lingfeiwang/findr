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
// This file contains OS specific routines
#ifndef _HEADER_LIB_OS_H_
#define _HEADER_LIB_OS_H_

// OS dependent new line
#ifdef __unix__
#ifdef _NEWLINE_
#undef _NEWLINE_
#endif
#define _NEWLINE_	"\n"
#endif
#ifdef __linux__
#ifdef _NEWLINE_
#undef _NEWLINE_
#endif
#define _NEWLINE_	"\n"
#endif
#ifdef __APPLE__
#ifdef _NEWLINE_
#undef _NEWLINE_
#endif
#define _NEWLINE_	"\n"
#endif
#ifndef _NEWLINE_
#error Unknown OS
#endif

#endif
