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
// This file contains the logger and error functions of different levels
#ifndef _HEADER_LIB_LOGGER_H_
#define _HEADER_LIB_LOGGER_H_
#include "config.h"
#include <stdlib.h>
#include <stdarg.h>
#include "os.h"
#ifdef __cplusplus
extern "C"
{
#endif

// global variable name for struct logger
#define LOGGER_VARIABLE logger_variable
// Logging macro. logs with significance level LV, and the rest are in printf format.
#define LOGS(LOGGERX,LV,...)	logger_log(LOGGERX,LV,__FILE__,__LINE__,__VA_ARGS__);
#define LOG(LV,...)	LOGS(&LOGGER_VARIABLE,LV,__VA_ARGS__)
#define LOGLV(LV)	LOGGER_VARIABLE.lv=LV
/* Logging levels:
 * CRITICAL(0),ERROR(1),ERROR(2),ERROR(3),WARNING(4),WARNING(5),WARNING(6),INFO(7),INFO(8),INFO(9),DEBUG(10),DEBUG(11),DEBUG(12)
 */

struct logger{
// 	logger output level. Only message levels<lv will be ouput: 0:critical, 1-3: error, 4-6: warning, 7-9: info, 10-12: debug
	size_t lv;
};

// Return the name of message level lv
const char* logger_mname(size_t lv);

// Outputs log with level lv to stderr, formatted as "Log level name:file name:line number: user defined format newline"
// file gives file name, line give line number, fmt gives user defined format, ... gives (printf) parameters of fmt
void logger_voutput(size_t lv,const char* file,size_t line,const char* fmt,va_list args);

// Similar with logger_voutput, ... version
void logger_output(size_t lv,const char* file,size_t line,const char* fmt,...);

// Logs message. Similar with logger_output, but checks if output level<=logger level.
// Return 0 on output, or 1 on not output because of high output level.
int logger_log(const struct logger* l,size_t lv,const char* file,size_t line,const char* fmt,...);

//	Initialize logger l to level lv. l must be unreferenced. Messages <=l are output to stderr.
int logger_init(struct logger* l,size_t lv);

//	Initialize default logger (LOGGER_VARIABLE) with logger_init
int logger_default_init(size_t lv);
//	Allocate new logger at level lv.
struct logger* logger_new(size_t lv);

extern struct logger LOGGER_VARIABLE;

#ifdef __cplusplus
}
#endif
#endif
