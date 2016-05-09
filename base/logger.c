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
#include "config.h"
#include <stdio.h>
#include <time.h>
#include "os.h"
#include "logger.h"

struct logger LOGGER_VARIABLE;

const char* logger_mname(size_t lv)
{
	static const char names[13][15]={"CRITICAL(0)","ERROR(1)","ERROR(2)","ERROR(3)","WARNING(4)","WARNING(5)","WARNING(6)","INFO(7)","INFO(8)","INFO(9)","DEBUG(10)","DEBUG(11)","DEBUG(12)"};
	if(lv>12)
		return 0;
	return names[lv];
}

void logger_voutput(size_t lv,const char* file,size_t line,const char* fmt,va_list args)
{
	char		timing[100];
	struct tm	str_time;
	time_t 		rawtime;
	
	time(&rawtime);
	localtime_r(&rawtime,&str_time);
	strftime(timing,99,"%F %T",&str_time);

	logprintf("%s:%s:%s:%lu: ",logger_mname(lv),timing,file,line);
	logvprintf(fmt,args);
	logprintf("%s",_NEWLINE_);
}

void logger_output(size_t lv,const char* file,size_t line,const char* fmt,...)
{
	va_list args;
	va_start (args, fmt);
	logger_voutput(lv,file,line,fmt,args);
}

int logger_log(const struct logger* l,size_t lv,const char* file,size_t line,const char* fmt,...)
{
	va_list args;
	va_start (args, fmt);
	if(lv>l->lv)
		return 1;
	logger_voutput(lv,file,line,fmt,args);
	return 0;
}

int logger_init(struct logger* l,size_t lv)
{
	if(!l)
	{
		logger_output(1,__FILE__,__LINE__,"NULL logger.");
		return 1;
	}
	l->lv=lv;
	return 0;
}

int logger_default_init(size_t lv)
{
	return logger_init(&LOGGER_VARIABLE,lv);
}

struct logger* logger_new(size_t lv)
{
	struct logger* l;
	l=(struct logger*)calloc(1,sizeof(struct logger));
	if(!l)
	{
		logger_output(1,__FILE__,__LINE__,"Logger allocation failed.");
		return 0;
	}
	if(logger_init(l,lv))
	{
		free(l);
		return 0;
	}
	return l;
}
