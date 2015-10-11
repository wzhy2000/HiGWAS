// fm_logger.h: interface for the logger
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FM_LOGGER_H_)
#define _FM_LOGGER_H_

#include <stdio.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

//output the contents into log file and console
void _log_prompt( const char* szSrc, int nSrcLine, const char* log_text, ... );
//output the contents into log file
void _log_info( const char* szSrc, int nSrcLine, const char* log_text, ... );
//only record the contents in the debug mode
void _log_debug( const char* szSrc, int nSrcLine, const char* log_text, ... );
//record the contents with th error mark.
void _log_error( const char* szSrc, int nSrcLine, const char* log_text, ... );
//record the contents and then close theapplication
void _log_fatal( const char* szSrc, int nSrcLine, const char* log_text, ... );
//unused now.
void _log_hexdump( const char *text, const char *buf, size_t len );

int start_log(int nDebug);


void log_set_pid( int pid );
FILE* get_log_stream();
const char* get_log_name(void);
int get_log_errorcount( int clear);

#define _HI_ __FILE__, __LINE__
#ifdef __cplusplus
}
#endif

#endif


