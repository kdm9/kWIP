/* Copyright © 2016 Kevin Murray <kdmfoss@gmail.com>
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef CLOGGED_H_UH1VXAPO
#define CLOGGED_H_UH1VXAPO

#include <stdlib.h>
#include <stdio.h>


typedef enum clg_loglevel {
    /* The idea is that the user can add values between these, if they need
     * to. a la Python's logging module. */
    CLG_LOG_DEBUG = 0, /* Debug gets compiled out according to NDEBUG */
    CLG_LOG_PROGRESS = 5,
    CLG_LOG_INFO = 10,
    CLG_LOG_WARNING = 20,
    CLG_LOG_ERROR = 30,
} clg_loglevel_t;

typedef struct clg_entry {
    char *message;
    char *formatted;
    clg_loglevel_t level;
} clg_entry_t;

typedef void (*clg_formatter_t)(clg_entry_t *);

typedef struct clg_dest {
    FILE *stream;
    clg_loglevel_t level;
    clg_formatter_t formatter;
} clg_dest_t;

typedef struct clg_logger {
    clg_dest_t *destinations;
    size_t n_destinations;
    clg_loglevel_t level;
    int lock;
} clg_logger_t;


/*******************************************************************************
*                              Logger Functions                               *
*******************************************************************************/

clg_logger_t *clg_logger_create(void);
int clg_logger_init(clg_logger_t *logger, clg_loglevel_t level);
int clg_logger_set_level(clg_logger_t *logger, clg_loglevel_t level);
int clg_logger_add_destination_formatted(clg_logger_t *logger,
                                         FILE *stream,
                                         clg_loglevel_t level,
                                         clg_formatter_t formatter);
#define clg_logger_add_destination(log, stream, level)                      \
    clg_logger_add_destination_formatted(log, stream, level,                \
                                         &clg_format_plain)

/* Two default formatters, pretty is colourised, plain is not */
void clg_formatter_plain(clg_entry_t *entry);
void clg_formatter_pretty(clg_entry_t *entry);

/* Default setup: log to stderr only, use colours if isatty(stderr) */
int clg_logger_default(clg_logger_t *logger, clg_loglevel_t level);

void _clg_logger_destroy(clg_logger_t *logger);
#define clg_logger_destroy(l) ({ _clg_logger_destroy(l); l = NULL; })


/*******************************************************************************
*                              User Entry Points                              *
*******************************************************************************/

int clg_log_msg(clg_logger_t *logger, clg_loglevel_t level,
                    const char *message);
#ifndef NDEBUG
#define clg_log_msg_debug(log, msg_) clg_log_msg(log, CLG_LOG_DEBUG, msg_)
#else
#define clg_log_msg_debug(log, msg_)
#endif
#define clg_log_msg_progress(log, msg_) clg_log_msg(log, CLG_LOG_PROGRESS, msg_)
#define clg_log_msg_info(log, msg_) clg_log_msg(log, CLG_LOG_INFO, msg_)
#define clg_log_msg_warning(log, msg_) clg_log_msg(log, CLG_LOG_WARNING, msg_)
#define clg_log_msg_error(log, msg_) clg_log_msg(log, CLG_LOG_ERROR, msg_)


int clg_log_fmt(clg_logger_t *logger, clg_loglevel_t level,
                   const char *format, ...);
#ifndef NDEBUG
#define clg_log_fmt_debug(log, fmt_, ...) \
        clg_log_fmt(log, CLG_LOG_DEBUG, fmt_, __VA_ARGS__)
#else
#define clg_log_fmt_debug(log, fmt_, ...)
#endif
#define clg_log_fmt_info(log, fmt_, ...) \
        clg_log_fmt(log, CLG_LOG_INFO, fmt_, __VA_ARGS__)
#define clg_log_fmt_progress(log, fmt_, ...) \
        clg_log_fmt(log, CLG_LOG_PROGRESS, fmt_, __VA_ARGS__)
#define clg_log_fmt_warning(log, fmt_, ...) \
        clg_log_fmt(log, CLG_LOG_WARNING, fmt_, __VA_ARGS__)
#define clg_log_fmt_error(log, fmt_, ...) \
        clg_log_fmt(log, CLG_LOG_ERROR, fmt_, __VA_ARGS__)

/*******************************************************************************
*                                 ANSI codes                                  *
*******************************************************************************/
/* If you want to make a custom formatter, use these colour codes */

#define CLG_ANSIBEG  "\033["
#define CLG_ANSIEND  "m"

#define CLG_ANSIRST  CLG_ANSIBEG "0" CLG_ANSIEND

#define CLG_CLRLINE  "\033[1G\033[2K"

#define CLG_SEP  ";"
#define CLG_ATNRM  "0"
#define CLG_ATBLD  "1"
#define CLG_ATDIM  "2"
#define CLG_ATULN  "3"
#define CLG_ATBNK  "5"
#define CLG_ATREV  "7"
#define CLG_ATHID  "8"

#define CLG_FGBLK  "30"
#define CLG_FGRED  "31"
#define CLG_FGGRN  "32"
#define CLG_FGYEL  "33"
#define CLG_FGBLU  "34"
#define CLG_FGMAG  "35"
#define CLG_FGCYN  "36"
#define CLG_FGWHT  "37"

#define CLG_BGBLK  "40"
#define CLG_BGRED  "41"
#define CLG_BGGRN  "42"
#define CLG_BGYEL  "43"
#define CLG_BGBLU  "44"
#define CLG_BGMAG  "45"
#define CLG_BGCYN  "46"
#define CLG_BGWHT  "47"


#endif /* end of include guard: CLOGGED_H_UH1VXAPO */
