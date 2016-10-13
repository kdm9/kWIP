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

#include <stdarg.h>
#include <unistd.h>
#include <string.h>

#include "clogged.h"

/*******************************************************************************
*                            Internal definitions                             *
*******************************************************************************/

clg_entry_t *clg_entry_create(void);
int clg_entry_init(clg_entry_t *entry, clg_loglevel_t level,
                       const char *message);
void clg_entry_clear(clg_entry_t *entry);


int clg_entry_format(clg_entry_t *entry, clg_loglevel_t level,
                         const char *format, ...);
int clg_entry_format_va(clg_entry_t *entry,
                            clg_loglevel_t level, const char *format,
                            va_list args);
int clg_logger_write_entry(clg_logger_t *logger,
                           clg_entry_t *entry);
void _clg_entry_destroy(clg_entry_t *log_entry);
#define clg_entry_destroy(l) ({ _clg_entry_destroy(l); l = NULL; })

#define clg_free(data)              \
    do {                            \
    if (data != NULL) {             \
        free(data);                 \
        data = NULL;                \
    }                               \
    } while (0)

static inline char *
clg_strdup(const char *src)
{
    size_t len = strlen(src);
    char *dst = malloc(len + 1);
    memcpy(dst, src, len);
    dst[len] = '\0';
    return dst;
}

/*******************************************************************************
*                                Logger Funcs                                 *
*******************************************************************************/

clg_logger_t *
clg_logger_create(void)
{
    return calloc(1, sizeof(clg_logger_t));
}

int
clg_logger_init(clg_logger_t   *logger,
                clg_loglevel_t  level)
{
    if (logger == NULL) return 1;

    logger->level = level;
    return 0;
}

int clg_logger_default(clg_logger_t *logger, clg_loglevel_t level)
{
    int ret = 0;
    if (logger == NULL) return 1;
    logger->level = level;

    clg_formatter_t fmtr;
    if (isatty(STDERR_FILENO)) {
        fmtr = &clg_formatter_pretty;
    } else {
        fmtr = &clg_formatter_plain;
    }

    ret = clg_logger_add_destination_formatted(logger, stderr, level, fmtr);

    return ret;
}

int clg_logger_set_level(clg_logger_t *logger, clg_loglevel_t level)
{
    if (logger == NULL) return 1;
    logger->level = level;
    return 0;
}

int
clg_logger_add_destination_formatted(clg_logger_t      *logger,
                                     FILE              *stream,
                                     clg_loglevel_t     level,
                                     clg_formatter_t    formatter)
{
    clg_dest_t *new = NULL;
    size_t new_sz = logger->n_destinations + 1;

    new = realloc(logger->destinations,
                  new_sz * sizeof(*logger->destinations));
    if (new == NULL) {
        return 1;
    }
    logger->destinations = new;
    logger->n_destinations = new_sz;
    /* For ease of reference, save the ptr to the (new) final struct in
     * the array */
    new = &new[new_sz - 1];
    new->stream = stream;
    new->level = level;
    new->formatter = formatter;
    return 0;
}

void
_clg_logger_destroy(clg_logger_t  *logger)
{
    if (logger != NULL) {
        clg_free(logger->destinations);
        clg_free(logger);
    }
}

clg_entry_t *
clg_log_entry_create(void)
{
    return calloc(1, sizeof(clg_entry_t));
}

int
clg_log_entry_init(clg_entry_t        *entry,
                   clg_loglevel_t      level,
                   const char         *message)
{
    if (entry == NULL || message == NULL) return -1;

    entry->level = level;
    entry->message = clg_strdup(message);
    entry->formatted = NULL;
    return 0;
}

int
clg_log_entry_format_va(clg_entry_t    *entry,
                        clg_loglevel_t  level,
                        const char     *format,
                        va_list         args)
{
    int res = 0;
    char message[4096] = "";

    /* Format the error message w/ user input */
    res = vsnprintf(message, 4096, format, args);
    if (res < 1) {
        return 1;
    }
    /* Make the entry struct */
    res = clg_log_entry_init(entry, level, message);
    return res;
}

int
clg_log_entry_format(clg_entry_t  *entry,
                     clg_loglevel_t     level,
                     const char            *format,
                     ...)
{
    va_list args;
    int res = 0;

    /* Format the error message w/ user input */
    va_start(args, format);
    res = clg_log_entry_format_va(entry, level, format, args);
    va_end(args);
    return res;
}

void clg_log_entry_clear(clg_entry_t *entry)
{
    if (entry != NULL) {
        clg_free(entry->message);
        clg_free(entry->formatted);
        entry->level = CLG_LOG_DEBUG;
    }
}

void
_clg_log_entry_destroy(clg_entry_t    *entry)
{
    if (entry != NULL) {
        clg_log_entry_clear(entry);
        clg_free(entry);
    }
}

int
clg_logger_write_entry(clg_logger_t      *logger,
                       clg_entry_t   *entry)
{
    size_t iii;
    int res;

    if (logger == NULL || entry == NULL) return 1;

    /* Message is to unimportant for this logger */
    if (logger->level > entry->level) return 0;

    for (iii = 0; iii < logger->n_destinations; iii++) {
        clg_dest_t *dest = &logger->destinations[iii];

        /* Message is to unimportant for this destination */
        if (dest->level > entry->level) continue;

        dest->formatter(entry);
        res = fprintf(dest->stream, "%s", entry->formatted);
        fflush(dest->stream);
        if (res < 0) return 1;
    }
    return 0;
}

int
clg_log_msg(clg_logger_t      *logger,
            clg_loglevel_t     level,
            const char        *message)
{
    clg_entry_t entry;
    int res = 0;

    res = clg_log_entry_format(&entry, level, "%s", message);
    if (res != 0) return res;
    res = clg_logger_write_entry(logger, &entry);
    clg_log_entry_clear(&entry);
    return res;
}

int
clg_log_fmt(clg_logger_t       *logger,
            clg_loglevel_t      level,
            const char         *format,
            ...)
{
    clg_entry_t entry;
    va_list args;
    int res = 0;

    va_start(args, format);
    res = clg_log_entry_format_va(&entry, level, format, args);
    va_end(args);
    if (res != 0) return res;
    res = clg_logger_write_entry(logger, &entry);
    clg_log_entry_clear(&entry);
    return res;
}

void
clg_formatter_plain(clg_entry_t *entry)
{
    /* In the plain-text case, we just pass the message as is. */
    if (entry == NULL || entry->message == NULL) return;
    entry->formatted = clg_strdup(entry->message);
}

void
clg_formatter_pretty(clg_entry_t *entry)
{
    char newmsg[4096] = "";
    const char *colour = CLG_ANSIRST;
    const char *reset = CLG_ANSIRST;
    int res = 0;

    if (entry == NULL || entry->message == NULL) return;

    if (entry->level <= CLG_LOG_DEBUG) {
        colour = CLG_ANSIBEG CLG_ATDIM CLG_SEP CLG_FGCYN CLG_ANSIEND;
    } else if (entry->level <= CLG_LOG_PROGRESS) {
        colour = CLG_ANSIBEG CLG_ATDIM CLG_SEP CLG_FGGRN CLG_ANSIEND;
    } else if (entry->level <= CLG_LOG_INFO) {
        colour = CLG_ANSIBEG CLG_ATNRM CLG_SEP CLG_FGCYN CLG_ANSIEND;
    } else if (entry->level <= CLG_LOG_WARNING) {
        colour = CLG_ANSIBEG CLG_ATULN CLG_SEP CLG_FGYEL CLG_ANSIEND;
    } else if (entry->level <= CLG_LOG_ERROR) {
        colour = CLG_ANSIBEG CLG_ATBLD CLG_SEP CLG_FGRED CLG_ANSIEND;
    }
    res = snprintf(newmsg, 4096, "%s%s%s", colour, entry->message, reset);
    if (res > 0) {
        //clg_free(entry->formatted);
        entry->formatted = clg_strdup(newmsg);
    }
}
