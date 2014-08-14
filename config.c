/*  config.c -- plugin utility functions.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <string.h>
#include <htslib/kstring.h>
#include "config.h"

char *config_get_string(const char *haystack, char *needle)
{
    kstring_t str = {0,0,0};
    ksprintf(&str,"%s=", needle);
    char *ret;
    while ( *haystack && (ret = strstr(haystack,str.s)) )
    {
        if ( !ret ) break;
        if ( ret!=haystack && ret[-1]!=':' )
        {
            // shared prefix
            haystack = ret+1;
            continue;
        }
        ret += str.l;
        char *se = ret;
        while ( *se && *se!=':' ) se++;
        str.l = 0;
        kputsn(ret,se-ret,&str);
        return str.s;
    }
    free(str.s);
    return NULL;
}

char **config_get_list(const char *opts, char *key, int *_n)
{
    char *list = config_get_string(opts,key);
    if ( !list ) return NULL;
    char *ss = list, **out;
    int n = 1;
    while ( *ss )
    {
        if ( *ss==',' ) n++;
        ss++;
    }
    out = (char**) malloc(sizeof(*out)*n);
    n = 0;
    out[n++] = list;
    ss = list;
    while ( *ss )
    {
        if ( *ss==',' )
        {
            *ss = 0;
            out[n++] = ss+1;
        }
        ss++;
    }
    *_n = n;
    return out;
}

