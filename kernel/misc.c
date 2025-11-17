// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double get_time ( void )
{
  struct timespec ts;
  clock_gettime ( CLOCK_REALTIME, &ts );
  return  ts.tv_sec + ( double ) ts.tv_nsec * 1.0e-9;
}

int remove_blank_destructive_for_csv ( char *buf )
{
  char *src = buf, *dst = buf;
  while ( *src && *src != '\n' ) { if ( *src == ' ' || *src == '\t' ) { char _ = *src++; } else { *dst++ = *src++; } }
  *dst = '\0';
  return strlen ( buf );
}

int strip_comment_destructive ( char *buf )
{
  char *src = buf, *dst = buf;
  while ( *src && *src != '\n' ) { if ( *src == '#' ) { break; } else { *dst++ = *src++; } }
  *dst = '\0';
  return strlen ( buf );
}

int get_lines ( const char *filename )
{
  FILE *file = fopen ( filename, "r" );
  if ( ! file ) { fprintf ( stderr, "Error: no such file %s\n", filename ); exit ( 1 ); }
  char buf [ 1024 ];
  int n = 0;
  while ( fgets ( buf, 1024, file ) ) {
    if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
    n++;
  }
  fclose ( file );
  return n;
}
