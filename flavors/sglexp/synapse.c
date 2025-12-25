// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "conn.h"
#include "synapse.h"
#include "config.h"

synapse_t *initialize_synapse ( conn_t *c )
{
  synapse_t *s = calloc ( 1, sizeof ( synapse_t ) );
  s -> sum0   = calloc ( c -> n_conn, sizeof ( double ) );
  s -> delay  = calloc ( c -> n_conn, sizeof ( int ) );
  return s;
}

void finalize_synapse ( synapse_t *s )
{
  if ( s -> sum0  != NULL ) { free ( s -> sum0  ); }
  if ( s -> delay != NULL ) { free ( s -> delay ); }
  free ( s );
}

void update_synapse ( const conn_t * __restrict__ c, synapse_t * __restrict__ s )
{
  for ( int i = 0; i < c -> n_conn; i++ ) { s -> sum0 [ i ] *= c -> decay [ i ]; }
}

void add_spike_to_synapse_per_ms ( const conn_t * __restrict__ c, synapse_t * __restrict__ s ) // each 1 ms
{
  for ( int i = 0; i < c -> n_conn; i++ ) {
    s -> delay [ i ] >>= 1;
    s -> sum0 [ i ] += ( s -> delay [ i ] == 1 ) ? 1 : 0;
  }
}
