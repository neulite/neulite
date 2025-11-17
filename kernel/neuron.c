// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include "popl.h"
#include "neuron.h"
#include "config.h"

extern double ca_init ( void );

static neuron_t *initialize ( const population_t *u )
{
  neuron_t *n = calloc ( 1, sizeof ( neuron_t ) );

  int nc = 0; for ( int i = 0; i < u -> n_popl; i++ ) { nc += u -> n_neuron [ i ] * u -> n_comp [ i ]; } // nc == number of all compartments
  
  n -> v  = calloc ( nc, sizeof ( double ) );
  n -> ca = calloc ( nc, sizeof ( double ) );
  n -> i_ext = calloc ( nc, sizeof ( double ) );

  int n_neuron = 0; for ( int i = 0; i < u -> n_popl; i++ ) { n_neuron += u -> n_neuron [ i ]; }
  n -> sid = calloc ( n_neuron, sizeof ( int ) );
  n -> pid = calloc ( n_neuron, sizeof ( int ) );
  n -> n_neuron = n_neuron;

  int idx = 0, offset = 0;
  for ( int pid = 0; pid < u -> n_popl; pid++ ) {
    const int n_comp = u -> n_comp [ pid ];
    const double *vl = &u -> vl [ u -> cid [ pid ] ];
    for ( int i = 0; i < u -> n_neuron [ pid ]; i++ ) {
      n -> sid [ idx ] = offset;
      n -> pid [ idx ] = pid;
      for ( int j = 0; j < n_comp; j++ ) {
	n -> v     [ offset + j ] = vl [ j ];
	n -> ca    [ offset + j ] = ca_init ( );
	n -> i_ext [ offset + j ] = 0.0;
      }
      idx++;
      offset += n_comp;
    }
  }

  return n;
}

neuron_t *initialize_neuron ( const population_t *u ) { return initialize ( u ); }

void finalize_neuron ( neuron_t *n )
{
  if ( n -> v     != NULL ) { free ( n -> v  ); }
  if ( n -> ca    != NULL ) { free ( n -> ca ); }
  if ( n -> i_ext != NULL ) { free ( n -> i_ext ); }
  if ( n -> sid   != NULL ) { free ( n -> sid ); }
  if ( n -> pid   != NULL ) { free ( n -> pid ); }
  if ( n -> sim   != NULL ) { finalize_population ( n -> sim ); }
  free ( n );
}
