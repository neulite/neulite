// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "popl.h"
#include "popl_func.h"
#include "ion.h"
#include "config.h"

extern int get_lines ( const char * );

static population_t *initialize ( const int n_popl, const int n_neuron [ ], const int n_comp [ ] )
{
  population_t *u = calloc ( 1, sizeof ( population_t ) );

  u -> n_popl = n_popl;
  u -> n_neuron = calloc ( u -> n_popl, sizeof ( int ) );
  u -> n_comp = calloc ( u -> n_popl, sizeof ( int ) );
  u -> cid = calloc ( ( u -> n_popl + 1 ), sizeof ( int ) ); // Cumulative compartment id
  u -> cid [ 0 ] = 0;

  for ( int i = 0; i < u -> n_popl; i++ ) {
    u -> n_neuron [ i ] = n_neuron [ i ];
    u -> n_comp   [ i ] = n_comp   [ i ];
    u -> cid [ i + 1 ]  = u -> cid [ i ] + u -> n_comp [ i ];
  }

  const int nc = u -> cid [ u -> n_popl ]; // total number of compartments

  u -> rad = calloc  ( nc, sizeof ( double ) );
  u -> len = calloc  ( nc, sizeof ( double ) );
  u -> area = calloc ( nc, sizeof ( double ) );
  u -> parent = calloc ( nc, sizeof ( int ) );
  u -> type   = calloc ( nc, sizeof ( int ) );

  u -> cm = calloc ( nc, sizeof ( double ) );
  u -> ra = calloc ( nc, sizeof ( double ) );
  u -> gl = calloc ( nc, sizeof ( double ) );
  u -> vl = calloc ( nc, sizeof ( double ) );

  u -> gamma = calloc ( n_popl * N_COMPTYPE, sizeof ( double ) );
  u -> decay = calloc ( n_popl * N_COMPTYPE, sizeof ( double ) );

  u -> gbar  = calloc ( ( ( ALLACTIVE == 1 ) ? nc : n_popl ) * N_GBAR, sizeof ( double ) ); // Default is perisomatic

  return u;
}

population_t *initialize_population ( const char *filename )
{
  int n_popl = get_lines ( filename );
  int n_neuron [ n_popl ], n_comp [ n_popl ];
  get_population_size ( filename, n_popl, n_neuron, n_comp );
  
  population_t *u = initialize ( n_popl, n_neuron, n_comp );
  
  FILE *file = fopen ( filename, "r" );
  if ( ! file ) { fprintf ( stderr, "Error: no such file %s\n", filename ); exit ( 1 ); }

  int pid = 0;
  char buf [ 1024 ];

  while ( fgets ( buf, 1024, file ) ) {
    if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
    if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
    int dn_neuron, dn_comp;
    char name [ 1024 ], swcfile [ 1024 ], ionfile [ 1024 ];
    const int nf = sscanf ( buf, "%d,%d,%[^,],%[^,],%[^,]", &dn_neuron, &dn_comp, name, swcfile, ionfile );
    assert ( nf == 5 );
    read_swc_file ( u, pid, swcfile );
    read_pas_file ( u, pid, ionfile );
    read_ion_file ( u, pid, ionfile );
    pid++;
  }

  fclose ( file );

  return u;
}

void finalize_population ( population_t *u )
{
  free ( u -> n_neuron );
  free ( u -> n_comp );
  free ( u -> cid );

  free ( u -> rad );
  free ( u -> len );
  free ( u -> area );
  free ( u -> parent );
  free ( u -> type );

  free ( u -> cm );
  free ( u -> ra );
  free ( u -> gl );
  free ( u -> vl );

  free ( u -> gbar );
  free ( u -> gamma );
  free ( u -> decay );

  free ( u );
}
