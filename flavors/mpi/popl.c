// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "popl.h"
#include "popl_func.h"
#include "ion.h"
#include "config.h"

#define MIN(a,b) ( ( (a) < (b) ) ? (a) : (b) )

extern int get_lines ( const char * );

int get_global_n_neurons ( const char *filename )
{
  int n_popl = get_lines ( filename );
  int n_neuron [ n_popl ], n_comp [ n_popl ];
  get_population_size ( filename, n_popl, n_neuron, n_comp );

  int global_n_neurons = 0;
  for ( int i = 0; i < n_popl; i++ ) { global_n_neurons += n_neuron [ i ]; }

  return global_n_neurons;
}

static population_t *initialize ( const int n_each, const int n_offset, const int n_popl, const int n_neuron [ ], const int n_comp [ ] )
{
  int global_n_neurons = 0;
  for ( int i = 0; i < n_popl; i++ ) { global_n_neurons += n_neuron [ i ]; }

  population_t *u = calloc ( 1, sizeof ( population_t ) );

  if ( n_offset >= global_n_neurons ) { u -> n_popl = 0; return u; }

  // Here, note that calculations of start_pid and end_pid are asymmetric.
  int start_pid = -1, end_pid = -1;
  {
    int acc = 0;
    for ( int i = 0; i < n_popl; i++ ) {
      start_pid = i;
      acc += n_neuron [ i ]; // Add first
      if ( n_offset < acc ) { break; }
    }
    for ( int i = start_pid; i < n_popl; i++ ) {
      end_pid = i;
      if ( MIN ( n_offset + n_each, global_n_neurons ) <= acc ) { break; } // Break first
      acc += n_neuron [ i + 1 ];
    }
  }
  assert ( 0 <= start_pid && start_pid < n_popl );
  assert ( 0 <= end_pid   && end_pid < n_popl );

  const int local_n_popl = end_pid - start_pid + 1;
  u -> n_popl = local_n_popl; //n_popl;
  u -> n_neuron = calloc ( u -> n_popl, sizeof ( int ) );
  u -> n_comp = calloc ( u -> n_popl, sizeof ( int ) );
  u -> start_pid = start_pid;
  u -> end_pid   = end_pid;

  {
    int acc = 0; for ( int i = 0; i < start_pid; i++ ) { acc += n_neuron [ i ]; }
    int n_rest = ( n_offset + n_each < global_n_neurons ) ? n_each : global_n_neurons - n_offset;
    int pid = start_pid;
    while ( n_rest > 0 ) {
      const int res = (pid == start_pid ) ? acc - n_offset + n_neuron [ pid ] : n_neuron [ pid ];
      u -> n_neuron [ pid - start_pid ] = MIN( res, n_rest );
      n_rest -= MIN( res, n_rest );
      acc += n_neuron [ pid ];
      pid++;
    }
  }

  u -> cid = calloc ( ( u -> n_popl + 1 ), sizeof ( int ) ); // Cumulative compartment id
  u -> cid [ 0 ] = 0;
  for ( int i = 0; i < u -> n_popl; i++ ) {
    u -> n_comp [ i ] = n_comp [ start_pid + i ];
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

population_t *initialize_population ( const int n_each, const int n_offset, const char *filename )
{
  int n_popl = get_lines ( filename );
  int n_neuron [ n_popl ], n_comp [ n_popl ];
  get_population_size ( filename, n_popl, n_neuron, n_comp );
  
  // Below, population size is of local
  population_t *u = initialize ( n_each, n_offset, n_popl, n_neuron, n_comp );

  if ( u -> n_popl == 0 ) { return u; }
  
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
    if ( u -> start_pid <= pid && pid <= u -> end_pid ) {
      const int local_pid = pid - u -> start_pid;
      read_swc_file ( u, local_pid, swcfile );
      read_pas_file ( u, local_pid, ionfile );
      read_ion_file ( u, local_pid, ionfile );
    }
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
