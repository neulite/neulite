// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <math.h> // isnan
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "config.h"

network_t *initialize_network ( const char *population_file, const char *connection_file )
{
  network_t *net = calloc ( 1, sizeof ( network_t ) );
  net -> u = initialize_population ( population_file );
  net -> n = initialize_neuron     ( net -> u );
  net -> i = initialize_ion        ( net -> n );
  net -> c = initialize_connection ( net -> u, net -> n, connection_file );
  net -> s = initialize_synapse    ( net -> c );

  net -> v_dat = fopen ( "v.dat", "w" );
  net -> s_dat = fopen ( "s.dat", "w" );

  net -> v_prev = calloc ( net -> n -> n_neuron, sizeof ( double ) );
  memset ( net -> v_prev, -100.0, net -> n -> n_neuron * sizeof ( double ) );

  net -> spike = calloc ( net -> n -> n_neuron, sizeof ( int ) );

  net -> tick = 0;
  net -> inv_dt = INV_DT;
  
  return net;
}

void finalize_network ( network_t *net )
{
  free ( net -> spike  );
  free ( net -> v_prev );
  fclose ( net -> s_dat );
  fclose ( net -> v_dat );
  finalize_synapse    ( net -> s );
  finalize_connection ( net -> c );
  finalize_ion        ( net -> i );
  finalize_neuron     ( net -> n );
  finalize_population ( net -> u );
  free ( net );
}

void output_v ( const double t, network_t *net )
{
  neuron_t *n = net -> n;

  fprintf ( net -> v_dat, "%f ", t );
  for ( int i = 0; i < n -> n_neuron; i++ ) {
    if ( isnan ( n -> v [ n -> sid [ i ] ] ) ) { fprintf ( stderr, "nan: %d\n", i ); exit ( 1 ); }
    fprintf ( net -> v_dat, "%f%s", n -> v [ n -> sid [ i ] ], ( i == n -> n_neuron - 1 ) ? "\n" : " " );
  }
}

void set_current ( const double t, network_t *net, double ( *current ) ( const double, const int ) )
{
  for ( int i = 0; i < net -> n -> n_neuron; i++ ) { net -> n -> i_ext [ net -> n -> sid [ i ] ] = current ( t, i ); }
}

void solve_network ( network_t *net, solver_t *solver )
{
  solve ( net -> u, net -> n, net -> i, net -> c, net -> s, solver );
}

void spike_detection ( network_t *net )
{ 
  const neuron_t *n = net -> n;

  for ( int i = 0; i < n -> n_neuron; i++ ) {
    net -> spike [ i ] += ( net -> v_prev [ i ] <= SPIKE_THRESHOLD && n -> v [ n -> sid [ i ] ] > SPIKE_THRESHOLD ); // This must be += but not =, because spike detection is made for every 1 ms.
    net -> v_prev [ i ] = n -> v [ n -> sid [ i ] ];
  }

  if ( net -> tick % net -> inv_dt == 0 ) {
    for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
      if ( net -> spike [ i ] ) { fprintf ( net -> s_dat, "%f %d\n", net -> tick * DT, i ); }
    }
  }
}

void spike_propagation ( const double t, network_t *net )
{
  if ( net -> tick % net -> inv_dt == 0 ) {
    add_spike_to_synapse_per_ms ( net -> c, net -> s ); // Add spike after delayed period is over

    // Set delay for propagation
    int size_spiking_neurons = 0;
    int *spiking_neurons = calloc ( net -> n -> n_neuron, sizeof ( int ) );
    
    for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
      if ( net -> spike [ i ] ) { spiking_neurons [ size_spiking_neurons++ ] = i; }
    }
    memset ( net -> spike, 0, net -> n -> n_neuron * sizeof ( int ) ); // net -> spike is no longer necessary
    
    const conn_t *c = net -> c;
    synapse_t *s = net -> s;
    int neuron_idx = 0, table_idx = 0;
    while ( neuron_idx < size_spiking_neurons && table_idx < c -> n_pre ) {
      if        ( spiking_neurons [ neuron_idx ] < c -> pre_table [ table_idx ] ) {
	neuron_idx++;
      } else if ( spiking_neurons [ neuron_idx ] > c -> pre_table [ table_idx ] ) {
	table_idx++;
      } else {
	for ( int j = c -> ptr_pre [ table_idx ]; j < c -> ptr_pre [ table_idx + 1 ]; j++ ) {
	  s -> delay [ c -> id [ j ] ] = ( 1 << c -> delay [ j ] );
	}
	table_idx++;
	neuron_idx++;
      }
    }
    free ( spiking_neurons );
  }
  net -> tick++;
}
