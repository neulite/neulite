// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

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

  net -> spike = calloc ( net -> n -> n_neuron, sizeof ( int ) );

  return net;
}

void finalize_network ( network_t *net )
{
  free ( net -> spike  );
  fclose ( net -> s_dat );
  fclose ( net -> v_dat );
  finalize_synapse    ( net -> s );
  finalize_connection ( net -> c );
  finalize_ion        ( net -> i );
  finalize_neuron     ( net -> n );
  finalize_population ( net -> u );
  free ( net );
}

void set_current ( const int t_ms, network_t *net, double ( *current ) ( const int, const int ) )
{
  for ( int i = 0; i < net -> n -> n_neuron; i++ ) { net -> n -> i_ext [ net -> n -> sid [ i ] ] = current ( t_ms, i ); }
}

void solve_network ( const int t_ms, network_t *net, solver_t *solver )
{
  const neuron_t *n = net -> n;

  double *v_hist = calloc ( n -> n_neuron * INV_DT, sizeof ( double ) );
  
  for ( int i = 0; i < n -> n_neuron; i++ ) {
    const int sid = n -> sid [ i ];
    double v_prev = n -> v [ sid ];
    int spike = 0;
    for ( int iter = 0; iter < INV_DT; iter++ ) {
      v_hist [ iter + INV_DT * i ] = n -> v [ sid ];
      solve ( i, net -> u, net -> n, net -> i, net -> c, net -> s, solver );
      spike += ( v_prev <= SPIKE_THRESHOLD && n -> v [ sid ] > SPIKE_THRESHOLD );
      v_prev = n -> v [ sid ];
    }
    net -> spike [ i ] = ( spike > 0 );
  }

  for ( int iter = 0; iter < INV_DT; iter++ ) {
    fprintf ( net -> v_dat, "%f ", t_ms + DT * iter );
    for ( int i = 0; i < n -> n_neuron; i++ ) {
      if ( isnan ( v_hist [ iter + INV_DT * i ] ) ) { fprintf ( stderr, "nan: %d\n", i ); exit ( 1 ); }
      fprintf ( net -> v_dat, "%f%s", v_hist [ iter + INV_DT * i ], ( i == n -> n_neuron - 1 ) ? "\n" : " " );
    }
  }
  
  free ( v_hist );
}

void spike_propagation ( const int t_ms, network_t *net )
{
  for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
    if ( net -> spike [ i ] ) { fprintf ( net -> s_dat, "%d %d\n", t_ms, i ); }
  }

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
