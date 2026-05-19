// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <math.h> // isnan
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "config.h"

#include <assert.h>
#include <stdint.h>

#define SPIKE_DATA_BITS (32U) // uint32_t
#define ID_BIT_WIDTH ( SPIKE_DATA_BITS - ITR_BIT_WIDTH ) // 24
#define ITR_BIT_OFFSET (ID_BIT_WIDTH)
#define ID_BIT_MASK  ( (1U << ID_BIT_WIDTH) - 1U )

#define PACK_ID_ITR( id, itr ) ( ((uint32_t) (itr) << ITR_BIT_OFFSET) | ((uint32_t) (id)) )
#define GET_ITR(x) ( (x) >> ITR_BIT_OFFSET )
#define GET_ID(x)  ( (x) & ID_BIT_MASK )

void assert_pack_spike_data( const int neuron_num ){
  static_assert( INV_DT < (1U << ITR_BIT_WIDTH), "INV_DT must be less than 2^ITR_BIT_WIDTH" );
  assert( neuron_num < (uint32_t) ID_BIT_MASK );
}

network_t *initialize_network ( const char *population_file, const char *connection_file )
{
  network_t *net = calloc ( 1, sizeof ( network_t ) );
  net -> u = initialize_population ( population_file );
  net -> n = initialize_neuron     ( net -> u ); assert_pack_spike_data( net -> n -> n_neuron );
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
      if ( v_prev <= SPIKE_THRESHOLD && n -> v [ sid ] > SPIKE_THRESHOLD ) {
        spike = iter + 1;
      }
      v_prev = n -> v [ sid ];
    }
    net -> spike [ i ] = ( spike );
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
    if ( net -> spike [ i ] ) { fprintf ( net -> s_dat, "%lf %d\n", (double)t_ms + DT*(double)(net->spike[i] - 1), i ); }
  }

  // Set delay for propagation
  int size_spiking_neurons = 0;
  uint32_t *spike_data =  calloc ( (unsigned int) net -> n -> n_neuron, sizeof ( uint32_t ) );

  for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
    if ( net -> spike [ i ] ) {
      spike_data [ size_spiking_neurons++ ] = PACK_ID_ITR( i, net->spike[i] );
    }
  }
  memset ( net -> spike, 0, net -> n -> n_neuron * sizeof ( int ) ); // net -> spike is no longer necessary

  const conn_t *c = net -> c;
  synapse_t *s = net -> s;
  int neuron_idx = 0, table_idx = 0;
  while ( neuron_idx < size_spiking_neurons && table_idx < c -> n_pre ) {
    if        ( GET_ID ( spike_data [ neuron_idx ] ) < c -> pre_table [ table_idx ] ) {
      neuron_idx++;
    } else if ( GET_ID ( spike_data [ neuron_idx ] ) > c -> pre_table [ table_idx ] ) {
      table_idx++;
    } else {
      uint32_t itr = GET_ITR ( spike_data [neuron_idx] );
      uint32_t id  = GET_ID  ( spike_data [neuron_idx] );
      for ( int j = c -> ptr_pre [ table_idx ]; j < c -> ptr_pre [ table_idx + 1 ]; j++ ) {
        s -> delay [ c -> id [ j ] ] = ( c -> delay [ j ] - (int)(1.f/DT) + itr );
      }
      table_idx++;
      neuron_idx++;
    }
  }
  free ( spike_data );
}
