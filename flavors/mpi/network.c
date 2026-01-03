// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <math.h> // isnan
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>
#include "network.h"
#include "config.h"

extern int get_global_n_neurons ( const char * );

network_t *initialize_network ( const int mpi_size, const int mpi_rank, const char *population_file, const char *connection_file )
{
  network_t *net = calloc ( 1, sizeof ( network_t ) );
  net -> global_n_neurons = get_global_n_neurons ( population_file );
  
  const int n_each   = ( net -> global_n_neurons + mpi_size - 1 ) / mpi_size;
  const int n_offset = n_each * mpi_rank;

  net -> u = initialize_population ( n_each, n_offset, population_file );
  net -> n = initialize_neuron     ( net -> u );
  net -> i = initialize_ion        ( net -> n );
  net -> c = initialize_connection ( n_each, n_offset, net -> u, net -> n, connection_file );
  net -> s = initialize_synapse    ( net -> c );
  
  {
    char filename [ 1024 ];
    snprintf ( filename, sizeof ( filename ), "v%d.dat", mpi_rank );
    net -> v_dat = fopen ( filename, "w" );
    snprintf ( filename, sizeof ( filename ), "s%d.dat", mpi_rank );
    net -> s_dat = fopen ( filename, "w" );
  }

  net -> spike = calloc ( net -> n -> n_neuron, sizeof ( int ) );

  net -> mpi_size = mpi_size;
  net -> mpi_rank = mpi_rank;

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
  // Detailed
  //const int n_each   = ( net -> global_n_neurons + net -> mpi_size - 1 ) / net -> mpi_size;
  //const int n_offset = n_each * net -> mpi_rank;
  //for ( int i = n_offset; i < MIN(n_offset + n_each, net -> global_n_neurons); i++ ) {
  //  net -> n -> i_ext [ net -> n -> sid [ i - n_offset ] ] = current ( t_ms, i );
  //}
  // Simple
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
  const int n_each   = ( net -> global_n_neurons + net -> mpi_size - 1 ) / net -> mpi_size;
  const int n_offset = n_each * net -> mpi_rank;
  const int mpi_size = net -> mpi_size;
  
  for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
    if ( net -> spike [ i ] ) { fprintf ( net -> s_dat, "%d %d\n", t_ms, n_offset + i ); }
  }

  add_spike_to_synapse_per_ms ( net -> c, net -> s ); // Add spike after delayed period is over

  //
  // Broadcast spikes across nodes
  //
  int *local_spiking_neurons = malloc ( n_each * sizeof ( int ) );
  int local_count = 0;
  for ( int i = 0; i < net -> n -> n_neuron; i++ ) {
    if ( net -> spike [ i ] ) { local_spiking_neurons [ local_count++ ] = n_offset + i; }
  }
  memset ( net -> spike, 0, net -> n -> n_neuron * sizeof ( int ) ); // net -> spike is no longer necessary
    
  // Gather # of neurons that emitted a spike on each process
  int *spike_counts = malloc ( mpi_size * sizeof ( int ) );
  MPI_Allgather ( &local_count, 1, MPI_INT, spike_counts, 1, MPI_INT, MPI_COMM_WORLD );
    
  // Calculate offsets
  int *displs = malloc ( mpi_size * sizeof ( int ) );
  displs [ 0 ] = 0;
  for ( int i = 1; i < mpi_size; i++ ) {
    displs [ i ] = displs [ i - 1 ] + spike_counts [ i - 1 ];
  }

  // Gather the neuron IDs
  int total_spikes = displs [ mpi_size - 1 ] + spike_counts [ mpi_size - 1 ];
  int *global_spiking_neurons = malloc ( total_spikes * sizeof ( int ) );
  MPI_Allgatherv ( local_spiking_neurons, local_count, MPI_INT, global_spiking_neurons, spike_counts, displs, MPI_INT, MPI_COMM_WORLD );

  // Copy to spiking_neurons for backward compatibility
  const int size_spiking_neurons = total_spikes;
  int *spiking_neurons = malloc ( ( size_spiking_neurons ? size_spiking_neurons : 1 ) * sizeof ( int ) );
  if ( size_spiking_neurons > 0 ) {
    memcpy(spiking_neurons, global_spiking_neurons, size_spiking_neurons * sizeof ( int ) );
  }
    
  // Free memory
  free ( local_spiking_neurons );
  free ( spike_counts );
  free ( displs );
  free ( global_spiking_neurons ); 

  //
  // Spike propagation
  //
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
