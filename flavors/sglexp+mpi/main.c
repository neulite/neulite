// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h> // for exit
#include <mpi.h>
#include "network.h"
#include "solver.h"
#include "config.h"

extern int get_global_n_neurons ( const char * );
extern double get_time ( void );

double constant_current ( const double t, const int i ) { return ( I_DELAY <= t && t < I_DELAY + I_DURATION ) ? I_AMP : 0.0; /* UNIT: pA [BMTK] */ }

int main ( int argc, char *argv [ ] )
{
  if ( argc < 3 ) { fprintf ( stderr, "usage: %s <population_csv> <connection_csv>\n", argv [ 0 ] ); exit ( 1 ); }
  
  MPI_Init ( &argc, &argv );
  int mpi_size; MPI_Comm_size ( MPI_COMM_WORLD, &mpi_size );
  int mpi_rank; MPI_Comm_rank ( MPI_COMM_WORLD, &mpi_rank );
  char mpi_host [ 1024 ]; { int len; MPI_Get_processor_name ( mpi_host, &len ); mpi_host [ len ] = '\0'; }
  fprintf ( stderr, "Hello from %s: rank %d of %d\n", mpi_host, mpi_rank, mpi_size );
  
  const int global_n_neurons = get_global_n_neurons ( argv [ 1 ] );
  if ( mpi_size > global_n_neurons ) {
    if ( mpi_rank == 0 ) { fprintf ( stderr, "Error: mpi_size (%d) > global_n_neurons (%d)\n", mpi_size, global_n_neurons ); }
    MPI_Finalize ( );
    exit ( 1 );
  }

  network_t *n = initialize_network ( mpi_size, mpi_rank, argv [ 1 ], argv [ 2 ] );
  solver_t *s  = initialize_solver  ( n -> u );
  
  const double timer_start = get_time ( );
  for ( int iter = 0; iter < TSTOP * INV_DT; iter++ ) {
    const double t = DT * iter;
    if ( iter % INV_DT == 0 ) { fprintf ( stderr, "t = %f\n", t ); }
    
    output_v ( t, n );
    set_current ( t, n, constant_current );
    solve_network ( n, s );
    spike_detection ( n );
    spike_propagation ( t, n );
  }
  const double timer_stop = get_time ( );
  if ( mpi_rank == 0 ) { fprintf ( stderr, "Elapsed time = %f sec.\n", timer_stop - timer_start ); }
  
  finalize_solver  ( s );
  finalize_network ( n );
  MPI_Finalize ( );
}
