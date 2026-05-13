// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h> // for exit
#include "network.h"
#include "solver.h"
#include "config.h"

extern double get_time ( void );

double constant_current ( const int t_ms, const int i ) { return ( I_DELAY <= t_ms && t_ms < I_DELAY + I_DURATION ) ? I_AMP : 0.0; /* UNIT: pA [BMTK] */ }

int main ( int argc, char *argv [ ] )
{
  if ( argc < 3 ) { fprintf ( stderr, "usage: %s <population_csv> <connection_csv>\n", argv [ 0 ] ); exit ( 1 ); }
  
  network_t *n = initialize_network ( argv [ 1 ], argv [ 2 ] );
  solver_t *s  = initialize_solver  ( n -> u );
  
  const double timer_start = get_time ( );
  for ( int t_ms = 0; t_ms < TSTOP; t_ms++ ) {
    fprintf ( stderr, "t = %d\n", t_ms );
    
    set_current ( t_ms, n, constant_current );
    solve_network ( t_ms, n, s );
    spike_propagation ( t_ms, n );
  }
  const double timer_stop = get_time ( );
  fprintf ( stderr, "Elapsed time = %f sec.\n", timer_stop - timer_start );
  
  finalize_solver  ( s );
  finalize_network ( n );
}
