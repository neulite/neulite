// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "popl.h"
#include "neuron.h"
#include "ion.h"
#include "conn.h"
#include "synapse.h"
#include "solver.h"
#include "hines.h"
#include "config.h"
//#include <omp.h>

// for Hines solver manipulation
extern hines_matrix_t *hines_matrix_initialize ( const population_t *, const int );
extern void           *hines_matrix_finalize ( hines_matrix_t * );

solver_t *initialize_solver ( const population_t *u )
{
  const int n_popl = u -> n_popl;
  int n_neuron = 0; for ( int i = 0; i < n_popl; i++ ) { n_neuron += u -> n_neuron [ i ]; }
  
  solver_t *solver = calloc ( 1, sizeof ( solver_t ) );
  solver -> linsys = calloc ( n_neuron, sizeof ( linsys_t ) );
  solver -> n_popl = n_popl;
  solver -> n_neuron = n_neuron;

  int offset = 0;
  for ( int pid = 0; pid < n_popl; pid++ ) {
    const int n_comp = u -> n_comp  [ pid ];

    // Matrix
    double *mat = calloc ( n_comp * n_comp, sizeof ( double ) );
    {
      const int offset = u -> cid [ pid ];
      double *rad = &u -> rad [ offset ];
      double *len = &u -> len [ offset ];
      double *ra  = &u -> ra  [ offset ];
      int *parent = &u -> parent [ offset ];

      for ( int i = 0; i < n_comp * n_comp; i++ ) { mat [ i ] = 0.0; }
      for ( int i = 0; i < n_comp; i++ ) {
	const int d = parent [ i ];
	if ( d >= 0 ) {
	  double r = ( 2.0 / ( ( ra [ i ] * len [ i ] ) / ( rad [ i ] * rad [ i ] * M_PI ) + ( ra [ d ] * len [ d ] ) / ( rad [ d ] * rad [ d ] * M_PI ) ) ); // -1 * [mS]
	  mat [ d + n_comp * i ] = r;
	  mat [ i + n_comp * d ] = mat [ d + n_comp * i ]; // i*NGO -> set Rows, +d -> set Columns
	}
      }
      for ( int i = 0; i < n_comp; i++ ) {
	double r = 0;
	for ( int j = 0; j < n_comp; j++ ) {
	  r += mat [ j + n_comp * i ];
	}
	mat [ i + n_comp * i ] = r;
      }
    }

    for ( int li = 0; li < u -> n_neuron [ pid ]; li++ ) {
      linsys_t *s = &solver -> linsys [ offset + li ];
      s -> H = hines_matrix_initialize( u, pid );
      for ( int i = 0; i < n_comp; i++ ) {
	const int parent_id = s -> H -> parent_id [ i ];
	s -> H -> Ad     [ i ] = mat [ i + n_comp * i ];
	s -> H -> Api    [ i ] = ( parent_id >= 0 ) ? -mat [ parent_id + n_comp * i ] : 0;
	s -> H -> bu_Ad  [ i ] = mat [ i + n_comp * i ];
	s -> H -> bu_Api [ i ] = ( parent_id >= 0 ) ? -mat [ parent_id + n_comp * i ] : 0;
      }
      s -> b = calloc ( n_comp, sizeof ( double ) ); // b value
    }

    offset += u -> n_neuron [ pid ];
    free ( mat );
  }

  return solver;
}

static void update_matrix ( const int id, const population_t * __restrict__ u, const neuron_t * __restrict__ n, const ion_t * __restrict__ i, const conn_t * __restrict__ c, const synapse_t * __restrict__ s, linsys_t * __restrict__ linsys, const double dt )
{
  const int sid = n -> sid [ id ];
  const int pid = n -> pid [ id ];
  const double *v      = &n -> v [ sid ];
  const double *ca     = &n -> ca [ sid ];
  const double *i_ext  = &n -> i_ext [ sid ];
  const double *cm     = &u -> cm [ u -> cid [ pid ] ];
  const double *g_leak = &u -> gl [ u -> cid [ pid ] ];
  const double *v_leak = &u -> vl [ u -> cid [ pid ] ];
  const int n_comp = u -> n_comp [ pid ];

  for ( int li = 0; li < linsys -> H -> n_comp; li++ ) { linsys -> H -> Ad [ li ] = linsys -> H -> bu_Ad [ li ]; linsys -> H  -> Api [ li ] = linsys -> H -> bu_Api[ li ]; }

  for ( int li = 0; li < n_comp; li++ ) {
    linsys -> H -> Ad[ li ] += ( cm [ li ] / dt ) + g_leak [ li ];
    linsys -> b [ li ]  = ( cm [ li ] / dt ) * v [ li ] + g_leak [ li ] * v_leak [ li ] + i_ext [ li ] * 1e-3; /* CONVERSION: 1e-3 from pA to nA */
  }
  
  double lhs = 0.0, rhs = 0.0;
  calc_lhs_and_rhs ( u, n, i, pid, id, &lhs, &rhs );
  linsys -> H -> Ad[ 0  ] += lhs;
  linsys -> b [ 0  ] += rhs;

  for ( int li = c -> ptr_post [ id ]; li < c -> ptr_post [ id + 1 ]; li++ ) {
    const int post_c    = c -> post_c [ li ];
    const double weight = c -> weight [ li ];
    const double erev   = c -> erev   [ li ];
    const double sum0   = s -> sum0   [ li ];
    linsys -> b [ post_c ] += - weight * 1e-6 * ( sum0 ) * ( v [ post_c ] - erev ); /* CONVERSION: 1e-6 (unitless?) */
  }
}

static void solve_matrix ( linsys_t * __restrict__ l )
{
  int n_comp = l -> H -> n_comp;
  double *Ad = l -> H -> Ad;
  double *Api = l -> H -> Api;
  int *parent_id = l -> H -> parent_id;
  double *b = l -> b;
  double *x = l -> b;
  
  // TRIANG
  for ( int i = n_comp - 1; i > 0; i-- ) {
    int pid = parent_id [ i ];
    Ad [ pid ] -= Api [ i ] * Api [ i ] / Ad [ i ]; // A(i,p) = A(p,i)
    b  [ pid ] -= b   [ i ] * Api [ i ] / Ad [ i ];
  }
  
  // FWSUB
  x [ 0 ] = b [ 0 ] / Ad [ 0 ];
  for ( int i = 1; i < n_comp; i++ ) {
    int pid = parent_id [ i ];
    x [ i ] = ( b [ i ] - x [ pid ] * Api [ i ] ) / Ad [ i ];
  }
}

void solve ( const population_t * __restrict__ u, neuron_t * __restrict__ n, ion_t * __restrict__ i, const conn_t * __restrict__ c, synapse_t * __restrict__ s, solver_t * __restrict__ solver )
{
  update_synapse ( c, s );

  //#pragma omp parallel for
  for ( int li = 0; li < n -> n_neuron; li++ ) {
    const int sid = n -> sid [ li ];
    const int pid = n -> pid [ li ];
    const int n_comp = u -> n_comp [ pid ];
    linsys_t *linsys = &solver -> linsys [ li ];
    update_matrix ( li, u, n, i, c, s, linsys, 0.5*DT );
    update_ion ( li, n, i, DT );
    update_ca ( li, u, i, n, DT );
    solve_matrix ( linsys );
    for ( int j = 0; j < n_comp; j++ ) { n -> v [ sid + j ] = 2 * linsys -> b [ j ] - n -> v [ sid + j ]; }
  }
}

void finalize_solver ( solver_t *solver )
{
  for ( int i = 0; i < solver -> n_neuron; i++ ) {
    linsys_t *s = &solver -> linsys [ i ];
    hines_matrix_finalize ( s -> H );
    free ( s -> H );
    free ( s -> b );
  }
  free ( solver -> linsys );
  free ( solver );
}
