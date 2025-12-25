// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "conn.h"
#include "neuron.h"
#include "config.h"

extern int strip_comment_destructive ( char * );
extern int remove_blank_destructive_for_csv ( char * );
extern int get_lines ( const char * );

conn_t *initialize_connection ( const population_t *u, const neuron_t *n, const char *filename )
{
  conn_t *c = calloc (1, sizeof ( conn_t ) );

  int *pre_ary  = calloc ( n -> n_neuron, sizeof ( int ) );
  int *post_ary = calloc ( n -> n_neuron, sizeof ( int ) );
  int n_conn = 0, n_pre = 0;
  {
    char buf [ 1024 ] ;
    FILE *file = fopen ( filename, "r" );
    if ( ! file ) { fprintf ( stderr, "Error: no such file %s\n", filename ); exit ( 1 ); }

    while ( fgets ( buf, 1024, file ) ) {
      if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
      if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
      int d_pre, d_post_i;
      const int nf = sscanf ( buf, "%d,%d", &d_pre, &d_post_i );
      assert ( nf == 2 );
      n_conn += 2;
      pre_ary  [ d_pre    ] += 2;
      post_ary [ d_post_i ] += 2;
    }
    fclose ( file );
    for ( int i = 0; i < n -> n_neuron; i++ ) { if ( pre_ary [ i ] > 0 ) { n_pre++; } }
  }

  c -> n_conn = n_conn;
  c -> n_pre = n_pre;
  c -> n_post = n -> n_neuron; // Note: RHS is X n_post, O n -> n_neuron

  c -> pre_table  = calloc ( c -> n_pre,  sizeof ( int ) ); // list of ids of presynaptic neurons
  {
    int j = 0;
    for ( int i = 0; i < n -> n_neuron; i++ ) {
      if ( pre_ary [ i ] > 0 ) { c -> pre_table [ j ] = i; j++; }
    }
    assert ( j == c -> n_pre );
  }
  c -> ptr_pre = calloc ( c -> n_pre + 1, sizeof ( int ) );
  {
    int j = 0;
    c -> ptr_pre [ 0 ] = 0;
    for ( int i = 0; i < n -> n_neuron; i++ ) {
      if ( pre_ary [ i ] > 0 ) { c -> ptr_pre [ j + 1 ] = c -> ptr_pre [ j ] + pre_ary [ i ]; j++; }
    }
    assert ( j == c -> n_pre );
  }
  c -> ptr_post = calloc ( c -> n_post + 1, sizeof ( int ) );
  {
    c -> ptr_post [ 0 ] = 0;
    for ( int i = 0; i < c -> n_post; i++ ) {
      c -> ptr_post [ i + 1 ] = c -> ptr_post [ i ] + post_ary [ i ];
    }
  }
  free (pre_ary);
  free (post_ary);
 
  c -> post_c = calloc ( c -> n_conn, sizeof ( int ) );
  c -> weight = calloc ( c -> n_conn, sizeof ( double ) );
  c -> erev   = calloc ( c -> n_conn, sizeof ( double ) );
  c -> decay  = calloc ( c -> n_conn, sizeof ( double ) );
  c -> delay = calloc ( c -> n_conn, sizeof ( int ) );
  c -> id    = calloc ( c -> n_conn, sizeof ( int ) );

  { 
    int *local_idx = calloc ( c -> n_post, sizeof ( int ) );
    FILE *file = fopen ( filename, "r" );
    char buf [ 1024 ];
    int idx = 0;
    while ( fgets ( buf, 1024, file ) ) {
      if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
      if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
      int d_pre, d_post_i, d_post_c, d_delay;
      double f_weight, f_decay, f_rise, f_erev;
      char c_type;
      const int nf = sscanf ( buf, "%d,%d,%d,%lf,%lf,%lf,%lf,%d,%c", &d_pre, &d_post_i, &d_post_c, &f_weight, &f_decay, &f_rise, &f_erev, &d_delay, &c_type );
      assert ( nf == 9 );

      assert ( d_post_c < u -> n_comp [ n -> pid [ d_post_i ] ] );
      assert ( d_delay > 0 );
      
      const double tau_prime = f_decay * f_rise / ( f_decay - f_rise );
      const double tau_diff  = f_rise / f_decay;
      const double norm_coef = 1.0 / ( pow ( tau_diff, ( tau_prime / f_decay ) ) - pow ( tau_diff, ( tau_prime / f_rise ) ) );
      const int solver_id1 = c -> ptr_post [ d_post_i ] + local_idx [ d_post_i ];
      c -> post_c [ solver_id1 ] = d_post_c;
      c -> weight [ solver_id1 ] = norm_coef * f_weight;
      c -> erev   [ solver_id1 ] = f_erev;
      c -> decay  [ solver_id1 ] = exp ( - DT / f_decay );
      local_idx [ d_post_i ]++;
      c -> delay [ idx ] = d_delay;
      c -> id    [ idx ] = solver_id1;
      idx++;
      const int solver_id2 = c -> ptr_post [ d_post_i ] + local_idx [ d_post_i ];
      c -> post_c [ solver_id2 ] = d_post_c;
      c -> weight [ solver_id2 ] = - norm_coef * f_weight;
      c -> erev   [ solver_id2 ] = f_erev;
      c -> decay  [ solver_id2 ] = exp ( - DT / f_rise );
      local_idx [ d_post_i ]++;
      c -> delay [ idx ] = d_delay;
      c -> id    [ idx ] = solver_id2;
      idx++;
    }
    free ( local_idx );
    fclose ( file );
  }
  
  return c;
}

void finalize_connection ( conn_t *c )
{
  if ( c -> pre_table != NULL ) { free ( c -> pre_table ); }
  if ( c -> ptr_pre   != NULL ) { free ( c -> ptr_pre   ); }
  if ( c -> ptr_post  != NULL ) { free ( c -> ptr_post  ); }
  if ( c -> post_c != NULL ) { free ( c -> post_c ); }
  if ( c -> weight != NULL ) { free ( c -> weight ); }
  if ( c -> erev   != NULL ) { free ( c -> erev   ); }
  if ( c -> decay  != NULL ) { free ( c -> decay  ); }
  if ( c -> delay != NULL ) { free ( c -> delay ); }
  if ( c -> id    != NULL ) { free ( c -> id    ); }
  free ( c );
}
