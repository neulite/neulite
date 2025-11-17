// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "popl.h"
#include "hines.h"

hines_matrix_t *hines_matrix_allocate ( int n )
{
  hines_matrix_t *H = ( hines_matrix_t * ) malloc ( sizeof ( hines_matrix_t ) );
  
  H -> n_comp    = n;
  H -> Ad        = calloc ( H -> n_comp, sizeof ( double ) );
  H -> Api       = calloc ( H -> n_comp, sizeof ( double ) );
  H -> bu_Ad     = calloc ( H -> n_comp, sizeof ( double ) );
  H -> bu_Api    = calloc ( H -> n_comp, sizeof ( double ) );
  H -> parent_id = calloc ( H -> n_comp, sizeof ( int ) );
  
  return H;
}

hines_matrix_t *hines_matrix_initialize ( const population_t *u, const int pid )
{
  hines_matrix_t *H = hines_matrix_allocate ( u -> n_comp [ pid ] );

  for ( int i = 0; i < H -> n_comp; i++ ) {
    H -> Ad        [ i ] = 0.0;
    H -> Api       [ i ] = 0.0;
    H -> bu_Ad     [ i ] = 0.0;
    H -> bu_Api    [ i ] = 0.0;
    H -> parent_id [ i ] = ( int ) u -> parent [ u -> cid [ pid ] + i ];
  }
  return H;
}

void hines_matrix_finalize ( hines_matrix_t * H )
{
    free ( H -> Ad );
    free ( H -> Api );
    free ( H -> bu_Ad );
    free ( H -> bu_Api );
    free ( H -> parent_id );
    H -> n_comp = 0;
}
