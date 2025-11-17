// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "popl.h"
#include "ion.h"
#include "config.h"

#define MAX_N_COMP ( 16384 )

extern int strip_comment_destructive ( char * );
extern int remove_blank_destructive_for_csv ( char * );

typedef struct {
  int id, type, parent;
  double x, y, z, r;
} swc_t;

typedef struct {
  int parent;
  swc_t proximal;
  swc_t distal;
  int type;
} segment_t;

typedef struct {
  segment_t *data;
  int size;
} segment_tree_t;

static segment_tree_t *read_in_neuron_style ( const char *filename, const int n_comp )
{
  FILE *file = fopen ( filename, "r" );
  if ( ! file ) { fprintf ( stderr, "Error: no such SWC file %s\n", filename ); exit ( 1 ); }

  swc_t swc [ n_comp ];
  int n_child [ n_comp ]; for ( int i = 0; i < n_comp; i++ ) { n_child [ i ] = 0; }

  char buf [ 1024 ];
  int i = 0;
  while ( fgets (buf, 1024, file ) ) {
    if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
    int id, d_type, d_parent;
    double f_r, f_x, f_y, f_z;
    const int nf = sscanf ( buf, "%d %d %lf %lf %lf %lf %d", &id, &d_type, &f_x, &f_y, &f_z, &f_r, &d_parent );
    assert ( nf == 7 );
    swc_t seg = { .id = id, .type = d_type, .x = f_x, .y = f_y, .z = f_z, .r = f_r, .parent = d_parent };
    swc [ id ] = seg;
    if ( d_parent >= 0 ) { n_child [ d_parent ]++; }
    assert ( id == i );
    i++;
  }
  fclose ( file );
  
  segment_tree_t *segment_tree = calloc ( 1, sizeof ( segment_tree_t ) );

  segment_tree -> data = calloc ( n_comp, sizeof ( segment_t ) );
  segment_tree -> size = 0;
  int segid_from_rid [ n_comp ]; for ( int i = 0; i < n_comp; i++ ) { segid_from_rid [ i ] = 0; }
  
  if ( swc [ 0 ] . type == SOMA ) {
    swc_t distal = swc [ 0 ];
    swc_t proximal = { .id = -1, .type = distal.type, .x = distal.x - 2*distal.r, .y = distal.y, .z = distal.z, .r = distal.r, .parent = distal.parent };
    segment_t seg = { .parent = -1, .proximal = proximal, .distal = distal, .type = SOMA };
    segment_tree -> data [ segment_tree -> size ] = seg;
    segid_from_rid [ 0 ] = segment_tree -> size;
    segment_tree -> size++;
  } else {
    fprintf ( stderr, "Error: swc [ 0 ] != SOMA.\n" );
    //print_swc ( &swc [ 0 ] );
    exit ( 1 );
  }

  for ( int rid = 1; rid < n_comp; rid++ ) {
    swc_t distal   = swc [ rid ];
    swc_t proximal = swc [ distal.parent ];
    const int grandparent_type = ( proximal.parent == -1 ) ? _DUMMY_ : swc [ proximal.parent ].type;

    // 1st block:
    //     The soma-dend0 segment is skipped if dend0 has one or more child records.
    // 2nd block:
    //     To match n_segment to n_record, a dendrite segment connected to the soma is split. 
    //     swc:(soma-)dend0-dend1 → seg:dend0-mid + mid-dend1
    // 3rd block:
    //     All other cases.
    //     - A soma-axon segment is created as a cylinder. 

    // In the case where dend0 (from soma-dend0) has multiple child records,
    // the condition (segid_from_rid[dend0.rid] == segid_from_rid[soma.rid]) detects
    // that the tree has not yet created a segment for dend0.rid.
    // Therefore, the 2nd block creates dend0-mid and sets segid_from_rid[dend0.rid] to dend0-mid.id.
    // dend0-dend2 is not split, and is directly connected to the soma.

    if ( proximal.type == SOMA && n_child [ rid ] != 0 && distal.type != AXON ) {
      segid_from_rid [ rid ] = segid_from_rid [ proximal.id ];
    } else if ( grandparent_type == SOMA && distal.type != AXON && segid_from_rid [ proximal.id ] == segid_from_rid [ proximal.parent ] ) {
      int pid = segid_from_rid [ proximal.id ];
      swc_t mid = {
	.id = -1,
	.type = distal.type,
	.x = ( distal.x + proximal.x ) * 0.5,
	.y = ( distal.y + proximal.y ) * 0.5,
	.z = ( distal.z + proximal.z ) * 0.5,
	.r = ( distal.r + proximal.r ) * 0.5,
	.parent = -1
      };

      segment_t first_half = {
	.parent = pid,
	.proximal = proximal,
	.distal = mid,
	.type = mid.type
      };
      segment_tree -> data [ segment_tree -> size ] = first_half;
      segid_from_rid [ proximal.id ] = segment_tree -> size;
      segment_tree -> size++;
      
      segment_t second_half = {
	.parent = pid,
	.proximal = mid,
	.distal = distal,
	.type = distal.type
      };
      segment_tree -> data [ segment_tree -> size ] = second_half;
      segid_from_rid [ distal.id ] = segment_tree -> size;
      segment_tree -> size++;
    } else {
      const int has_skipped_parent = ( grandparent_type == SOMA && distal.type != AXON );
      const int pid = ( has_skipped_parent ) ? segid_from_rid [ proximal.parent ] : segid_from_rid [ proximal.id ];

      if ( proximal.type == SOMA && distal.type == AXON ) {
	proximal.r = distal.r;
      }
      segment_t seg = {
	.parent = pid,
	.proximal = proximal,
	.distal = distal,
	.type = distal.type
      };
      segment_tree -> data [ segment_tree -> size ] = seg;
      segid_from_rid [ rid ] = segment_tree -> size;
      segment_tree -> size++;
    }
  }

  return segment_tree;
}

static void get_population_size ( const char *filename, const int n_popl, int n_neuron [ n_popl ], int n_comp [ n_popl ] )
{
  FILE *file = fopen ( filename, "r" );
  if ( ! file ) { fprintf ( stderr, "Error: no such file %s\n", filename ); exit ( 1 ); }
  int i = 0;
  char buf [ 1024 ];
  while ( fgets ( buf, 1024, file ) ) {
    if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
    if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
    int dn_neuron, dn_comp;
    const int nf = sscanf ( buf, "%d,%d", &dn_neuron, &dn_comp );
    assert ( nf == 2 );
    n_neuron [ i ] = dn_neuron;
    n_comp [ i ] = dn_comp;
    i++;
  }
  fclose ( file );
}

static void read_swc_file ( population_t *p, const int pid, const char *filename )
{
  segment_tree_t *st = read_in_neuron_style ( filename, p -> n_comp [ pid ] );

  int n_comp = st -> size;
  p -> n_comp [ pid ] = n_comp; // update n_comp
  
  const int offset = p -> cid [ pid ];
  double *rad  = &p -> rad  [ offset ];
  double *len  = &p -> len  [ offset ];
  double *area = &p -> area [ offset ];
  int *parent = &p -> parent [ offset ];
  int *type   = &p -> type   [ offset ];

  for(int i = 0; i < n_comp; i++){
    segment_t *seg = &st -> data [ i ];
    parent [ i ] = seg -> parent;
    type   [ i ] = seg -> type;
    const double dx = seg -> distal.x - seg -> proximal.x;
    const double dy = seg -> distal.y - seg -> proximal.y;
    const double dz = seg -> distal.z - seg -> proximal.z;
    const double dr = seg -> distal.r - seg -> proximal.r;
    const double _len = sqrt ( dx * dx + dy * dy + dz * dz ) * 1.0e-4; // [mum -> cm]
    area [ i ] = M_PI * ( seg -> proximal.r + seg -> distal.r ) * 1.0e-4 * sqrt ( dr * dr * 1.0e-8 + _len * _len );
    len  [ i ] = _len;
    rad  [ i ] = seg -> distal.r * 1.0e-4; // TODO
  }
  
  free ( st -> data );
  st -> size = 0;
  free ( st );
}

static void read_pas_file ( population_t *u, const int pid, const char *filename )
{
  double _cm [ N_COMPTYPE ], _ra [ N_COMPTYPE ], _gl [ N_COMPTYPE ], _vl [ N_COMPTYPE ];
  {
    FILE *file = fopen ( filename, "r" );
    if ( ! file ) { fprintf ( stderr, "Error: no such PAS file %s\n", filename ); exit ( 1 ); }
    char buf [ 1024 ];
    while ( fgets ( buf, 1024, file ) ) {
      if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
      if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
      
      int d_type;
      double f_cm, f_ra, f_gl, f_vl;
      const int nf = sscanf ( buf, "%d,%lf,%lf,%lf,%lf", &d_type, &f_cm, &f_ra, &f_gl, &f_vl );
      assert ( nf == 5 );
      _cm [ d_type ] = f_cm;
      _ra [ d_type ] = f_ra;
      _gl [ d_type ] = f_gl;
      _vl [ d_type ] = f_vl;
    }
    fclose ( file );
  }

  const int offset = u -> cid [ pid ];
  double *cm = &u -> cm [ offset ];
  double *ra = &u -> ra [ offset ];
  double *gl = &u -> gl [ offset ];
  double *vl = &u -> vl [ offset ];

  for ( int i = 0; i < u -> n_comp [ pid ]; i++ ) {
    const int type    = u -> type   [ offset + i ];
    const double area = u -> area   [ offset + i ];
    cm [ i ] = _cm [ type ] * area; // [muF]
    ra [ i ] = _ra [ type ] * 1e-3;
    gl [ i ] = _gl [ type ] * area * 1e3; /* CONVERSION: 1e3 from S to mS */
    vl [ i ] = _vl [ type ];
  }
}

static void read_ion_file ( population_t *u, const int pid, const char *filename )
{
  FILE *file = fopen ( filename, "r" );
  if ( ! file ) { fprintf ( stderr, "Error: no such ION file %s\n", filename ); exit ( 1 ); }
  char buf [ 1024 ];
  while ( fgets ( buf, 1024, file ) ) {
    if ( strip_comment_destructive ( buf ) == 0 ) { continue; }
    if ( remove_blank_destructive_for_csv ( buf ) == 0 ) { continue; }
      
    int d_type;
    double f_cm, f_ra, f_gl, f_vl, f_gamma, f_decay;
    double f [ N_GBAR ] = { 0.0 };
    const int nf = sscanf ( buf, "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
			    &d_type, &f_cm, &f_ra, &f_gl, &f_vl, &f_gamma, &f_decay, &f[0], &f[1], &f[2], &f[3], &f[4], &f[5], &f[6], &f[7], &f[8], &f[9], &f[10], &f[11], &f[12], &f[13], &f[14] );
    assert ( nf == 5 || nf == 22 );
    if ( nf == 22 ) {
      u -> gamma [ d_type + N_COMPTYPE * pid ] = f_gamma;
      u -> decay [ d_type + N_COMPTYPE * pid ] = f_decay;
    }
    if ( nf == 22 && ALLACTIVE == 1 ) {
      const int n_comp = u -> n_comp [ pid ];
      for ( int i = 0; i < n_comp; i++ ) {
	const double area = u -> area [ u -> cid [ pid ] + i ];
	for ( int j = 0; j < N_GBAR; j++ ) {
	  u -> gbar [ j + N_GBAR * ( i + n_comp * pid ) ] = f [ j ] * area * 1e3; /* CONVERSION: 1e3 from S to mS */
	}
      }
    } else { // Default is perisomatic
      if ( d_type == SOMA ) {
	const double area = u -> area [ u -> cid [ pid ] + 0 ]; // 0 == SOMA
	for ( int i = 0; i < N_GBAR; i++ ) { u -> gbar [ i + N_GBAR * pid ] = f[i] * area * 1e3; /* CONVERSION: 1e3 from S to mS */ }
      }
    }
  }
  fclose ( file );
}
