// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

#include "popl.h"
#include "neuron.h"
#include "ion.h"
#include "conn.h"
#include "synapse.h"
#include "solver.h"

typedef struct {
  population_t *u;
  neuron_t     *n;
  ion_t        *i;
  conn_t       *c;
  synapse_t    *s;
  FILE *v_dat, *s_dat;
  double *v_prev;
  int *spike;
  int tick, inv_dt;
} network_t;

extern network_t *initialize_network ( const char *, const char * );
extern void finalize_network ( network_t * );
extern void output_v ( const double, network_t * );
extern void set_current ( const double, network_t *, double ( *current ) ( const double, const int ) );
extern void solve_network ( network_t *, solver_t * );
extern void spike_detection ( network_t * );
extern void spike_propagation ( const double, network_t * );
