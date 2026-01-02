// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

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
  int *spike;
} network_t;

extern network_t *initialize_network ( const char *, const char * );
extern void finalize_network ( network_t * );
extern void set_current ( const int, network_t *, double ( *current ) ( const int, const int ) );
extern void solve_network ( const int, network_t *, solver_t * );
extern void spike_propagation ( const int, network_t * );
