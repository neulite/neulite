// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

#include "neuron.h"
#include "ion.h"
#include "conn.h"
#include "synapse.h"
#include "hines.h"

typedef struct {
  hines_matrix_t *H;
  double *b;
} linsys_t;

typedef struct {
  linsys_t *linsys;
  int n_popl, n_neuron;
} solver_t;

extern solver_t *initialize_solver ( const population_t * );
extern void solve ( const population_t *, neuron_t *, ion_t *, const conn_t *, synapse_t *, solver_t *solver );
extern void finalize_solver ( solver_t * );
