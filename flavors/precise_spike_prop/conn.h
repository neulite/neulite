// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

#include "popl.h"
#include "neuron.h"

typedef struct {
  int *post_c; // for solver
  double *weight, *erev, *decay; // for solver
  int *delay, *id; // for synapse
  int n_pre, n_post, n_conn;
  int *pre_table;
  int *ptr_pre, *ptr_post; // cumulative connection id
} conn_t;

extern conn_t *initialize_connection ( const population_t *, const neuron_t *, const char * );
extern void finalize_connection ( conn_t * );
