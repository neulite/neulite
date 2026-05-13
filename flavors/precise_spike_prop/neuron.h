// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

#include "popl.h"

typedef struct {
  double *v, *ca, *i_ext; // size == # neurons * # compartments
  int *sid, *pid;         // size == # neurons; sid == soma id, pid == population id
  int n_neuron;
  void *sim; // used for single neuron simulation
} neuron_t;

extern neuron_t *initialize_neuron ( const population_t * );
extern void finalize_neuron ( neuron_t * );
