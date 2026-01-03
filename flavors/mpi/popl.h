// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

typedef enum { _DUMMY_, SOMA, AXON, APICAL, DEND, N_COMPTYPE } comptype_t;

typedef struct {

  // Morphology
  double *rad, *len, *area; // size == # populations * # compartments
  int *parent, *type; // size == # populations * # compartments

  // Passive
  double *cm, *ra, *gl, *vl; // size == # populations * # compartments

  // Size
  int *n_neuron, *n_comp, *cid; // size == # populations; cid = compartment id
  int n_popl;

  // Conductances, Ca2+ params (gamma, decay)
  double *gbar, *gamma, *decay; // size == # populations * N_GBAR (gbar, perisomatic); size == # populations * N_COMPTYPE (gamma, decay)

  // Local population id
  int start_pid, end_pid;
} population_t;

extern population_t *initialize_population ( const int, const int, const char * );
extern void finalize_population ( population_t * );
