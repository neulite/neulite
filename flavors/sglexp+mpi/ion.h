// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

#include <stdio.h>
#include "popl.h"
#include "neuron.h"

#define V_NA    (  53.0  )
#define V_K     ( -107.0 )
#define V_HCN   ( -45.0  )

#define N_STATE_NAV 12
typedef enum { G_NAV, G_NATS, G_NATA, G_NAP, G_KV2, G_KV3, G_KP, G_KT, G_KD, G_IM, G_IMV2, G_IH, G_SK, G_CAHVA, G_CALVA, N_GBAR } ion_gbar_t;
typedef enum { M_NATS, H_NATS, M_NATA, H_NATA, H_NAP, M_KV2, H1_KV2, H2_KV2, M_KV3, M_KP, H_KP, M_KT, H_KT, M_KD, H_KD, M_IM, M_IMV2, M_IH, Z_SK, M_CAHVA, H_CAHVA, M_CALVA, H_CALVA, OO_NaV, C1_NaV, C2_NaV, C3_NaV, C4_NaV, C5_NaV, I1_NaV, I2_NaV, I3_NaV, I4_NaV, I5_NaV, I6_NaV, N_GATEVAL } ion_gateval_t;

typedef struct {
  double *gate; // size == # neurons * N_GATEVAL // perisomatic
  int n_neuron;
} ion_t;

extern ion_t *initialize_ion ( const neuron_t * );
extern void finalize_ion ( ion_t * );
extern void update_ion ( const int, const neuron_t *, ion_t *, const double );
extern void update_ca ( const int, const population_t *, const ion_t *, neuron_t *, const double );
extern void calc_lhs_and_rhs ( const population_t *, const neuron_t *, const ion_t *, const int, const int, double *, double * );
