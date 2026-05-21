// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025,2026 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

typedef struct {
  double *sum0;
  int *delay;
} synapse_t;

extern synapse_t *initialize_synapse ( conn_t * );
extern void update_synapse ( const int, const conn_t *, synapse_t * );
extern void finalize_synapse ( synapse_t * );
