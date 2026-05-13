// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#pragma once

typedef struct {
    int n_comp;
    double *Ad, *Api;
    double *bu_Ad, *bu_Api;
    int *parent_id;
} hines_matrix_t;
