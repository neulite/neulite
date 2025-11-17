// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include "ion.h"
#include "ion_func.h"
#include "popl.h"
#include "neuron.h"
#include "config.h"

ion_t *initialize_ion ( const neuron_t *n )
{
  ion_t *i = calloc ( 1, sizeof ( ion_t ) );

  i -> n_neuron = n -> n_neuron;

  i -> gate = calloc ( N_GATEVAL * i -> n_neuron, sizeof ( double ) );

  for ( int li = 0; li < i -> n_neuron; li++ ) {
    const double _v  = n ->  v [ n -> sid [ li ] + 0 ]; // compartment id 0 == SOMA
    const double _ca = n -> ca [ n -> sid [ li ] + 0 ]; // compartment id 0 == SOMA
    double *ion = &i -> gate [ N_GATEVAL * li ];
    ion [ M_NATS ] = inf_m_NaTs  ( _v );
    ion [ H_NATS ] = inf_h_NaTs  ( _v );
    ion [ M_NATA ] = inf_m_NaTa  ( _v );
    ion [ H_NATA ] = inf_h_NaTa  ( _v );
    ion [ H_NAP  ] = inf_h_Nap   ( _v );
    ion [ M_KV2  ] = inf_m_Kv2   ( _v );
    ion [ H1_KV2 ] = inf_h_Kv2   ( _v );
    ion [ H2_KV2 ] = inf_h_Kv2   ( _v );
    ion [ M_KV3  ] = inf_m_Kv3   ( _v );
    ion [ M_KP   ] = inf_m_KP    ( _v );
    ion [ H_KP   ] = inf_h_KP    ( _v );
    ion [ M_KT   ] = inf_m_KT    ( _v );
    ion [ H_KT   ] = inf_h_KT    ( _v );
    ion [ M_KD   ] = inf_m_Kd    ( _v );
    ion [ H_KD   ] = inf_h_Kd    ( _v );
    ion [ M_IM   ] = inf_m_Im    ( _v );
    ion [ M_IMV2 ] = inf_m_Imv2  ( _v );
    ion [ M_IH   ] = inf_m_Ih    ( _v );
    ion [ Z_SK   ] = inf_z_SK    ( _v, _ca );
    ion [ M_CAHVA ] = inf_m_CaHVA ( _v );
    ion [ H_CAHVA ] = inf_h_CaHVA ( _v );
    ion [ M_CALVA ] = inf_m_CaLVA ( _v );
    ion [ H_CALVA ] = inf_h_CaLVA ( _v );

    { // NaV
      double vca [ N_STATE_NAV ] [ N_STATE_NAV ] = {};
      double vcb [ N_STATE_NAV ] = {};
      Init_Nav_param ( _v, vca );
      for ( int col = 0; col < N_STATE_NAV; col++ ) { vca [ N_STATE_NAV - 1 ] [ col ] = 1.0; } // I6 channel 
      vcb [ N_STATE_NAV - 1 ] = 1.0;
      gaussian_elimination ( N_STATE_NAV, vca, vcb );
      ion [ C1_NaV ] = vcb [ 0  ];
      ion [ C2_NaV ] = vcb [ 1  ];
      ion [ C3_NaV ] = vcb [ 2  ];
      ion [ C4_NaV ] = vcb [ 3  ];
      ion [ C5_NaV ] = vcb [ 4  ];
      ion [ I1_NaV ] = vcb [ 5  ];
      ion [ I2_NaV ] = vcb [ 6  ];
      ion [ I3_NaV ] = vcb [ 7  ];
      ion [ I4_NaV ] = vcb [ 8  ];
      ion [ I5_NaV ] = vcb [ 9  ];
      ion [ I6_NaV ] = vcb [ 10 ];
      ion [ OO_NaV ] = vcb [ 11 ];
    }
  }
  
  return i;
}

void finalize_ion ( ion_t *i )
{
  if ( i -> gate != NULL ) { free ( i -> gate ); }
  free ( i );
}

void update_ion ( const int id, const neuron_t * __restrict__ n, ion_t * __restrict__ i, const double dt )
{
  const double _v  = n -> v  [ n -> sid [ id ] + 0 ];
  const double _ca = n -> ca [ n -> sid [ id ] + 0 ];
  double *ion = &i -> gate [ N_GATEVAL * id ];

  Nav_update ( _v, &ion [ OO_NaV ], &ion [ C1_NaV ], &ion [ C2_NaV ], &ion [ C3_NaV ], &ion [ C4_NaV ], &ion [ C5_NaV ],
	           &ion [ I1_NaV ], &ion [ I2_NaV ], &ion [ I3_NaV ], &ion [ I4_NaV ], &ion [ I5_NaV ], &ion [ I6_NaV ] );

  ion [ M_NATS ]  = inf_m_NaTs ( _v )      + ( ion [ M_NATS ]  - inf_m_NaTs ( _v ) )    * exp ( - dt / tau_m_NaTs   ( _v ) );
  ion [ H_NATS ]  = inf_h_NaTs ( _v )      + ( ion [ H_NATS ]  - inf_h_NaTs ( _v ) )    * exp ( - dt / tau_h_NaTs   ( _v ) );
  ion [ M_NATA ]  = inf_m_NaTa ( _v )      + ( ion [ M_NATA ]  - inf_m_NaTa ( _v ) )    * exp ( - dt / tau_m_NaTa   ( _v ) );
  ion [ H_NATA ]  = inf_h_NaTa ( _v )      + ( ion [ H_NATA ]  - inf_h_NaTa ( _v ) )    * exp ( - dt / tau_h_NaTa   ( _v ) );
  ion [ H_NAP  ]  = inf_h_Nap  ( _v )      + ( ion [ H_NAP  ]  - inf_h_Nap  ( _v ) )    * exp ( - dt / tau_h_Nap    ( _v ) );
  ion [ M_KV2  ]  = inf_m_Kv2  ( _v )      + ( ion [ M_KV2  ]  - inf_m_Kv2  ( _v ) )    * exp ( - dt / tau_m_Kv2    ( _v ) );
  ion [ H1_KV2 ]  = inf_h_Kv2  ( _v )      + ( ion [ H1_KV2 ]  - inf_h_Kv2  ( _v ) )    * exp ( - dt / tau_h1_Kv2   ( _v ) );
  ion [ H2_KV2 ]  = inf_h_Kv2  ( _v )      + ( ion [ H2_KV2 ]  - inf_h_Kv2  ( _v ) )    * exp ( - dt / tau_h2_Kv2   ( _v ) );
  ion [ M_KV3  ]  = inf_m_Kv3  ( _v )      + ( ion [ M_KV3  ]  - inf_m_Kv3  ( _v ) )    * exp ( - dt / tau_m_Kv3    ( _v ) );
  ion [ M_KP   ]  = inf_m_KP   ( _v )      + ( ion [ M_KP   ]  - inf_m_KP   ( _v ) )    * exp ( - dt / tau_m_KP     ( _v ) );
  ion [ H_KP   ]  = inf_h_KP   ( _v )      + ( ion [ H_KP   ]  - inf_h_KP   ( _v ) )    * exp ( - dt / tau_h_KP     ( _v ) );
  ion [ M_KT   ]  = inf_m_KT   ( _v )      + ( ion [ M_KT   ]  - inf_m_KT   ( _v ) )    * exp ( - dt / tau_m_KT     ( _v ) );
  ion [ H_KT   ]  = inf_h_KT   ( _v )      + ( ion [ H_KT   ]  - inf_h_KT   ( _v ) )    * exp ( - dt / tau_h_KT     ( _v ) );
  ion [ M_KD   ]  = inf_m_Kd   ( _v )      + ( ion [ M_KD   ]  - inf_m_Kd   ( _v ) )    * exp ( - dt / tau_m_Kd     ( _v ) );
  ion [ H_KD   ]  = inf_h_Kd   ( _v )      + ( ion [ H_KD   ]  - inf_h_Kd   ( _v ) )    * exp ( - dt / tau_h_Kd     ( _v ) );
  ion [ M_IM   ]  = inf_m_Im   ( _v )      + ( ion [ M_IM   ]  - inf_m_Im   ( _v ) )    * exp ( - dt / tau_m_Im     ( _v ) );
  ion [ M_IMV2 ]  = inf_m_Imv2 ( _v )      + ( ion [ M_IMV2 ]  - inf_m_Imv2 ( _v ) )    * exp ( - dt / tau_m_Imv2   ( _v ) );
  ion [ M_IH   ]  = inf_m_Ih   ( _v )      + ( ion [ M_IH   ]  - inf_m_Ih   ( _v ) )    * exp ( - dt / tau_m_Ih     ( _v ) );
  ion [ Z_SK   ]  = inf_z_SK   ( _v, _ca ) + ( ion [ Z_SK   ]  - inf_z_SK ( _v, _ca ) ) * exp ( - dt / tau_z_SK     ( _v ) );
  ion [ M_CAHVA ] = inf_m_CaHVA ( _v )     + ( ion [ M_CAHVA ] - inf_m_CaHVA ( _v ) )   * exp ( - dt / tau_m_CaHVA  ( _v ) );
  ion [ H_CAHVA ] = inf_h_CaHVA ( _v )     + ( ion [ H_CAHVA ] - inf_h_CaHVA ( _v ) )   * exp ( - dt / tau_h_CaHVA  ( _v ) );
  ion [ M_CALVA ] = inf_m_CaLVA ( _v )     + ( ion [ M_CALVA ] - inf_m_CaLVA ( _v ) )   * exp ( - dt / tau_m_CaLVA  ( _v ) );
  ion [ H_CALVA ] = inf_h_CaLVA ( _v )     + ( ion [ H_CALVA ] - inf_h_CaLVA ( _v ) )   * exp ( - dt / tau_h_CaLVA  ( _v ) );
}

void update_ca ( const int id, const population_t * __restrict__ u, const ion_t * __restrict__ i, neuron_t * __restrict__ n, const double dt )
{
  const int sid = n -> sid [ id ];
  const int pid = n -> pid [ id ];
  const double *gbar = &u -> gbar [ N_GBAR * pid ]; // perisomatic
  const double *ion  = &i -> gate [ N_GATEVAL * id ];

  {
    const double v  = n -> v  [ sid ];
    const double ca = n -> ca [ sid ];
    const double area  = u -> area [ u -> cid [ pid ] ];
    const double gamma = u -> gamma [ SOMA + N_COMPTYPE * pid ]; // perisomatic
    const double decay = u -> decay [ SOMA + N_COMPTYPE * pid ]; // perisomatic
    const double i_ca =  (1e-3 * ( v - rev_ca ( ca ) ) * ( gbar [ G_CAHVA ] * ion [ M_CAHVA ] * ion [ M_CAHVA ] * ion [ H_CAHVA ]
							   + gbar [ G_CALVA ] * ion [ M_CALVA ] * ion [ M_CALVA ] * ion [ H_CALVA ] ) ) / area;
    n -> ca [ sid ] += dt * dcadt ( ca, i_ca, gamma, decay );
  }
}

void calc_lhs_and_rhs ( const population_t * __restrict__ u, const neuron_t * __restrict__ n, const ion_t * __restrict__ i, const int pid, const int id, double * __restrict__ lhs, double * __restrict__ rhs )
{
  const double _v  = n -> v  [ n -> sid [ id ] + 0 ];
  const double _ca = n -> ca [ n -> sid [ id ] + 0 ];
  double *ion  = &i -> gate [ N_GATEVAL * id ];
  double *gbar = &u -> gbar [ N_GBAR * pid ]; // perisomatic
  double _l = 0.0, _r = 0.0;
  double \
  _c = gbar [ G_NAV   ] * ion [ OO_NaV ];                                                                  _l += _c; _r += _c * V_NA;
  _c = gbar [ G_NATS  ] * ion [ M_NATS ] * ion [ M_NATS ] * ion [ M_NATS ] * ion [ H_NATS ];               _l += _c; _r += _c * V_NA;
  _c = gbar [ G_NATA  ] * ion [ M_NATA ] * ion [ M_NATA ] * ion [ M_NATA ] * ion [ H_NATA ];               _l += _c; _r += _c * V_NA;
  _c = gbar [ G_NAP   ] * inf_m_Nap ( _v ) * ion [ H_NAP ];                                                _l += _c; _r += _c * V_NA;
  _c = gbar [ G_KV2   ] * ion [ M_KV2 ] * ion [ M_KV2 ] * ( 0.5 * ion [ H1_KV2 ] + 0.5 * ion [ H2_KV2 ] ); _l += _c; _r += _c * V_K;
  _c = gbar [ G_KV3   ] * ion [ M_KV3 ];                                                                   _l += _c; _r += _c * V_K;
  _c = gbar [ G_KP    ] * ion [ M_KP ] * ion [ M_KP ] * ion [ H_KP ];                                      _l += _c; _r += _c * V_K;
  _c = gbar [ G_KT    ] * ion [ M_KT ] * ion [ M_KT ] * ion [ M_KT ] * ion [ M_KT ] * ion [ H_KT ];        _l += _c; _r += _c * V_K;
  _c = gbar [ G_KD    ] * ion [ M_KD ] * ion [ H_KD ];                                                     _l += _c; _r += _c * V_K;
  _c = gbar [ G_IM    ] * ion [ M_IM ];                                                                    _l += _c; _r += _c * V_K;
  _c = gbar [ G_IMV2  ] * ion [ M_IMV2 ];                                                                  _l += _c; _r += _c * V_K;
  _c = gbar [ G_IH    ] * ion [ M_IH ];                                                                    _l += _c; _r += _c * V_HCN;
  _c = gbar [ G_SK    ] * ion [ Z_SK ];                                                                    _l += _c; _r += _c * V_K;
  _c = gbar [ G_CAHVA ] * ion [ M_CAHVA ] * ion [ M_CAHVA ] * ion [ H_CAHVA ];                             _l += _c; _r += _c * rev_ca ( _ca );
  _c = gbar [ G_CALVA ] * ion [ M_CALVA ] * ion [ M_CALVA ] * ion [ H_CALVA ];                             _l += _c; _r += _c * rev_ca ( _ca );
  *lhs = _l;
  *rhs = _r;
}
