// When config.h does not exist at compile, this file is copied to config.h and used.

#pragma once

#undef DEBUG

// Simulation parameters
#define TSTOP ( 2000.0 )
#define DT ( 0.1 )
#define INV_DT ( ( int ) ( 1.0 / ( DT ) ) )

// Neuron parameters
#define SPIKE_THRESHOLD ( -15.0 )
#define ALLACTIVE ( 0 ) // Set to 1 for allactive models

// Current injection parameters
#define I_AMP ( 0.12 )
#define I_DELAY ( 500.0 )
#define I_DURATION ( 1000.0 )
