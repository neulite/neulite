// SPDX-License-Identifier: GPL-2.0-only
// Copyright (C) 2024,2025 Neulite Core Team <neulite-core@numericalbrain.org>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"

#define MIN(a,b) ( ( ( a ) < ( b ) ) ? ( a ) : ( b ) )
#define MAX(a,b) ( ( ( a ) > ( b ) ) ? ( a ) : ( b ) )

#define Ca1OUT      ( 2.0 )  // [mM]
#define F           ( 9.6485e4 ) // Faraday constant [s*A/mol]
#define RVAL        ( 8.31446261815324 ) // [J/K*mol]
#define MIN_CA      ( 1e-4 ) // [mM]
#define DEPTH       ( 0.1e-4 ) // [ 0.1mu m -> 0.1e-4cm]

double ca_init ( void ) { return MIN_CA; }
double rev_ca ( const double ca ) { const double celcius = 34.0; return 1000 * (( RVAL * (273.0 + celcius )) / (2 * F)) * log ( Ca1OUT / ca ); }
double dcadt ( const double ca, const double i_ca, const double gamma, const double decay )
{
  return ( - (10000) * ( i_ca * gamma / ( 2 * F * DEPTH * 1e4 ) ) - ( ca - MIN_CA ) / decay );
}

static double vtrap (const double x, const double y){
  if ( fabs ( x / y ) < 1e-6 ) {
    return y * ( 1.0 - ( x / y ) / 2.0 );
  }else{
    return x / ( exp ( x / y ) - 1.0 );
  }
}

static double alpha_m_NaTs ( const double v ) {
  const double malphaF = 0.182;
  const double mvhalf  = -40.0; // (mV)
  const double mk      = 6.0; // (mV)
  return malphaF * vtrap ( - ( v - mvhalf ), mk ); }

static double beta_m_NaTs  ( const double v ) {
  const double mbetaF  = 0.124;
  const double mvhalf  = -40.0; // (mV)
  const double mk      = 6.0; // (mV)
  return mbetaF * vtrap ( ( v - mvhalf ), mk ); }

static double alpha_h_NaTs ( const double v ) {
  const double halphaF = 0.015;
  const double hvhalf = -66.0; // (mV)
  const double hk = 6.0; // (mV)
  return halphaF * vtrap ( ( v - hvhalf ), hk ); }

static double beta_h_NaTs  ( const double v ) {
  const double hbetaF = 0.015;
  const double hvhalf = -66.0; // (mV)
  const double hk = 6.0; // (mV)
  return hbetaF * vtrap ( - ( v - hvhalf ), hk ); }

static double inf_m_NaTs   ( const double v ) { return alpha_m_NaTs(v) / (alpha_m_NaTs(v) + beta_m_NaTs(v)); }
static double inf_h_NaTs   ( const double v ) { return alpha_h_NaTs(v) / (alpha_h_NaTs(v) + beta_h_NaTs(v)); }
static double tau_m_NaTs   ( const double v ) {
  const double celsius = 34.0;
  const double qt = pow(2.3, (celsius-23.0)/10.0);
  return (1.0 / (alpha_m_NaTs(v) + beta_m_NaTs(v))) / qt; }

static double tau_h_NaTs   ( const double v ) {
  const double celsius = 34.0;
  const double qt = pow(2.3, (celsius-23.0)/10.0);
  return (1.0 / (alpha_h_NaTs(v) + beta_h_NaTs(v))) / qt; }


static double alpha_m_NaTa ( const double v ){
  const double malphaF = 0.182;
  const double mvhalf = -48.0;
  const double mk = 6.0;
  return malphaF * vtrap(-(v - mvhalf), mk);
}

static double beta_m_NaTa ( const double v ){
  const double mbetaF= 0.124;
  const double mvhalf = -48.0;
  const double mk = 6.0;
  return mbetaF * vtrap((v - mvhalf), mk);
}

static double alpha_h_NaTa ( const double v ){
  const double halphaF = 0.015;
  const double hvhalf = -69.0;
  const double hk     = 6.0;
  return halphaF * vtrap(v - hvhalf, hk);
}

static double beta_h_NaTa ( const double v ){
  const double hbetaF = 0.015;
  const double hvhalf = -69.0;
  const double hk     = 6.0;
  return hbetaF * vtrap(-(v - hvhalf), hk);
}

static double inf_m_NaTa ( const double v ){
  return alpha_m_NaTa(v) / (alpha_m_NaTa(v) + beta_m_NaTa(v));
}

static double inf_h_NaTa ( const double v ){
  return alpha_h_NaTa(v) / (alpha_h_NaTa(v) + beta_h_NaTa(v));
}

static double tau_m_NaTa ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 23.0 ) / 10.0 );
  return (1.0 / (alpha_m_NaTa(v) + beta_m_NaTa(v))) / qt;
}

static double tau_h_NaTa ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 23.0 ) / 10.0 );
  return (1.0 / (alpha_h_NaTa(v) + beta_h_NaTa(v)))/ qt;
}


static double alpha_h_Nap ( const double v ){
  return 2.88e-6 * vtrap(v + 17.0, 4.63);
  //return 2.88e-3 * vtrap(v + 17.0, 4.63);
}

static double beta_h_Nap ( const double v ){
  return 6.94e-6 * vtrap(-(v + 64.4), 2.63);
  //return 6.94e-3 * vtrap(-(v + 64.4), 2.63);
}

extern double inf_m_Nap ( const double v ){
  return 1.0/(1.0+exp((v- (-52.6))/-4.6));
}

static double inf_h_Nap ( const double v ){
  return 1.0/(1.0+exp((v- (-48.8))/10.0));
}

static double tau_h_Nap ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3 , (celcius - 21.0)/10.0 );
  return (1.0 / (alpha_h_Nap(v) + beta_h_Nap(v))) / qt;
}

static double alpha_m_Kv2 ( const double v ){
  return 0.12 * vtrap( -(v - 43.0), 11.0);
}

static double beta_m_Kv2 ( const double v ){
  return 0.02 * exp(-(v + 1.27) / 120.0);
}

static double inf_m_Kv2 ( const double v ){
  return alpha_m_Kv2(v) / (alpha_m_Kv2(v) + beta_m_Kv2(v));
}

static double inf_h_Kv2 ( const double v ){
  return 1.0/(1.0 + exp((v + 58.0) / 11.0));
}

static double tau_m_Kv2 ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0) / 10.0);
  return 2.5 * (1.0 / (qt * (alpha_m_Kv2(v) + beta_m_Kv2(v))));
}

static double tau_h1_Kv2 ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0) / 10.0 );
  return (360 + (1010 + 23.7 * (v + 54)) * exp(-((v + 75) / 48) *((v + 75) / 48))) / qt;
}

static double tau_h2_Kv2 ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0) / 10.0 );
  return (2350 + 1380 * exp(-0.011 * v) - 210 * exp(-0.03 * v)) / qt;
}

static double inf_m_Kv3 ( const double v ){
  const double vshift = 0.0;
  return  1.0/(1.0+exp(((v -(18.700 + vshift))/(-9.700))));
}

static double tau_m_Kv3 ( const double v ){
  const double vshift = 0.0;
  return 0.2 * 20.000 / (1+exp(((v -(-46.560 + vshift))/(-44.140))));
}


static double inf_m_KP ( const double v ){
  const double vshift = 0.0;
  return 1.0 / (1.0 + exp(-(v - (-14.3 + vshift)) / 14.6));
}

static double inf_h_KP ( const double v ) {
  const double vshift = 0.0; //-1.1;
  return 1.0/(1.0 + exp(-(v - (-54.0 + vshift))/-11.0));
}

static double tau_m_KP ( const double v ){

  const double vshift = 0.0;
  const double tauF = 1.0;
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0)/ 10.0 ); 
    if (v < -50.0 + vshift){
      return tauF * (1.25+175.03*exp(-(v - vshift) * -0.026))/qt;
    } else {
      return tauF * (1.25+13*exp(-(v - vshift) * 0.026))/qt;
    }
}

static double tau_h_KP ( const double v ){
  const double vshift = 0.0;
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0 )/ 10.0 );
  return (360.0+(1010.0+24.0*(v - (-55.0 + vshift)))*exp(-((v - (-75.0 + vshift))/48)*((v - (-75.0 + vshift))/48.0)))/qt;
}

static double inf_m_KT ( const double v ){
  const double vshift = 0.0;
  return 1.0/(1.0 + exp(-(v - (-47 + vshift)) / 29.0));
}

static double inf_h_KT ( const double v ){
  const double vshift = 0.0;
  return 1.0/(1.0 + exp(-(v + 66.0 - vshift) / -10.0));
}

static double tau_m_KT ( const double v ){
  double mTauF = 1.0;
  double vshift = 0.0;
  double celcius = 34.0;
  double qt = pow ( 2.3, (celcius - 21.0)/10.0);
  return (0.34 + mTauF * 0.92 * exp(- ((v + 71.0 - vshift) / 59.0) * ((v + 71.0 - vshift) / 59.0))) / qt;
}

static double tau_h_KT ( const double v ){
  double hTauF = 1.0;
  double vshift = 0.0;
  double celcius = 34.0;
  double qt = pow ( 2.3, (celcius - 21.0) / 10.0 );
  return (8.0 + hTauF * 49.0 * exp(- ((v + 73.0 - vshift) / 23.0) * ((v + 73.0 - vshift) / 23.0 ))) / qt;
}

static double inf_m_Kd ( const double v ){
  return 1.0 - 1.0 / (1.0 + exp((v - (-43.0)) / 8.0));
}

static double inf_h_Kd ( const double v ){
  return 1.0 / (1.0 + exp((v - (-67.0)) / 7.3));
}

static double tau_m_Kd ( const double v ){
  return 1.0;
}

static double tau_h_Kd ( const double v ){
  return 1500.0;
}


static double alpha_m_Im ( const double v ){
  return 3.3e-3*exp(2.5*0.04*(v - (-35.0)));
}

static double beta_m_Im ( const double v ){
  return 3.3e-3*exp(-2.5*0.04*(v - (-35.0)));
}

static double inf_m_Im ( const double v ){
  return alpha_m_Im(v)/(alpha_m_Im(v)+beta_m_Im(v));
}

static double tau_m_Im ( const double v ){
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius - 21.0)/10.0 );
  return (1.0/(alpha_m_Im(v)+beta_m_Im(v)))/qt;
}


static double alpha_m_Imv2 ( const double v ){
  return 0.007 * exp( (6.0 * 0.4 * (v - (-48.0))) / 26.12 );
}

static double beta_m_Imv2 ( const double v ){
  return 0.007 * exp( (-6.0 * (1.0 - 0.4) * (v - (-48.0))) / 26.12 );
}

static double inf_m_Imv2 ( const double v ){
  return alpha_m_Imv2(v) / (alpha_m_Imv2(v) + beta_m_Imv2(v));
}

static double tau_m_Imv2 ( const double v ){
  const double celcius = 34.0;
  const double qt = pow (2.3, (celcius-30.0) / 10.0 );
  return  (15.0 + 1.0 / (alpha_m_Imv2(v) + beta_m_Imv2(v))) / qt;
}

static double alpha_m_Ih ( const double v ){
  return 0.001 * 6.43 * vtrap(v + 154.9, 11.9);
  //return 6.43 * vtrap(v + 154.9, 11.9);
}

static double beta_m_Ih ( const double v ){
  return 0.001 * 193.0 * exp( v /33.1);
  //return 193.0 * exp( v /33.1);
}

static double inf_m_Ih ( const double v ){
  return alpha_m_Ih(v) / (alpha_m_Ih(v) + beta_m_Ih(v));
}

static double tau_m_Ih ( const double v ){
  return 1.0 / (alpha_m_Ih(v) + beta_m_Ih(v) );
}


static double inf_z_SK ( const double v, double ca ){
  if(ca < 1e-7){
    ca = ca + 1e-07;
  }
  return 1.0 / (1.0 + pow ((0.00043 / ca), 4.8));
}

static double tau_z_SK ( const double v ){
  return 1.0;
}


static double alpha_m_CaHVA ( const double v ){
  return 0.055 * vtrap(-27.0 - v, 3.8);
}
static double beta_m_CaHVA ( const double v ){
  return (0.94 * exp((-75.0 - v) / 17.0));
}
static double inf_m_CaHVA ( const double v ){
  return alpha_m_CaHVA(v) / (alpha_m_CaHVA(v) + beta_m_CaHVA(v));
}
static double tau_m_CaHVA ( const double v ){
  return 1.0 / (alpha_m_CaHVA(v) + beta_m_CaHVA (v));
}

static double alpha_h_CaHVA ( const double v ){
  return (0.000457 * exp((-13.0 - v) / 50.0));
}
static double beta_h_CaHVA ( const double v ){
  return (0.0065 / (exp((-v - 15.0) / 28.0) + 1.0));
}
static double inf_h_CaHVA ( const double v ){
  return alpha_h_CaHVA(v) / (alpha_h_CaHVA(v) + beta_h_CaHVA(v));
}
static double tau_h_CaHVA ( const double v ){
  return 1.0 / (alpha_h_CaHVA(v) + beta_h_CaHVA(v));
}


static double inf_m_CaLVA ( const double v ){
  const double v_new = v + 10.0;
  return 1.0 / (1.0+ exp((v_new - (-30.0)) / -6.0));
}

static double inf_h_CaLVA ( const double v ){
  const double v_new = v + 10.0;
  return 1.0 / (1.0+ exp((v_new - (-80.0))/6.4));
}

static double tau_m_CaLVA ( const double v ){
  const double v_new = v + 10.0;
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius-21.0) / 10.0 );
  return (5.0 + 20.0 / (1.0 + exp((v_new - (-25.0)) / 5.0))) / qt;
}

static double tau_h_CaLVA ( const double v ){
  const double v_new = v + 10.0;
  const double celcius = 34.0;
  const double qt = pow ( 2.3, (celcius-21.0) / 10.0 );
  return (20.0 + 50.0 / (1.0 + exp((v_new - (-40.0)) / 7.0))) / qt;
}

//
// Nav
//
#define DT_NaV MIN( DT, 0.01 ) // [ms]
#define TMP_ITER_NaV ( DT / DT_NaV )
#define ITER_NaV MAX(1, TMP_ITER_NaV)
/*
static void swap_row( int rA, int rB, int n, double A [ ] [ n ], double b [ ] ){ 
  for ( int i = 0; i < n; i++ ) {
    double tmp = A [ rA ] [ i ];
    A [ rA ] [ i ] = A [ rB ] [ i ];
    A [ rB ] [ i ] = tmp;
  }
  double tmp = b [ rA ];
  b [ rA ] = b [ rB ];
  b [ rB ] = tmp;     
}

static int pivoting( int p, int n, const double A [ ] [ n ] ){
  int max = p;
  for ( int i = p + 1; i < n; i++ ) {
    if( fabs( A [ i ] [ p ] ) > fabs( A [ max ] [ p ] ) ){
      max = i;
    }
  }
  return max;
}
*/
static void forward_elimination ( int n, double A [ ] [ n ], double b [ ] ) {
  for ( int i = 0; i < n; i++ ) {
    //int ip = pivoting ( i, n, A );
    //if ( i != ip ) { 
    //  swap_row (ip, i, n, A, b ); 
    //}  
    //if ( fabs ( A [ i ] [ i ] ) < 1.0e-10 ) { exit ( 1 ); }

    for ( int j = i + 1; j < n; j++ ) {
      double factor = A [ j ] [ i ] / A [ i ] [ i ];
      for ( int k = i; k < n; k++ ) {
        A [ j ] [ k ] -= factor * A [ i ] [ k ];
      }
      b [ j ] -= factor * b [ i ];
    }
  }
}

static void backward_substitution ( int n, double A [ ] [ n ], double b [ ] ) {
  for ( int i = n - 1; i >= 0; i-- ) {
    b [ i ] /= A [ i ] [ i ];
    A [ i ] [ i ] = 1.0;

    for( int j = i - 1; j >= 0; j-- ) {
      double factor = A [ j ] [ i ];
      A [ j ] [ i ] = 0.0;
      b [ j ] -= factor * b [ i ];
    }
  }
}

static void gaussian_elimination ( int n, double A [ ] [ n ], double b [ ] ) {
  forward_elimination ( n, A, b );
  backward_substitution ( n, A, b );
}

static void Init_Nav_param ( const double l_v, double matrix [ N_STATE_NAV ] [ N_STATE_NAV ] ) {
  const double qt   = 0.77889990199; // qt = q10^((celsius-37(degC))/10(degC))=pow(2.3,((34-37))/10) = 0.77889990199      
  double alphaShift =  qt * 400.0   * exp ( l_v /  24.0   ); // alpha   * exp ( v / x1 )
  double betaShift  =  qt * 12.0    * exp ( l_v / -24.0   ); // beta    * exp ( v / x2 ) * qt
  const double gammaShift =  qt * 250.0                          ; // gamma   * exp ( v / x3 ) * qt
  const double deltaShift =  qt * 60.0                           ; // delta * exp ( v / x4 ) * qt
  const double a = 2.51; //pow ( ( Oon  / Con  ), 0.25 );
  const double b = 5.32; //pow ( ( Ooff / Coff ), 0.25 );

  double c1c2 = 4.0  * alphaShift;
  double c2c3 = 3.0  * alphaShift;
  double c3c4 = 2.0  * alphaShift;
  double c4c5 = 1.0  * alphaShift;
  double c5oo = 1.0  * gammaShift;

  double i1i2 = 4.0  * alphaShift * a;
  double i2i3 = 3.0  * alphaShift * a;
  double i3i4 = 2.0  * alphaShift * a;
  double i4i5 = 1.0  * alphaShift * a;
  double i5i6 = 1.0  * gammaShift;

  double c1i1 = 0.01 * qt; // Con*qt
  double c2i2 = 0.01 * qt * a; 
  double c3i3 = 0.01 * qt * a * a;
  double c4i4 = 0.01 * qt * a * a * a;
  double c5i5 = 0.01 * qt * a * a * a * a;
  double ooi6 = 8.0  * qt; // Oon  [ 1 / ms ] * qt

  double ooc5 = 1.0  * deltaShift;
  double c5c4 = 4.0  * betaShift;
  double c4c3 = 3.0  * betaShift;
  double c3c2 = 2.0  * betaShift;
  double c2c1 = 1.0  * betaShift;

  double i6i5 = 1.0 * deltaShift; 
  double i5i4 = 4.0 * betaShift / b;
  double i4i3 = 3.0 * betaShift / b;
  double i3i2 = 2.0 * betaShift / b;
  double i2i1 = 1.0 * betaShift / b;

  double i1c1 = 40.0  * qt; // Coff*qt
  double i2c2 = 40.0  * qt / ( b );
  double i3c3 = 40.0  * qt / ( b * b );
  double i4c4 = 40.0  * qt / ( b * b * b );
  double i5c5 = 40.0  * qt / ( b * b * b * b );
  double i6oo = 0.05  * qt; // Ooff [ 1 / ms ] * qt;

  matrix[0][0] = 0-(c1c2+c1i1); matrix[0][1] = c2c1;               matrix[0][5] = i1c1;
  matrix[1][0] = c1c2;          matrix[1][1] = 0-(c2c1+c2c3+c2i2); matrix[1][2] = c3c2; matrix[1][6] = i2c2;
  matrix[2][1] = c2c3;          matrix[2][2] = 0-(c3c2+c3c4+c3i3); matrix[2][3] = c4c3; matrix[2][7] = i3c3;
  matrix[3][2] = c3c4;          matrix[3][3] = 0-(c4c3+c4c5+c4i4); matrix[3][4] = c5c4; matrix[3][8] = i4c4;
  matrix[4][3] = c4c5;          matrix[4][4] = 0-(c5c4+c5oo+c5i5);                      matrix[4][9] = i5c5; matrix[4][11] = ooc5;
  matrix[5][0] = c1i1;          matrix[5][5] = 0-(i1c1+i1i2);      matrix[5][6] = i2i1;
  matrix[6][1] = c2i2;          matrix[6][5] = i1i2;               matrix[6][6] = 0-(i2i1+i2i3+i2c2); matrix[6][7] = i3i2;
  matrix[7][2] = c3i3;          matrix[7][6] = i2i3;               matrix[7][7] = 0-(i3i2+i3i4+i3c3); matrix[7][8] = i4i3;
  matrix[8][3] = c4i4;          matrix[8][7] = i3i4;               matrix[8][8] = 0-(i4i3+i4i5+i4c4); matrix[8][9] = i5i4;
  matrix[9][4] = c5i5;          matrix[9][8] = i4i5;               matrix[9][9] = 0-(i5i4+i5i6+i5c5); matrix[9][10] =i6i5;
  matrix[10][9]= i5i6;          matrix[10][10] = 0-(i6i5+i6oo);    matrix[10][11]= ooi6;
  matrix[11][4]= c5oo;          matrix[11][10] = i6oo;             matrix[11][11] = 0-(ooc5+ooi6);
}

static void Set_Nav_param ( const double l_v, double matrix [ N_STATE_NAV ] [ N_STATE_NAV ] ) {
  double qt   = 0.77889990199; // qt = q10^((celsius-37(degC))/10(degC))=pow(2.3,((34-37))/10) = 0.77889990199      
  double alphaShift =  qt * 400.0   * exp ( l_v /  24.0   ); // alpha   * exp ( v / x1 )
  double betaShift  =  qt * 12.0    * exp ( l_v / -24.0   ); // beta    * exp ( v / x2 ) * qt
  double gammaShift =  qt * 250.0                          ; // gamma   * exp ( v / x3 ) * qt
  double deltaShift =  qt * 60.0                           ; // delta * exp ( v / x4 ) * qt
  double a = 2.51; //pow ( ( Oon  / Con  ), 0.25 );
  double b = 5.32; //pow ( ( Ooff / Coff ), 0.25 );

  double c1c2 = - DT_NaV * 4.0  * alphaShift;
  double c2c3 = - DT_NaV * 3.0  * alphaShift;
  double c3c4 = - DT_NaV * 2.0  * alphaShift;
  double c4c5 = - DT_NaV * 1.0  * alphaShift;
  double c5oo = - DT_NaV * 1.0  * gammaShift;

  double i1i2 = - DT_NaV * 4.0  * alphaShift * a;
  double i2i3 = - DT_NaV * 3.0  * alphaShift * a;
  double i3i4 = - DT_NaV * 2.0  * alphaShift * a;
  double i4i5 = - DT_NaV * 1.0  * alphaShift * a;
  double i5i6 = - DT_NaV * 1.0  * gammaShift;

  double c1i1 = - DT_NaV * 0.01 * qt; // Con*qt
  double c2i2 = - DT_NaV * 0.01 * qt * a; 
  double c3i3 = - DT_NaV * 0.01 * qt * a * a;
  double c4i4 = - DT_NaV * 0.01 * qt * a * a * a;
  double c5i5 = - DT_NaV * 0.01 * qt * a * a * a * a;
  double ooi6 = - DT_NaV * 8.0  * qt; // Oon  [ 1 / ms ] * qt

  double ooc5 = - DT_NaV * 1.0  * deltaShift;
  double c5c4 = - DT_NaV * 4.0  * betaShift;
  double c4c3 = - DT_NaV * 3.0  * betaShift;
  double c3c2 = - DT_NaV * 2.0  * betaShift;
  double c2c1 = - DT_NaV * 1.0  * betaShift;

  double i6i5 = - DT_NaV * 1.0 * deltaShift; 
  double i5i4 = - DT_NaV * 4.0 * betaShift / b;
  double i4i3 = - DT_NaV * 3.0 * betaShift / b;
  double i3i2 = - DT_NaV * 2.0 * betaShift / b;
  double i2i1 = - DT_NaV * 1.0 * betaShift / b;

  double i1c1 = - DT_NaV * 40.0  * qt; // Coff*qt
  double i2c2 = - DT_NaV * 40.0  * qt / ( b );
  double i3c3 = - DT_NaV * 40.0  * qt / ( b * b );
  double i4c4 = - DT_NaV * 40.0  * qt / ( b * b * b );
  double i5c5 = - DT_NaV * 40.0  * qt / ( b * b * b * b );
  double i6oo = - DT_NaV * 0.05  * qt; // Ooff [ 1 / ms ] * qt;

  //                  c1,                c2,                c3,                c4,                c5,                i1,                i2,                i3,                i4,               i5,               i6,                 oo
  double l_mat [ N_STATE_NAV ] [ N_STATE_NAV ] = {
  /*c1*/  {1-(c1c2+c1i1),              c2c1,               0.0,               0.0,               0.0,              i1c1,               0.0,               0.0,               0.0,               0.0,              0.0,               0.0 },
  /*c2*/  {         c1c2,1-(c2c1+c2c3+c2i2),              c3c2,               0.0,               0.0,               0.0,              i2c2,               0.0,               0.0,               0.0,              0.0,               0.0 },
  /*c3*/  {          0.0,              c2c3,1-(c3c2+c3c4+c3i3),              c4c3,               0.0,               0.0,               0.0,              i3c3,               0.0,               0.0,              0.0,               0.0 },
  /*c4*/  {          0.0,               0.0,              c3c4,1-(c4c3+c4c5+c4i4),              c5c4,               0.0,               0.0,               0.0,              i4c4,               0.0,              0.0,               0.0 },
  /*c5*/  {          0.0,               0.0,               0.0,              c4c5,1-(c5c4+c5oo+c5i5),               0.0,               0.0,               0.0,               0.0,              i5c5,              0.0,              ooc5 },
  /*i1*/  {         c1i1,               0.0,               0.0,               0.0,               0.0,     1-(i1c1+i1i2),              i2i1,               0.0,               0.0,               0.0,              0.0,               0.0 },
  /*i2*/  {          0.0,              c2i2,               0.0,               0.0,               0.0,              i1i2,1-(i2i1+i2i3+i2c2),              i3i2,               0.0,               0.0,              0.0,               0.0 },
  /*i3*/  {          0.0,               0.0,              c3i3,               0.0,               0.0,               0.0,              i2i3,1-(i3i2+i3i4+i3c3),              i4i3,               0.0,              0.0,               0.0 },
  /*i4*/  {          0.0,               0.0,               0.0,              c4i4,               0.0,               0.0,               0.0,              i3i4,1-(i4i3+i4i5+i4c4),              i5i4,              0.0,               0.0 },
  /*i5*/  {          0.0,               0.0,               0.0,               0.0,              c5i5,               0.0,               0.0,               0.0,              i4i5,1-(i5i4+i5i6+i5c5),             i6i5,               0.0 },
  /*i6*/  {          0.0,               0.0,               0.0,               0.0,               0.0,               0.0,               0.0,               0.0,               0.0,              i5i6,    1-(i6i5+i6oo),              ooi6 },
  /*oo*/  {          0.0,               0.0,               0.0,               0.0,              c5oo,               0.0,               0.0,               0.0,               0.0,               0.0,             i6oo,     1-(ooc5+ooi6) }
  };
  for ( int i = 0; i < N_STATE_NAV; i++ ) { 
    for ( int j = 0; j < N_STATE_NAV; j++ ) { 
      matrix [ i ] [ j ] = l_mat [ i ] [ j ]; 
    }
  }   
}


// First order
static void Nav_update ( const double l_v, double *oo, double *c1, double *c2, double *c3, double *c4, double *c5, double *i1, double *i2, double *i3, double *i4, double *i5, double *i6 ) {
  
  double vca [ N_STATE_NAV ] [ N_STATE_NAV ] = {};
  double vcb [ N_STATE_NAV ] = { c1 [ 0 ] , c2 [ 0 ], c3 [ 0 ], c4 [ 0 ], c5 [ 0 ], i1 [ 0 ], i2 [ 0 ], i3 [ 0 ], i4 [ 0 ], i5 [ 0 ], i6 [ 0 ], oo [ 0 ] };  
  double matrixA [ N_STATE_NAV ] [ N_STATE_NAV ] = {};
  //GR Na channel parameters
  Set_Nav_param ( l_v, matrixA );
  for ( int iter = 0; iter < ( int ) ITER_NaV; iter++ ) {
    for ( int i = 0; i < N_STATE_NAV; i++ ) { 
      for ( int j = 0; j < N_STATE_NAV; j++ ) { 
        vca [ i ] [ j ] = matrixA [ i ] [ j ]; 
      }
    }     

    // X(t+dt) = X(t)    + dt * vca * X(t+dt), ( X:= {c1, c2, ..., o} = vcb )
    // X(t)    = X(t+dt) - dt * vca * X(t+dt), 
    // X(t)    =   ( I - dt * vca ) * X(t+dt) 
    //for ( int i = 0; i < N_STATE_NAV; i++ ) {
    //  for ( int j = 0; j < N_STATE_NAV; j++ ) { 
    //    vca [ i ] [ j ] *= - DT_NaV; // *= - DT
    //  }
    //}
    //for ( int i = 0; i < N_STATE_NAV; i++ ) { vca [ i ] [ i ] += 1.0; }
    
    for ( int col = 0; col < N_STATE_NAV; col++ ) { vca [ N_STATE_NAV - 1 ] [ col ] = 1.0; } // I6 channel 
    vcb [ N_STATE_NAV - 1 ] = 1.0;
    gaussian_elimination ( N_STATE_NAV, vca, vcb );
  }      
  c1 [ 0 ] = vcb [ 0 ];
  c2 [ 0 ] = vcb [ 1 ]; 
  c3 [ 0 ] = vcb [ 2 ];
  c4 [ 0 ] = vcb [ 3 ];
  c5 [ 0 ] = vcb [ 4 ];
  i1 [ 0 ] = vcb [ 5 ];
  i2 [ 0 ] = vcb [ 6 ];
  i3 [ 0 ] = vcb [ 7 ];
  i4 [ 0 ] = vcb [ 8 ];
  i5 [ 0 ] = vcb [ 9 ];
  i6 [ 0 ] = vcb [ 10 ];
  oo [ 0 ] = vcb [ 11 ];

  // Debug
  //double sum = ( c1 [ 0 ] + c2 [ 0 ] + c3 [ 0 ] + c4 [ 0 ] + c5 [ 0 ] + i1 [ 0 ] + i2 [ 0 ] + i3 [ 0 ] + i4 [ 0 ] + i5 [ 0 ] + i6 [ 0 ] + oo [ 0 ] );
  //if ( ( sum > 1.0000001 ) || ( sum < 0.9999999 ) ) {
  //  //double error = 1.0 - sum;
  //  //i6 [ 0 ] -= error;
  //  //if ( i6 [ 0 ] < 0.0 ) {
  //    printf ( "Nav11 error %.15f\n", 1.0 - sum ); 
  //    printf ("c1 [ 0 ] = %f\nc2 [ 0 ] = %f\nc3 [ 0 ] = %f\nc4 [ 0 ] = %f\nc5 [ 0 ] = %f\ni1 [ 0 ] = %f\ni2 [ 0 ] = %f\ni3 [ 0 ] = %f\ni4 [ 0 ] = %f\ni5 [ 0 ] = %f\ni6 [ 0 ] = %f\noo [ 0 ] = %f\n",c1 [ 0 ],c2 [ 0 ], c3 [ 0 ], c4 [ 0 ], c5 [ 0 ], i1 [ 0 ], i2 [ 0 ], i3 [ 0 ], i4 [ 0 ], i5 [ 0 ], i6 [ 0 ], oo [ 0 ] ); 
  //    exit ( 1 );
  //  //}
  //}
}
