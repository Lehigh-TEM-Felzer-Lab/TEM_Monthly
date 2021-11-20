//
/* **************************************************************
*****************************************************************
ATMS45.CPP - object describes physical characteristics of the
	         atmosphere

Modifications:

20060126 - DWK created by modifying atms50b5.cpp
20060126 - DWK changed include from atms50b5.h to atms437.h
20060126 - DWK changed Atmosphere50:: to class Atmsflux::
20070105 - TWC changed Atmsflux:: to class Atms45
2007 - TWC/BSF: Summary
  using std:sin
  Add three new function: daylength

****************************************************************
************************************************************* */

#include<cmath>
  using std::exp;
  using std::sin;
  using std::cos;
  using std::pow;
  using std::acos;
  using std::tan;
  using std::fabs;

#include<iostream>

    using std::cout;
    using std::ios;
    using std::cerr;
    using std::endl;


/* *************************************************************
************************************************************* */

#include "atms45_ndep.h"

Atms45::Atms45()
{
  // Number of days per month

  ndays[0] = 31;
  ndays[1] = 28;
  ndays[2] = 31;
  ndays[3] = 30;
  ndays[4] = 31;
  ndays[5] = 30;
  ndays[6] = 31;
  ndays[7] = 31;
  ndays[8] = 30;
  ndays[9] = 31;
  ndays[10] = 30;
  ndays[11] = 31;

};

/* **************************************************************
************************************************************** */

double Atms45::daylength( const float& lat, const int& dm )
{

//  double yd[CYCLE],ampl;
  double decldeg, sunrise_angle, sumdayl, sumdaymonth;
  int ymd[CYCLE];
  int initday, mday, jday;
  
  double pi = 3.14159265359;

/*  yd[0] = 15;
  yd[1] = 46;
  yd[2] = 74;
  yd[3] = 105;
  yd[4] = 135;
  yd[5] = 166;
  yd[6] = 196;
  yd[7] = 227;
  yd[8] = 258;
  yd[9] = 288;
  yd[10] = 319;
  yd[11] = 349; 
  
  ampl = exp(7.42+0.045*lat)/3600.;
  dayl = ampl*(sin((yd[dm]-79.)*0.01721))+12.0;
  //old formula from Running & Coughlan 1988, technically incorrect
*/ 
  
  ymd[0] = 0;
  ymd[1] = 31;
  ymd[2] = 59;
  ymd[3] = 90;
  ymd[4] = 120;
  ymd[5] = 151;
  ymd[6] = 181;
  ymd[7] = 212;
  ymd[8] = 243;
  ymd[9] = 273;
  ymd[10] = 304;
  ymd[11] = 334;
  
  initday = ymd[dm];
  sumdayl = ZERO;
  sumdaymonth = ZERO;

  for(mday = 0; mday < ndays[dm]; ++mday)
  {
    sumdaymonth += 1.0;
    jday = initday + mday;
    decldeg = -23.4856*cos(2.0*pi*(jday + 10.0)/365.25);
    if( (90.0 - fabs(lat)) > fabs(decldeg) ) // then not in region where day is 0 or 24 hours
    {
      sunrise_angle = acos( -tan(pi*lat/180.0)*tan(pi*decldeg/180.0) );
      sumdayl += sunrise_angle*(24.0/pi);
    }
    else // then in region where day is 0 or 24 hours
    {
      if( lat * decldeg >= 0.0 ) { sumdayl += 23.9999; }
      else { sumdayl += 0.0001; }
      // daylengths actually 24 and zero, but make very slightly different for comp. purposes
    }     
  }
  
//  dayl = sumdayl / (double) ndays[dm];
  dayl = sumdayl / sumdaymonth;         
  

//  cout << "diag = " << dm << " " << dayl << " " << lat << endl;

  return dayl;

};

/* *************************************************************
************************************************************* */
/* **************************************************************
************************************************************** */
double Atms45::lwrad( const double& tair,
                      const double& vpr,
                      const double& nirr,
                      const double& girr )

{

    double epsc, ccor;
    double lwout;
    
    epsc = 1.24 * pow((vpr/(tair+273.15)),1.0/7.0);

    ccor = 1.6 * (nirr/girr) - 0.2;

    if (ccor < 0.2) {ccor = 0.2;}
    if (ccor > 1.0) {ccor = 1.0;}

    lwout = (1.0 - epsc) * ccor * 5.67e-8 * pow((tair+273.15),4.);
    
    return lwout;

};

/* **************************************************************
************************************************************** */

/* *************************************************************
************************************************************* */

void Atms45::petjh( const double& nirr,
                          const double& tair,
                          const int& pdm )
{

  double f;
  double rt;

  f = ((9.0/5.0) * tair) + 32.0;
  rt = nirr * 0.016742;
  pet = ((0.014*f) - 0.37) * rt * ndays[pdm];

  if ( pet <= ZERO ) { pet = 0.001; }

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atms45::precsplt( const double& prec,
                             const double& tair,
                             double& rain,
                             double& snowfall )
{


/* *************************************************************
	Willmott's assumptions on snow/rain split:
************************************************************** */

  if ( tair >= -1.0 )
  {
    rain = prec;
    snowfall = ZERO;
  }
  else
  {
    rain = ZERO;
    snowfall = prec;
  }

};

/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atms45::resetMonthlyFluxes( void )
{
  // Reset monthly fluxes to zero

  pet = ZERO;

};
/* *************************************************************
************************************************************* */

/* *************************************************************
************************************************************* */

void Atms45::resetYrFluxes( void )
{
  // Reset annual fluxes and summary variables to zero

  yrrain = ZERO;
  yrsnowfall = ZERO;
  yrpet = ZERO;

};


