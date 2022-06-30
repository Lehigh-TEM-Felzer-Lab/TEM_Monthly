/* *************************************************************
disturbdat.CPP - object to read and write the structure of
                   land use/land cover data from/to files

Modifications:

20060422 - DWK created by modifying lulcdat437.cpp
20060422 - DWK changed include from lulcdat437.h to lulcdat44.h
20060422 - DWK changed Lulcdata43:: to Lulcdata44::
20060422 - DWK added int disturbflag and int disturbmonth to
           functions
20070105 - TWC renamed to lulcdat45

************************************************************* */

#include<cstdio>

  using std::fscanf;
  using std::FILE;

#include<iostream>

  using std::ios;
  using std::endl;

#include<fstream>

  using std::ifstream;
  using std::ofstream;

#include<iomanip>

  using std::setprecision;

#include<string>

  using std::string;


#include "disturbdat.h"


Disturbdat::Disturbdat( void )
{

  disturbend = 1;
  lagpos = -99;
  curpos = 0;

};

/* **************************************************************
                    Public Functions
************************************************************** */

/* *************************************************************
************************************************************* */


/* *************************************************************
************************************************************* */

int Disturbdat::getdel( FILE* infile )
{
  char tmpvarname[40];

  disturbend = fscanf( infile,
                    "%f,%f, %s ,%lf",
                    &col,
                    &row,
                    tmpvarname,
                    &retint );

  varname = tmpvarname;

  return disturbend;

};

/* *************************************************************
************************************************************* */
