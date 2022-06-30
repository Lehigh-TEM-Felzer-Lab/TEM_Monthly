/* *************************************************************
TCLMDAT45.H - object to read and write the structure of the
                climate data from files used by the Terrestrial
                Ecosystem Model (TEM)

20060114 - DWK created by modifying tclmdat425.h
20060114 - DWK changed class Clmdata to class Clmdata43
20060114 - DWK changed long year to int year in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed char varname[9] to string varname in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed char contnent[9] to string contnent in
           out(), outdel(), pctout() and pctoutdel()
20060114 - DWK changed public variable char varname[9] to
           string varname
20060114 - DWK changed public variable char contnent[9] to
           string contnent
20060114 - DWK deleted tclmdat425.cpp at bottom of file

************************************************************* */

#ifndef TCLMDAT45_H
#define TCLMDAT45_H

#include "temconsts45.hpp"

class Clmdat45
{

  public:

     Clmdat45( void );

/* *************************************************************
		      Public Functions
************************************************************* */

     // read data structure.

     int get( ifstream& ifile );
     int getdel( FILE* ifile );

     //write data structure.

     void out( ofstream& ofile,
               const float& col,
               const float& row,
               const string& varname,
               const int& carea,
               const int& year,
               double mon[CYCLE],
               const string& contnent );

     void outdel( ofstream& ofile,
                  const float& col,
                  const float& row,
                  const string& varname,
                  const int& carea,
                  const int& year,
                  double mon[CYCLE],
                  const string& contnent );

     void pctout( ofstream& ofile,
                  const float& col,
                  const float& row,
                  const string& varname,
                  const int& carea,
                  const int& year,
                  double mon[CYCLE],
                  const string& contnent );

     void poutdel( ofstream& ofile,
                   const float& col,
                   const float& row,
                   const string& varname,
                   const int& carea,
                   const int& year,
                   double mon[CYCLE],
                   const string& contnent );


/* *************************************************************
		     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)
     float col;

     // row or latitude of grid cell (degrees)
     float row;

     // climate variable name
     string varname;

     // area covered by grid cell (sq. km)
     int carea;

      // date (year) of data
     //long year;
     int year;

     // annual sum of monthly data for grid cell
     double total;

      // maximum monthly value for grid cell
     double max;

     // mean annual value for grid cell
     double ave;

     // minimum monthly value for grid cell
     double min;

     // monthly values for the grid cell
     double mon[CYCLE];

      // name of continent containing grid cell
     string contnent;


  private:

/* *************************************************************
		      Private Variables
************************************************************* */

     int clmend;
     long curpos;
     long lagpos;

};

#endif
