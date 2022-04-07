/* *************************************************************
TTEMDAT45A.H - object to read and write the structure of the
               output data from files generated by the
               Terrestrial Ecosystem Model (TEM)

Modifications:

20031019 - DWK created by modifying ttemdat425.h
20031019 - DWK changed class Temdata to Temdata43
20031019 - DWK changed char varname[9] to string varname
20031019 - DWK changed char contnent[9] to string contnent
20031019 - DWK changed long year to int year
20040531 - DWK changed public int carea to long carea
20040531 - DWK added public long subarea
20040830 - DWK added public int currentveg and int subtype
20040830 - DWK renamed public int tveg to be int potveg
20040830 - DWK renamed public int subtveg to be int cmnt
20040830 - DWK renamed string contnent to be string region
20051117 - DWK added include temconsts43.hpp
20051117 - DWK deleted ttemdat425.cpp from bottom of file
20060414 - DWK added public int icohort
20070105 - TWC changed name to ttemdat45

************************************************************* */

#ifndef TTEMDAT45A_H
#define TTEMDAT45A_H

#include "temconsts45.hpp"

class Temdat45
{

public:
     Temdat45(void);

     /* *************************************************************
                Public Functions
     ************************************************************* */

     int get(ifstream &ifile);

     int getdel(FILE *ifile);

     void getyrsum(double months[CYCLE]);

     void out(ofstream &ofile,
              const float &col,
              const float &row,
              const string &varname,
              const int &icohort,
              const int &standage,
              const int &potveg,
              const int &currentveg,
              const int &subtype,
              const int &cmnt,
              const double &psiplusc,
              const int &qlcon,
              const long &carea,
              const long &subarea,
              const int &year,
              double mon[CYCLE],
              const string &region);

     void outdel(ofstream &ofile,
                 const float &col,
                 const float &row,
                 const string &varname,
                 const int &icohort,
                 const int &standage,
                 const int &potveg,
                 const int &currentveg,
                 const int &subtype,
                 const int &cmnt,
                 const double &psiplusc,
                 const int &qlcon,
                 const long &carea,
                 const long &subarea,
                 const int &year,
                 double mon[CYCLE],
                 const string &region);

     void pctout(ofstream &ofile,
                 const float &col,
                 const float &row,
                 const string &varname,
                 const int &icohort,
                 const int &standage,
                 const int &potveg,
                 const int &currentveg,
                 const int &subtype,
                 const int &cmnt,
                 const double &psiplusc,
                 const int &qlcon,
                 const long &carea,
                 const long &subarea,
                 const int &year,
                 double mon[CYCLE],
                 const string &region);

     void poutdel(ofstream &ofile,
                  const float &col,
                  const float &row,
                  const string &varname,
                  const int &icohort,
                  const int &standage,
                  const int &potveg,
                  const int &currentveg,
                  const int &subtype,
                  const int &cmnt,
                  const double &psiplusc,
                  const int &qlcon,
                  const long &carea,
                  const long &subarea,
                  const int &year,
                  double mon[CYCLE],
                  const string &region);

     /* *************************************************************
                Public Variables
     ************************************************************* */

     /* NOTE:  substitute "daily" for "monthly" if CYCLE = 365 rather
               than CYCLE = 12 */

     // mean annual value for the grid cell
     double ave;

     // area covered by grid cell (square kilometers)
     long carea;

     // Index for community type
     int cmnt;

     // column or longitude of grid cell (degrees)
     float col;

     // Current vegetation type in grid cell
     int currentveg;

     // Cohort index
     int icohort;

     // maximum monthly value for the grid cell
     double max;

     // minimum monthly value for the grid cell
     double min;

     // monthly values of the TEM output variable for the
     //   grid cell

     double mon[CYCLE];

     // Potential vegetation type of grid cell
     //   (categorical data)
     int potveg;

     // percent silt plus percent clay of the grid cell
     double psiplusc;

     // quality control flag
     int qlcon;

     // name of region containing grid cell
     string region;

     // row or latitude of grid cell (degrees)
     float row;

     // Stand age of cohort
     int standage;

     // area covered by cohort in a grid cell
     //   (square kilometers)
     long subarea;

     // vegetation subtype of grid cell (categorical data)
     int subtype;

     // annual sum of monthly data for the grid cell
     double total;

     // name of TEM output variable
     string varname;

     // date (year) of data (transient version of TEM)
     //   or number of years until equilibrium are met
     //  (equilibrium version of TEM)
     int year;

private:
     /* ************************************************************
                     Private Variables
     ************************************************************ */

     int temend;
     long curpos;
     long lagpos;
};

#endif
