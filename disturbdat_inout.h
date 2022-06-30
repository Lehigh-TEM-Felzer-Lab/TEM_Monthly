/* *************************************************************
DISTURB.H - object to read and write the structure of land
                  use / land cover cohort data from/to files

20060422 - DWK created by modifying lulcdat437.h
20060422 - DWK changed class Lulcdata43 to class Lulcdata44
20060422 - DWK deleted public double RAP
20060422 - DWK added public int disturbflag and public int
           disturbmonth
20070105 - TWC changed name to lulcdat45

************************************************************* */

#ifndef DISTURBDAT_INOUT_H
#define DISTURBDAT_INOUT_H

class Disturbdat
{

  public:

     Disturbdat( void );

/* *************************************************************
                      Public Functions
************************************************************* */

// read data structure.
     int getdel( FILE* infile );

     void outdel( ofstream& ofile,
           const float& col,
           const float& row,
           const string& varname,
           const float& retint );


/* *************************************************************
                     Public Variables
************************************************************* */

     // column or longitude of grid cell (degrees)
     float col;

     // row or latitude of grid cell (degrees)
     float row;

     // varname = "STORM" for this data set
     string varname;

     // return interval
     double retint;


  private:

/* *************************************************************
                      Private Variables
************************************************************* */

     int disturbend;
     long curpos;
     long lagpos;

};

#endif

