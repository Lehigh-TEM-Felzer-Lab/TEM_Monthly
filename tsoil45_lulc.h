/* **************************************************************
*****************************************************************
TSOIL45_LULC.H - object describing characteristics of soil used by
	        the Terrestrial Ecosystem Model (TEM)

Modifications:

20060126 - DWK created by modifying tsoil50b5.h
20060126 - DWK changed include from temconts51.hpp to
           temconsts43.hpp
20060126 - DWK changed include from tprocessXML51.h to
           tprocessXML437.h
20060126 - DWK deleted include qsoiltemp50b5.h
20060126 - DWK changed class Tsoil50 to class Tsoilflux
20060126 - DWK changed inheritance from ProcessXML50 to
           ProcessXML43
20060126 - DWK deleted public function double
           updateActiveLayerRootZ(), inline double
           getACTLAYER(), inline double getDST10(), inline
           double getNEXTDST10(), inline double getPREVDST10()
           and inline double getTSOIL()
20060126 - DWK deleted  public Soilthermal50 stm, int stmflg,
           double yrdst10 and double yrtsoil
20060126 - DWK deleted private double activeLayer, double dst10,
           double nextdst10, double prevdst10 and double tsoil
20060129 - DWK added public inline getVSM() and setVSM()
20060320 - DWK added private double surfrun and double drainage
20070105 - TWC changed name to tsoil45
2007:  TWC/BSF
     setSWP public function
	 get/set: evap
	 Public var.: yrevap
	 Private var: evap

****************************************************************
************************************************************* */

#ifndef TSOIL45_LULC_H
#define TSOIL45_LULC_H

// Tsoil45 uses the global constant MAXCMNT
#include "temconsts45.hpp"

#include "tprocessXML45.h"

#include "bioms45.hpp"      // Tsoilflux uses Biomass class


class Tsoil45 : public ProcessXML45
{

  public:

     Tsoil45( void );

/* **************************************************************
		 Public Functions
************************************************************** */


     void getrootz( const string& ecd );

     void lake( const double& tair,
                const double& prec,
                double& rain,
                double& snowfall,
                const double& pet,
                double& eet );


     void resetMonthlyFluxes( void );

     void resetYrFluxes( void );

     void setKH2O( const double& vsm,
                   const int& moistlim );
                        
     void setSWP( void );

     void showecd( void );

     double snowmelt( const double& elev,
                      const double& tair,
                      const double& prevtair,
                      const double& snowpack );


     void updateHydrology( const double& elev,
                           const double& tair,
                           const double& prevtair,
                           const double& prev2tair,
                           const double& rain,
                           const double& pet,
                           const double& avlh2o,
                           const double& rgrndh2o,
                           const double& sgrndh2o,
                           const int& irrgflag,
                           double& irrigate,
                           const int& pdm );

     void updateDOCLEACH( const double& doc, 
                          const double& sh2o );

     void updateNLosses( const int& pdcmnt,
                         const double& h2oloss,
                         const double& availn,
                         const double& soilh2o );


     double updateRootZ( const int& pdcmnt,
                         const double& sh2o,
                         const double& finerootc );



     void xtext( const int& pdcmnt,
                 const double& pctsilt,
                 const double& pctclay );


     // "Get" and "Set" private variables and parameters

     // availn *************************************************

     inline double getAVLN( void ) { return availn; }

     inline void setAVLN( const double& pavln )
     {
       availn = pavln;
     }


     // avlh2o *************************************************

     inline double getAVLH2O( void ) { return avlh2o; }

     inline void setAVLH2O( const double& pavlh2o )
     {
       avlh2o = pavlh2o;
     }


     // awcapmm ************************************************

     inline double getAWCAPMM( void ) { return awcapmm; }

     inline void setAWCAPMM( const double& pawcapmm )
     {
       awcapmm = pawcapmm;
     }

     // lchdon ************************************************

     inline double getLCHDON( void ) { return lchdon; }

     inline void setLCHDON( const double& plchdon )
     {
       lchdon = plchdon;
     }

     // lchdin ************************************************
     
           inline double getLCHDIN( void ) { return lchdin; }
     
           inline void setLCHDIN( const double& plchdin )
           {
               lchdin = plchdin;
           }
     

     // drainage ***********************************************

     inline double getDRAINAGE( void ) { return drainage; }

     inline void setDRAINAGE( const double& pdrainage )
     {
       drainage = pdrainage;
     }


     // eet ****************************************************

     inline double getEET( void ) { return eet; }

     inline void setEET( const double& peet )
     {
       eet = peet;
     }

     // fldcap *************************************************

     inline double getFLDCAP( void ) { return fldcap; }


     // fldcapa ************************************************

     inline double getFLDCAPA( void ) { return fldcapa; }

     inline void setFLDCAPA( const double& pfldcapa )
     {
       fldcapa = pfldcapa;
     }


     // fldcapb ************************************************

     inline double getFLDCAPB( void ) { return fldcapb; }

     inline void setFLDCAPB( const double& pfldcapb )
     {
       fldcapb = pfldcapb;
     }
     
     // gm ***************************************************

     inline double getGM( void ) { return gm; }


     // h2oyld *************************************************

     inline double getH2OYLD( void ) { return h2oyld; }

     inline void setH2OYLD( const double& ph2oyld )
     {
       h2oyld = ph2oyld;
     }


     // ineet **************************************************

     inline double getINEET( void ) { return ineet; }

     inline void setINEET( const double& pineet )
     {
       ineet = pineet;
     }


     // kh2o ***************************************************

     inline double getKH2O( void ) { return kh2o; }

     // lchdoc *************************************************

     inline double getLCHDOC( void ) { return lchdoc; }

     inline void setLCHDOC( const double& plchdoc )
     {
       lchdoc = plchdoc;
     }


     // minrootz ***********************************************

     inline double getMINROOTZ( const int& pcmnt )
     {
       return minrootz[pcmnt];
     }

     inline void setMINROOTZ( const double& pminrootz,
                              const int& pcmnt )
     {
       minrootz[pcmnt] = pminrootz;
     }

     // z0 ***********************************************

     inline double getZ0( const int& pcmnt )
     {
       return z0[pcmnt];
     }

     inline void setZ0( const double& pz0,
                              const int& pcmnt )
     {
       z0[pcmnt] = pz0;
     }



     // moist **************************************************

     inline double getMOIST( void ) { return moist; }

     inline void setMOIST( const double& psh2o )
     {
       moist = psh2o;
     }


     // ninput *************************************************

     inline double getNINPUT( void ) { return ninput; }

     inline void setNINPUT( const double& pninput )
     {
       ninput = pninput;
     }


     // nloss and denitr ***************************************

     inline double getNLOSS( const int& pcmnt )
     {
       return nloss[pcmnt];
     }

     inline void setNLOSS( const double& pnloss,
                           const int& pcmnt )
     {
       nloss[pcmnt] = pnloss;
     }

     inline double getDENITR( const int& pcmnt )
     {
       return denitr[pcmnt];
     }

     inline void setDENITR( const double& pdenitr,
                           const int& pcmnt )
     {
       denitr[pcmnt] = pdenitr;
     }


     // nlost **************************************************

     inline double getNLOST( void ) { return nlost; }

     inline void setNLOST( const double& pnlst )
     {
       nlost = pnlst;
     }


     // pctclay ************************************************

     inline double getPCTCLAY( void ) { return pctclay; }

     inline void setPCTCLAY( const double& ppctclay )
     {
       pctclay = ppctclay;
     }


     // pcfldcap ***********************************************

     inline double getPCTFLDCAP( void ) { return pcfldcap; }

     inline void setPCTFLDCAP( const double& ppctfldcap )
     {
       pcfldcap = ppctfldcap;
     }


     // pctp ***************************************************

     inline double getPCTP( void ) { return pctp; }

     inline void setPCTP( const double& ppctp )
     {
       pctp = ppctp;
     }

     // pctpor *************************************************

     inline double getPCTPOR( void ) { return pctpor; }


     // pctpora ************************************************

     inline double getPCTPORA( void ) { return pctpora; }

     inline void setPCTPORA( const double& ppctpora )
     {
       pctpora = ppctpora;
     }


     // pctporb ************************************************

     inline double getPCTPORB( void ) { return pctporb; }

     inline void setPCTPORB( const double& ppctporb )
     {
       pctporb = ppctporb;
     }


     // pctsand ************************************************

     inline double getPCTSAND( void ) { return pctsand; }

     inline void setPCTSAND( const double& ppctsand )
     {
       pctsand = ppctsand;
     }


     // pctsilt ************************************************

     inline double getPCTSILT( void ) { return pctsilt; }

     inline void setPCTSILT( const double& ppctsilt )
     {
       pctsilt = ppctsilt;
     }
     
     // pctwiltpt ************************************************

     inline double getPCTWILTPT( void ) { return pcwiltpt; }


     // prevspack **********************************************

     inline double getPREVSPACK( void ) { return prevspack; }

     inline void setPREVSPACK( const double& pprvspack )
     {
       prevspack = pprvspack;
     }


     // psiplusc ***********************************************

     inline double getPSIPLUSC( void ) { return psiplusc; }

     inline void setPSIPLUSC( const double& ppsiplusc )
     {
       psiplusc = ppsiplusc;
     }

     // reet **************************************************

     inline double getREET( void ) { return reet; }

     inline void setREET( const double& peet )
     {
       reet = peet;
     }


     // rootz **************************************************

     inline double getROOTZ( void ) { return rootz; }

     inline void setROOTZ( const double& prootz )
     {
       rootz = prootz;
     }


     // rootza *************************************************

     inline double getROOTZA( const int& pcmnt )
     {
       return rootza[pcmnt];
     }

     inline void setROOTZA( const double& prootza,
                            const int& pcmnt )
     {
       rootza[pcmnt] = prootza;
     }


     // rootzb *************************************************

     inline double getROOTZB( const int& pcmnt )
     {
       return rootzb[pcmnt];
     }

     inline void setROOTZB( const double& prootzb,
                            const int& pcmnt )
     {
       rootzb[pcmnt] = prootzb;
     }


     // rootzc *************************************************

     inline double getROOTZC( const int& pcmnt )
     {
       return rootzc[pcmnt];
     }

     inline void setROOTZC( const double& prootzc,
                            const int& pcmnt )
     {
       rootzc[pcmnt] = prootzc;
     }


     // rperc **************************************************

     inline double getRPERC( void ) { return rperc; }


     // rrun ***************************************************

     inline double getRRUN( void ) { return rrun; }

     inline void setRRUN( const double& prrun )
     {
       rrun = prrun;
     }
     
     // snowinf ************************************************

     inline double getSNOWINF( void ) { return snowinf; }

     inline void setSNOWINF( const double& psnwinf )
     {
       snowinf = psnwinf;
     }


     // snowpack ***********************************************

     inline double getSNOWPACK( void ) { return snowpack; }

     inline void setSNOWPACK( const double& psnwpck )
     {
       snowpack = psnwpck;
     }


     // sperc **************************************************

     inline double getSPERC( void ) { return sperc; }

     // Nonsymbiotic N fix ***************************************************

     inline double getSONINP( void ) { return soninp; }

     inline void setSONINP( const double& psoninp )
     {
       soninp = psoninp;
     }


     // srun ***************************************************

     inline double getSRUN( void ) { return srun; }

     inline void setSRUN( const double& psrun )
     {
       srun = psrun;
     }


     // surfrun ************************************************

     inline double getSURFRUN( void ) { return surfrun; }

     inline void setSURFRUN( const double& psurfrun )
     {
       surfrun = psurfrun;
     }
     
     // swp ****************************************************
     
     inline double getSWP( void ) { return swp; }
     
     // totpor *************************************************

     inline double getTOTPOR( void ) { return totpor; }

     inline void setTOTPOR( const double& ptotpor )
     {
       totpor = ptotpor;
     }


     // vsm *************************************************

     inline double getVSM( void ) { return vsm; }

     inline void setVSM( const double& pvsm )
     {
       vsm = pvsm;
     }

     // wiltpt *************************************************

     inline double getWILTPT( void ) { return wiltpt; }

     inline void setWILTPT( const double& pwiltpt )
     {
       wiltpt = pwiltpt;
     }


     // wiltpta ************************************************

     inline double getWILTPTA( void ) { return wiltpta; }

     inline void setWILTPTA( const double& pwiltpta )
     {
       wiltpta = pwiltpta;
     }


     // wiltptb ************************************************

     inline double getWILTPTB( void ) { return wiltptb; }

     inline void setWILTPTB( const double& pwiltptb )
     {
       wiltptb = pwiltptb;
     }


     // wsoil **************************************************

     inline int getWSOIL( void ) { return wsoil; }

     inline void setWSOIL( const int& pwsoil )
     {
       wsoil = pwsoil;
     }

/* *************************************************************
		 Public Variables
************************************************************* */

     // Annual sum of avlh2o
     double yravlh2o;

     // Annual sum of availn
     double yravln;

     // Ratio of soil reactive organic carbon to
     //   soil reactive organic nitrogen
     double yrc2n;

     // Annual estimated actual evapotranspiration (mm / year)
     double yreet;

     // Annual sum of h2oyld (mm / year)
     double yrh2oyld;

     // Annual initial estimated actual evapotranspiration
     //   (mm / year)
     double yrineet;

     // Annual sum of ninput
     double yrnin;     // (g N / (sq. meter * year))

     // Annual sum of nlost
     double yrnlost;   // (g N / (sq. meter * year))
     double yrlchdin;   // (g N / (sq. meter * year))

     // Annual sum of org.carbon
     double yrorgc;

     // Annual sum of org.nitrogen
     double yrorgn;

     // Annual sum of pctp
     double yrpctp;

     // Annual sum of rgrdnh2o
     double yrrgrndh2o;

     // Annual sum of rperc
     double yrrperc;   // (mm / year)

     // Annual sum of rrun
     double yrrrun;  // (mm / year)

     // Annual sum of sgrndh2o
     double yrsgrndh2o;

     // Annual sum of moist
     double yrsmoist;

     // Annual sum of snowinf
     double yrsnowinf;      // (mm / year)

     // Annual sum of snowpack
     double yrsnowpack;

     // Annual sum of sperc
     double yrsperc; // (mm / year)

     // Annual sum of srun
     double yrsrun;  // (mm / year)

     // Annual sum of vsm
     double yrvsm;


   private:

/* *************************************************************
		 Private Functions
************************************************************* */

     void percol( const double& rain );

     double rrunoff( const double& rgrndh2o );

     double srunoff( const double& elev,
                     const double& tair,
                     const double& prevtair,
                     const double& prev2tair,
		             const double& sgrndh2o );


/* **************************************************************
		 Private Variables
************************************************************** */

     // Available inorganic nitrogen (g N / sq. meter)
     double availn;

     // Available water (mm)
     double avlh2o;

     // Available water capacity (mm)
     double awcapmm;

     // Dissolved Organic Nitrogen leaching
     double lchdon;

     // Dissolved Inorganic Nitrogen leaching
     double lchdin;
     
     // Monthly drainage (mm / month)
     double drainage;

     // Monthly estimated actual Evapotranspiration (mm / month)
     double eet;
     
     // Volume of water at field capacity (mm)
     double fldcap;

     // Soil dryness function, similar in shape to soil water potential
     double gm;

      // Water yield (mm / month)
     double h2oyld;

     // Initial Estimated Actual Evapotranspiration (mm / month)
     double ineet;

     // Relative hydraulic conductivity through soil profile
     double kh2o;

    // DOC leaching flux
     double lchdoc;

     // Mean annual volumetric soil moisture (%)
     double meanvsm;

     // Soil moisture (mm)
     double moist;

     double ndays[CYCLE];

     // Total nitrogen input to soils
     double ninput;    // (g N / (sq. meter * month))

     // Total nitrogen lost from soils
     double nlost;     // (g N / (sq. meter * month))

      // Reactive soil organic matter
     Biomass org;      //  (g C or g N / sq. meter)

     // Soil moisture as %field capacity
     double pcfc;

     // Percent clay in soil texture
     double pctclay;

     // Soil moisture as %total pore space
     double pctp;

     // Percent sand in soil texture
     double pctsand;

     // Percent silt in soil texture
     double pctsilt;

     // Previous month's snow pack
     double prevspack;

     // Proportion silt and clay in soil texture
     double psiplusc;

     // Rain runoff storage (mm / month)
     double rgrndh2o;

     // Rain percolation (excess)
     double rperc;   // (mm / month)

     // Rain Runoff (mm / month)
     double rrun;

     // Snowmelt runoff storage (mm)
     double sgrndh2o;

     // Snow melt infiltration (mm / month)
     double snowinf;

     // Snowpack (mm)
     double snowpack;

     // Snow melt percolation (excess)
     double sperc;   // (mm / month)

     // Nonsymbiotic N fix
     double soninp;

     // Snow runoff
     double srun;    // (mm / month)

     // Surface runoff (mm / month)
     double surfrun;
     
     // Soil water potential (MPa)
     double swp; 

     // Soil texture (categorical data)
     int text;

     // volume of total pore space (mm)
     double totpor;

     // Volumetric soil moisture (as %rooting depth)
     double vsm;

     // Volume of water at wilting point (mm)
     double wiltpt;

     // wetland soil type designation (categorical data)
     int wsoil;


/* *************************************************************
		 Private Parameters
************************************************************* */

     // Field capacity (%soil volume)

     double fldcapa;
     double fldcapb;
     double pcfldcap;


     // Proportion of available nitrogen lost from soil
     //   (g N / (square meter))

     double nloss[MAXCMNT];
     double denitr[MAXCMNT];


     // Porosity of soil (%soil volume)

     double pctpor;
     double pctpora;
     double pctporb;


     // Effective rooting depth (m)

     double minrootz[MAXCMNT];
     double z0[MAXCMNT];
     double rootz;
     double rootza[MAXCMNT];
     double rootzb[MAXCMNT];
     double rootzc[MAXCMNT];
     double effwiltpt;

     double reet;

     // Wilting point (%soil volume)

     double pcwiltpt;
     double wiltpta;
     double wiltptb;
     
     // Soil Surface resistance (s m^-1)
     
     double rssmin;
     double rssslope;

};

#endif
