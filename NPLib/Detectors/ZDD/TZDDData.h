#ifndef __ZDDDATA__
#define __ZDDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <iostream>
using namespace std;

// ROOT
#include "TObject.h"


class TZDDData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    // IC
    vector<UShort_t>   fZDD_IC_E;
    vector<UShort_t>   fZDD_IC_N;
    vector<UShort_t>   fZDD_IC_TS;

    // Plastic
    vector<UShort_t>   fZDD_PM_E;
    vector<UShort_t>   fZDD_PM_N;
    vector<UShort_t>   fZDD_PM_TS;
    
    // DC
    vector<UShort_t>   fZDD_DC_E;
    vector<UShort_t>   fZDD_DC_N;
    vector<UShort_t>   fZDD_DC_TS;
    
    // EXOZDD
    vector<UShort_t>   fZDD_EXO_E;
    vector<UShort_t>   fZDD_EXO_N;
    vector<UShort_t>   fZDD_EXO_TS;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TZDDData();
    ~TZDDData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // Energy
    inline void SetZDDIC(const UShort_t& DetNumb,const UShort_t& Energy, const UShort_t& TimeStamp){
      fZDD_IC_E.push_back(Energy);
      fZDD_IC_N.push_back(DetNumb);
      fZDD_IC_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDPM(const UShort_t& DetNumb,const UShort_t& Energy, const UShort_t& TimeStamp){
      fZDD_PM_E.push_back(Energy);
      fZDD_PM_N.push_back(DetNumb);
      fZDD_PM_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDDC(const UShort_t& DetNumb,const UShort_t& Energy, const UShort_t& TimeStamp){
      fZDD_DC_E.push_back(Energy);
      fZDD_DC_N.push_back(DetNumb);
      fZDD_DC_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDEXO(const UShort_t& DetNumb,const UShort_t& Energy, const UShort_t& TimeStamp){
      fZDD_EXO_E.push_back(Energy);
      fZDD_EXO_N.push_back(DetNumb);
      fZDD_EXO_TS.push_back(TimeStamp);
    };//!
    
    //////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetZDD_ICE(UShort_t& i) { return fZDD_IC_E[i]; }
    inline UShort_t GetZDD_ICN(UShort_t& i) { return fZDD_IC_N[i]; }
    inline UShort_t GetZDD_ICTS(UShort_t& i) { return fZDD_IC_TS[i]; }
    inline UShort_t GetZDD_PME(UShort_t& i) { return fZDD_PM_E[i]; }
    inline UShort_t GetZDD_PMN(UShort_t& i) { return fZDD_PM_N[i]; }
    inline UShort_t GetZDD_PMTS(UShort_t& i) { return fZDD_PM_TS[i]; }
    inline UShort_t GetZDD_DCE(UShort_t& i) { return fZDD_DC_E[i]; }
    inline UShort_t GetZDD_DCN(UShort_t& i) { return fZDD_DC_N[i]; }
    inline UShort_t GetZDD_DCTS(UShort_t& i) { return fZDD_DC_TS[i]; }
    inline UShort_t GetZDD_EXOE(UShort_t& i) { return fZDD_EXO_E[i]; }
    inline UShort_t GetZDD_EXON(UShort_t& i) { return fZDD_EXO_N[i]; }
    inline UShort_t GetZDD_EXOTS(UShort_t& i) { return fZDD_EXO_TS[i]; }

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TZDDData,1)  // ZDDData structure
};

#endif