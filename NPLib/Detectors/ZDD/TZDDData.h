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
    vector<ULong64_t>   fZDD_IC_TS;

    // Plastic
    vector<UShort_t>   fZDD_PM_E;
    vector<UShort_t>   fZDD_PM_N;
    vector<ULong64_t>   fZDD_PM_TS;
    
    // DC
    vector<UShort_t>   fZDD_DC_E;
    vector<UShort_t>   fZDD_DC_N;
    vector<ULong64_t>   fZDD_DC_TS;
    
    // EXOZDD
    vector<UShort_t>   fZDD_EXO_E;
    vector<UShort_t>   fZDD_EXO_N;
    vector<ULong64_t>   fZDD_EXO_TS;


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
    inline void SetZDDIC(const UShort_t& DetNumb,const UShort_t& Energy, const ULong64_t& TimeStamp){
      fZDD_IC_E.push_back(Energy);
      fZDD_IC_N.push_back(DetNumb);
      fZDD_IC_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDPM(const UShort_t& DetNumb,const UShort_t& Energy, const ULong64_t& TimeStamp){
      fZDD_PM_E.push_back(Energy);
      fZDD_PM_N.push_back(DetNumb);
      fZDD_PM_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDDC(const UShort_t& DetNumb,const UShort_t& Energy, const ULong64_t& TimeStamp){
      fZDD_DC_E.push_back(Energy);
      fZDD_DC_N.push_back(DetNumb);
      fZDD_DC_TS.push_back(TimeStamp);
    };//!
    inline void SetZDDEXO(const UShort_t& DetNumb,const UShort_t& Energy, const ULong64_t& TimeStamp){
      fZDD_EXO_E.push_back(Energy);
      fZDD_EXO_N.push_back(DetNumb);
      fZDD_EXO_TS.push_back(TimeStamp);
    };//!
    
    //////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetZDD_ICE(const UShort_t& i)const  { return fZDD_IC_E[i]; }
    inline UShort_t GetZDD_ICN(const UShort_t& i)const  { return fZDD_IC_N[i]; }
    inline ULong64_t GetZDD_ICTS(const UShort_t& i)const  { return fZDD_IC_TS[i]; }
    inline UShort_t GetZDD_PME(const UShort_t& i)const  { return fZDD_PM_E[i]; }
    inline UShort_t GetZDD_PMN(const UShort_t& i)const  { return fZDD_PM_N[i]; }
    inline ULong64_t GetZDD_PMTS(const UShort_t& i)const  { return fZDD_PM_TS[i]; }
    inline UShort_t GetZDD_DCE(const UShort_t& i)const  { return fZDD_DC_E[i]; }
    inline UShort_t GetZDD_DCN(const UShort_t& i)const  { return fZDD_DC_N[i]; }
    inline ULong64_t GetZDD_DCTS(const UShort_t& i)const  { return fZDD_DC_TS[i]; }
    inline UShort_t GetZDD_EXOE(const UShort_t& i)const  { return fZDD_EXO_E[i]; }
    inline UShort_t GetZDD_EXON(const UShort_t& i)const  { return fZDD_EXO_N[i]; }
    inline ULong64_t GetZDD_EXOTS(const UShort_t& i)const  { return fZDD_EXO_TS[i]; }

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TZDDData,1)  // ZDDData structure
};

#endif