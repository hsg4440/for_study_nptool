#ifndef __TACDATA__
#define __TACDATA__
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TTACData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
    vector<UShort_t>   fTAC_Channel;
    vector<UShort_t>   fTAC_Nbr;
    vector<ULong64_t>  fTAC_TS;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TTACData();
    ~TTACData();
    

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
    inline void SetTAC(const UShort_t& TACNbr,const UShort_t& Channel, const ULong64_t& TS){
      fTAC_Nbr.push_back(TACNbr);
      fTAC_Channel.push_back(Channel);
      fTAC_TS.push_back(TS);
    };//!



    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetTAC_Mult() const
      {return fTAC_Channel.size();}
    inline UShort_t GetTAC_Channel(const UShort_t &i) const 
      {return fTAC_Channel[i];}//!
    inline UShort_t GetTAC_Nbr(const UShort_t &i) const 
      {return fTAC_Nbr[i];}//!
    inline ULong64_t GetTAC_TS(const ULong64_t &i) const 
      {return fTAC_TS[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TTACData,1)  // TACData structure
};

#endif
