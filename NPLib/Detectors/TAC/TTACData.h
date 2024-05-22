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
  vector<unsigned int> fTAC_Time;
  vector<unsigned int> fTAC_N;
  vector<std::string> fTAC_Name;
  vector<unsigned long long> fTAC_TS;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
 public:
  TTACData();
  ~TTACData();

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
 public:
  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
 public:
  //////////////////////    SETTERS    ////////////////////////
  // Channel -> N
  inline void SetTAC(const unsigned int& Channel, const unsigned int& Time, const unsigned long long& TS,
                     const std::string& Name) {
    fTAC_N.push_back(Channel);
    fTAC_Name.push_back(Name);
    fTAC_Time.push_back(Time);
    fTAC_TS.push_back(TS);
  }; //!

  //////////////////////    GETTERS    ////////////////////////
  // Energy
  inline unsigned int GetTAC_Mult() const { return fTAC_Time.size(); }
  inline unsigned int GetTAC_Time(const unsigned int& i) const { return fTAC_Time[i]; } //!
  inline unsigned int GetTAC_N(const unsigned int& i) const { return fTAC_N[i]; }        //!
  inline std::string GetTAC_Name(const unsigned int& i) const { return fTAC_Name[i]; }        //!
  inline unsigned long long GetTAC_TS(const unsigned long long& i) const { return fTAC_TS[i]; }    //!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TTACData, 1) // TACData structure
};

#endif
