#ifndef __EXOGAMDATA__
#define __EXOGAMDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2009                                               *
 * Last update    : july 2019                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: Added vectors for real energy/time values (double) (T.Goigoux CEA) *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include "TObject.h"
#include <map>
#include <vector>

class TExogamData : public TObject {

 public:
  TExogamData();
  virtual ~TExogamData();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

 public:
  std::vector<unsigned int> fExo_Flange;
  std::vector<unsigned int> fExo_Crystal;
  std::vector<unsigned int> fExo_E;
  std::vector<unsigned int> fExo_E_HG; // High gain x20
  std::vector<unsigned long long> fExo_TS;
  std::vector<unsigned int> fExo_TDC;
  std::vector<unsigned int> fExo_BGO;
  std::vector<unsigned int> fExo_CsI;
  std::vector<unsigned int> fExo_Outer1;
  std::vector<unsigned int> fExo_Outer2;
  std::vector<unsigned int> fExo_Outer3;
  std::vector<unsigned int> fExo_Outer4;


  /////////////////////           SETTERS           ////////////////////////
  inline void SetExoFlange(unsigned int& Flange) { fExo_Flange.push_back(Flange); }
  inline void SetExoCrystal(unsigned int& Crystal) { fExo_Crystal.push_back(Crystal); }
  inline void SetExoE(unsigned int& Energy) { fExo_E.push_back(Energy); }
  inline void SetExoEHG(unsigned int& Energy) { fExo_E_HG.push_back(Energy); }
  inline void SetExoTS(unsigned long long& TS) { fExo_TS.push_back(TS); }
  inline void SetExoTDC(unsigned int& TDC) { fExo_TDC.push_back(TDC); }
  inline void SetExoBGO(unsigned int& BGO) { fExo_BGO.push_back(BGO); }
  inline void SetExoCsI(unsigned int& CsI) { fExo_CsI.push_back(CsI); }
  inline void SetExoOuter1(unsigned int& Outer1) { fExo_Outer1.push_back(Outer1); }
  inline void SetExoOuter2(unsigned int& Outer2) { fExo_Outer2.push_back(Outer2); }
  inline void SetExoOuter3(unsigned int& Outer3) { fExo_Outer3.push_back(Outer3); }
  inline void SetExoOuter4(unsigned int& Outer4) { fExo_Outer4.push_back(Outer4); }

  /////////////////////           GETTERS           ////////////////////////
  inline unsigned int GetExoFlange(unsigned int& i) { return fExo_Flange[i]; }
  inline unsigned int GetExoCrystal(unsigned int& i) { return fExo_Crystal[i]; }
  inline unsigned int GetExoE(unsigned int& i) { return fExo_E[i]; }
  inline unsigned int GetExoEHG(unsigned int& i) { return fExo_E_HG[i]; }
  inline unsigned long long GetExoTS(unsigned int& i) { return fExo_TS[i]; }
  inline unsigned int GetExoTDC(unsigned int& i) { return fExo_TDC[i]; }
  inline unsigned int GetExoBGO(unsigned int& i) { return fExo_BGO[i]; }
  inline unsigned int GetExoCsI(unsigned int& i) { return fExo_CsI[i]; }
  inline unsigned int GetExoOuter1(unsigned int& i) { return fExo_Outer1[i]; }
  inline unsigned int GetExoOuter2(unsigned int& i) { return fExo_Outer2[i]; }
  inline unsigned int GetExoOuter3(unsigned int& i) { return fExo_Outer3[i]; }
  inline unsigned int GetExoOuter4(unsigned int& i) { return fExo_Outer4[i]; }

  ClassDef(TExogamData, 1) // ExogamData structure
};

#endif
