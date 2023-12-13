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
  inline void SetExo(unsigned int& Crystal, unsigned int& Energy,
  unsigned int& Energy_HG, unsigned long long& TS, unsigned int& TDC,
  unsigned int& BGO, unsigned int& CsI, unsigned int& Outer1,
  unsigned int& Outer2, unsigned int& Outer3,unsigned int& Outer4) { 
  fExo_Crystal.push_back(Crystal);
  fExo_E.push_back(Energy);
  fExo_E_HG.push_back(Energy_HG);
  fExo_TS.push_back(TS);
  fExo_TDC.push_back(TDC);
  fExo_BGO.push_back(BGO);
  fExo_CsI.push_back(CsI);
  fExo_Outer1.push_back(Outer1);
  fExo_Outer2.push_back(Outer2);
  fExo_Outer3.push_back(Outer3);
  fExo_Outer4.push_back(Outer4);
  }
  /////////////////////           GETTERS           ////////////////////////
  inline unsigned int GetExoMult(unsigned int& i) { return fExo_Crystal[i]; }
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
