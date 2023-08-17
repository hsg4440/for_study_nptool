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
#include <vector>

class TExogamData : public TObject {
 private:
  std::vector<UShort_t> fExoE;
  std::vector<UShort_t> fExoE_CrystalNbr;
  std::vector<ULong64_t> fExoE_TS;

  std::vector<UShort_t> fExoHG; // Same as Energy but with High Gain (for Low Energy events) 
  std::vector<UShort_t> fExoHG_CrystalNbr;
  std::vector<ULong64_t> fExoHG_TS;

  std::vector<UShort_t> fExoTDC; // Internal TDC of EXOGAM
  std::vector<UShort_t> fExoTDC_CrystalNbr;
  std::vector<ULong64_t> fExoTDC_TS;

  std::vector<UShort_t> fExoOuter;
  std::vector<UShort_t> fExoOuter_SubCrystalNbr; 

  std::vector<UShort_t> fExoBGO;
  std::vector<UShort_t> fExoBGO_CrystalNbr;

  std::vector<UShort_t> fExoCsI;
  std::vector<UShort_t> fExoCsI_CrystalNbr;

 public:
  TExogamData();
  virtual ~TExogamData();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  /////////////////////           SETTERS           ////////////////////////
  inline void SetExo(const UShort_t& Energy, const UShort_t& DetNumb, const ULong64_t& TimeStamp) {
    fExoE.push_back(Energy);
    fExoE_CrystalNbr.push_back(DetNumb);
    fExoE_TS.push_back(TimeStamp);
  }
  inline void SetExoHG(const UShort_t& Energy, const UShort_t& DetNumb, const ULong64_t& TimeStamp) {
    fExoHG.push_back(Energy);
    fExoHG_CrystalNbr.push_back(DetNumb);
    fExoHG_TS.push_back(TimeStamp);
  }
  inline void SetExoDelta(const UShort_t& Time, const UShort_t& DetNumb, const ULong64_t& TimeStamp) {
    fExoTDC.push_back(Time);
    fExoTDC_CrystalNbr.push_back(DetNumb);
    fExoTDC_TS.push_back(TimeStamp);
  }
  inline void SetExoOuter(const UShort_t& Energy, const UShort_t& OutersNumb) {
    fExoOuter.push_back(Energy);
    fExoOuter_SubCrystalNbr.push_back(OutersNumb);
  }
  inline void SetExoBGO(const UShort_t& Energy, const UShort_t& BGONumb) {
    fExoBGO.push_back(Energy);
    fExoBGO_CrystalNbr.push_back(BGONumb);
  }
  inline void SetExoCsI(const UShort_t& Energy, const UShort_t& CsINumb) {
    fExoCsI.push_back(Energy);
    fExoCsI_CrystalNbr.push_back(CsINumb);
  }
  /////////////////////           GETTERS           ////////////////////////
  inline UShort_t GetExoE(UShort_t& i) { return fExoE[i]; }
  inline UShort_t GetExoE_CrystalNbr(UShort_t& i) { return fExoE_CrystalNbr[i]; }
  inline ULong64_t GetExoE_TS(UShort_t& i) { return fExoE_TS[i]; }
  inline UShort_t GetExoHG(UShort_t& i) { return fExoHG[i]; }
  inline UShort_t GetExoHG_CrystalNbr(UShort_t& i) { return fExoHG_CrystalNbr[i]; }
  inline ULong64_t GetExoHG_TS(UShort_t& i) { return fExoHG_TS[i]; }
  inline UShort_t GetExoTDC(UShort_t& i) { return fExoTDC[i]; }
  inline UShort_t GetExoTDC_CrystalNbr(UShort_t& i) { return fExoTDC_CrystalNbr[i]; }
  inline ULong64_t GetExoTDC_TS(UShort_t& i) { return fExoTDC_TS[i]; }
  inline UShort_t GetExoOuter(UShort_t& i) { return fExoOuter[i]; }
  inline UShort_t GetExoOuter_SubCrystalNbr(UShort_t& i) { return fExoOuter_SubCrystalNbr[i]; }
  inline UShort_t GetExoBGO(UShort_t& i) { return fExoBGO[i]; }
  inline UShort_t GetExoBGO_CrystalNbr(UShort_t& i) { return fExoBGO_CrystalNbr[i]; }
  inline UShort_t GetExoCsI(UShort_t& i) { return fExoCsI[i]; }
  inline UShort_t GetExoCsI_CrystalNbr(UShort_t& i) { return fExoCsI_CrystalNbr[i]; }

  ClassDef(TExogamData, 1) // ExogamData structure
};

#endif
