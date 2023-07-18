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
  std::vector<UShort_t> fEXO_E;
  std::vector<UShort_t> fEXO_E_CrystalNbr;
  std::vector<UShort_t> fEXO_E_TS;

  std::vector<UShort_t> fEXO_HG; // Same as Energy but with High Gain (for Low Energy events) 
  std::vector<UShort_t> fEXO_HG_CrystalNbr;
  std::vector<UShort_t> fEXO_HG_TS;

  std::vector<UShort_t> fEXO_TDC; // Internal TDC of EXOGAM
  std::vector<UShort_t> fEXO_TDC_CrystalNbr;
  std::vector<UShort_t> fEXO_TDC_TS;

  std::vector<UShort_t> fEXO_Outer;
  std::vector<UShort_t> fEXO_Outer_SubCrystalNbr; 

  std::vector<UShort_t> fEXO_BGO;
  std::vector<UShort_t> fEXO_BGO_CrystalNbr;

  std::vector<UShort_t> fEXO_CSI;
  std::vector<UShort_t> fEXO_CSI_CrystalNbr;

 public:
  TExogamData();
  virtual ~TExogamData();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  /////////////////////           SETTERS           ////////////////////////
  inline void SetInner6MV(UShort_t& Energy, UShort_t& DetNumb, UShort_t& TimeStamp) {
    fEXO_E.push_back(Energy);
    fEXO_E_CrystalNbr.push_back(DetNumb);
    fEXO_E_TS.push_back(TimeStamp);
  }
  inline void SetInner20MV(UShort_t& Energy, UShort_t& DetNumb, UShort_t& TimeStamp) {
    fEXO_HG.push_back(Energy);
    fEXO_HG_CrystalNbr.push_back(DetNumb);
    fEXO_HG_TS.push_back(TimeStamp);
  }
  inline void SetDeltaTV(UShort_t& Time, UShort_t& DetNumb, UShort_t& TimeStamp) {
    fEXO_TDC.push_back(Time);
    fEXO_TDC_CrystalNbr.push_back(DetNumb);
    fEXO_TDC_TS.push_back(TimeStamp);
  }
  inline void SetOutersV(UShort_t& Energy, UShort_t& OutersNumb) {
    fEXO_Outer.push_back(Energy);
    fEXO_Outer_SubCrystalNbr.push_back(OutersNumb);
  }
  inline void SetBGOV(UShort_t& Energy, UShort_t& BGONumb) {
    fEXO_BGO.push_back(Energy);
    fEXO_BGO_CrystalNbr.push_back(BGONumb);
  }
  inline void SetCSIV(UShort_t& Energy, UShort_t& CSINumb) {
    fEXO_CSI.push_back(Energy);
    fEXO_CSI_CrystalNbr.push_back(CSINumb);
  }
  /////////////////////           GETTERS           ////////////////////////
  inline UShort_t GetEXO_E(UShort_t& i) { return fEXO_E[i]; }
  inline UShort_t GetEXO_E_CrystalNbr(UShort_t& i) { return fEXO_E_CrystalNbr[i]; }
  inline UShort_t GetEXO_E_TS(UShort_t& i) { return fEXO_E_TS[i]; }
  inline UShort_t GetEXO_HG(UShort_t& i) { return fEXO_HG[i]; }
  inline UShort_t GetEXO_HG_CrystalNbr(UShort_t& i) { return fEXO_HG_CrystalNbr[i]; }
  inline UShort_t GetEXO_HG_TS(UShort_t& i) { return fEXO_HG_TS[i]; }
  inline UShort_t GetEXO_TDC(UShort_t& i) { return fEXO_TDC[i]; }
  inline UShort_t GetEXO_TDC_CrystalNbr(UShort_t& i) { return fEXO_TDC_CrystalNbr[i]; }
  inline UShort_t GetEXO_TDC_TS(UShort_t& i) { return fEXO_TDC_TS[i]; }
  inline UShort_t GetEXO_Outer(UShort_t& i) { return fEXO_Outer[i]; }
  inline UShort_t GetEXO_Outer_SubCrystalNbr(UShort_t& i) { return fEXO_Outer_SubCrystalNbr[i]; }
  inline UShort_t GetEXO_BGO(UShort_t& i) { return fEXO_BGO[i]; }
  inline UShort_t GetEXO_BGO_CrystalNbr(UShort_t& i) { return fEXO_BGO_CrystalNbr[i]; }
  inline UShort_t GetEXO_CSI(UShort_t& i) { return fEXO_CSI[i]; }
  inline UShort_t GetEXO_CSI_CrystalNbr(UShort_t& i) { return fEXO_CSI_CrystalNbr[i]; }

  ClassDef(TExogamData, 1) // ExogamData structure
};

#endif
