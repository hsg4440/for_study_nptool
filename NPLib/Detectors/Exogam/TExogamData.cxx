/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2009                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
using namespace std;

#include "TExogamData.h"


ClassImp(TExogamData)

TExogamData::TExogamData()
{
   // Default constructor
   Clear();
}



TExogamData::~TExogamData()
{
}



void TExogamData::Clear()
{
  fEXO_E.clear();
  fEXO_E_CrystalNbr.clear();
  fEXO_E_TS.clear();
  fEXO_HG.clear(); 
  fEXO_HG_CrystalNbr.clear();
  fEXO_HG_TS.clear();
  fEXO_TDC.clear();
  fEXO_TDC_CrystalNbr.clear();
  fEXO_TDC_TS.clear();
  fEXO_Outer.clear();
  fEXO_Outer_SubCrystalNbr.clear(); 
  fEXO_BGO.clear();
  fEXO_BGO_CrystalNbr.clear();
  fEXO_CSI.clear();
  fEXO_CSI_CrystalNbr.clear();
}



void TExogamData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event XXXXXXXXXXXXXXXXX" << endl;
   
   cout << "Inner6MV Mult = " << fEXO_E.size() << endl;
   for (UShort_t i = 0; i < fEXO_E.size(); i++) {
      cout << "Energy: " << fEXO_E[i] << " Cristal Numb: " << fEXO_E_CrystalNbr[i] << " TimeStamp: " << fEXO_E_TS[i] << endl;
   }
   cout << "Inner20MV Mult = " << fEXO_HG.size() << endl;
   for (UShort_t i = 0; i < fEXO_HG.size(); i++) {
      cout << "Energy: " << fEXO_HG[i] << " Cristal Numb: " << fEXO_HG_CrystalNbr[i] << " TimeStamp: " << fEXO_HG_TS[i] << endl;
   }
   cout << "OutersV Mult = " << fEXO_Outer.size() << endl;
   for (UShort_t i = 0; i < fEXO_Outer.size(); i++) {
      cout << "Energy: " << fEXO_Outer[i] << " Cristal Numb: " << fEXO_Outer_SubCrystalNbr[i] << endl;
   }
   cout << "DeltaTV Mult = " << fEXO_TDC.size() << endl;
   for (UShort_t i = 0; i < fEXO_TDC.size(); i++) {
      cout << "Energy: " << fEXO_TDC[i] << " Cristal Numb: " << fEXO_TDC_CrystalNbr[i] << " TimeStamp: " << fEXO_TDC_TS[i] << endl;
   }
   cout << "BGOV Mult = " << fEXO_BGO.size() << endl;
   for (UShort_t i = 0; i < fEXO_BGO.size(); i++) {
      cout << "Energy: " << fEXO_BGO[i] << " Cristal Numb: " << fEXO_BGO_CrystalNbr[i] << endl;
   }
   cout << "CSIV Mult = " << fEXO_CSI.size() << endl;
   for (UShort_t i = 0; i < fEXO_CSI.size(); i++) {
      cout << "Energy: " << fEXO_CSI[i] << " Cristal Numb: " << fEXO_CSI_CrystalNbr[i] << endl;
   }
}
