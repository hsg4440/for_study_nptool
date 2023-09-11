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
  fExoE.clear();
  fExoE_CrystalNbr.clear();
  fExoE_TS.clear();
  fExoHG.clear(); 
  fExoHG_CrystalNbr.clear();
  fExoHG_TS.clear();
  fExoTDC.clear();
  fExoTDC_CrystalNbr.clear();
  fExoTDC_TS.clear();
  fExoOuter.clear();
  fExoOuter_SubCrystalNbr.clear(); 
  fExoBGO.clear();
  fExoBGO_CrystalNbr.clear();
  fExoCsI.clear();
  fExoCsI_CrystalNbr.clear();
}



void TExogamData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event XXXXXXXXXXXXXXXXX" << endl;
   
   cout << "Inner6MV Mult = " << fExoE.size() << endl;
   for (UShort_t i = 0; i < fExoE.size(); i++) {
      cout << "Energy: " << fExoE[i] << " Cristal Numb: " << fExoE_CrystalNbr[i] << " TimeStamp: " << fExoE_TS[i] << endl;
   }
   cout << "Inner20MV Mult = " << fExoHG.size() << endl;
   for (UShort_t i = 0; i < fExoHG.size(); i++) {
      cout << "Energy: " << fExoHG[i] << " Cristal Numb: " << fExoHG_CrystalNbr[i] << " TimeStamp: " << fExoHG_TS[i] << endl;
   }
   cout << "OutersV Mult = " << fExoOuter.size() << endl;
   for (UShort_t i = 0; i < fExoOuter.size(); i++) {
      cout << "Energy: " << fExoOuter[i] << " Cristal Numb: " << fExoOuter_SubCrystalNbr[i] << endl;
   }
   cout << "DeltaTV Mult = " << fExoTDC.size() << endl;
   for (UShort_t i = 0; i < fExoTDC.size(); i++) {
      cout << "Energy: " << fExoTDC[i] << " Cristal Numb: " << fExoTDC_CrystalNbr[i] << " TimeStamp: " << fExoTDC_TS[i] << endl;
   }
   cout << "BGOV Mult = " << fExoBGO.size() << endl;
   for (UShort_t i = 0; i < fExoBGO.size(); i++) {
      cout << "Energy: " << fExoBGO[i] << " Cristal Numb: " << fExoBGO_CrystalNbr[i] << endl;
   }
   cout << "CsIV Mult = " << fExoCsI.size() << endl;
   for (UShort_t i = 0; i < fExoCsI.size(); i++) {
      cout << "Energy: " << fExoCsI[i] << " Cristal Numb: " << fExoCsI_CrystalNbr[i] << endl;
   }
}
