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
#include "TZDDData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TZDDData)


//////////////////////////////////////////////////////////////////////
TZDDData::TZDDData() {
}



//////////////////////////////////////////////////////////////////////
TZDDData::~TZDDData() {
}



//////////////////////////////////////////////////////////////////////
void TZDDData::Clear() {
  // Drift_Chambers
  fZDD_Drift_DetectorNbr.clear();
  fZDD_DriftTime.clear();
  // Energy
  fZDD_E_DetectorNbr.clear();
  fZDD_Energy.clear();
  // Time
  fZDD_T_DetectorNbr.clear();
  fZDD_Time.clear();

  // IC Energy
  fZDD_E_IC_Nbr.clear();
  fZDD_IC_Energy.clear();
  // IC Time
  fZDD_T_IC_Nbr.clear();
  fZDD_IC_Time.clear();
  
  // Plastic Energy
  fZDD_E_Plastic_Nbr.clear();
  fZDD_Plastic_Energy.clear();
  // Plastic Time
  fZDD_T_Plastic_Nbr.clear();
  fZDD_Plastic_Time.clear();
}



//////////////////////////////////////////////////////////////////////
void TZDDData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TZDDData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // Energy
  size_t mysize = fZDD_E_DetectorNbr.size();
  cout << "ZDD_E_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fZDD_E_DetectorNbr[i]
         << " Energy: " << fZDD_Energy[i];
  }
  
  // Time
  mysize = fZDD_T_DetectorNbr.size();
  cout << "ZDD_T_Mult: " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "DetNbr: " << fZDD_T_DetectorNbr[i]
         << " Time: " << fZDD_Time[i];
  }
}
