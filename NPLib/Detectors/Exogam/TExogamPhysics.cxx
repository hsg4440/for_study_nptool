/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sandra GIRON  contact address: giron@ipno.in2p3.fr       *
 *                  Benjamin LE CROM		   lecrom@ipno.in2p3.fr              *
 * Creation Date  : march 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold exogam treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TExogamPhysics.h"
using namespace EXOGAM_LOCAL;

//	STL
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <functional>

//	NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPVDetector.h"
#include "RootInput.h"
#include "RootOutput.h"
//	ROOT
#include "TChain.h"

///////////////////////////////////////////////////////////////////////////
ClassImp(TExogamPhysics)
    ///////////////////////////////////////////////////////////////////////////
    TExogamPhysics::TExogamPhysics() {
  EventMultiplicity = 0;
  ECC_Multiplicity = 0;
  GOCCE_Multiplicity = 0;
  NumberOfHitClover = 0;
  NumberOfHitCristal = 0;
  m_Spectra = NULL;
  NumberOfClover = 0;
  m_EXO_E_RAW_Threshold = 0;
  m_EXO_E_Threshold = 0;
  m_EXO_EHG_RAW_Threshold = 0;
  m_EXO_TDC_RAW_Threshold = 0;
  m_EXO_TDC_RAW_Threshold = 0;

  m_PreTreatedData = new TExogamData;
  m_EventData = new TExogamData;
  m_EventPhysics = this;
  NumberOfClover = 0;
  CloverMult = 0;
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::PreTreat() {
  // Clearing PreTreat TExogamData
  ClearPreTreatedData();
  // Clearing local variables for pretreat

  //E
  m_EXO_Mult = m_EventData->GetExoMult();

  for (unsigned int i = 0; i < m_EXO_Mult; ++i) {
  ResetPreTreatVariable();
    // std::cout << "test1 " << EXO_E << std::endl;
    if (m_EventData->GetExoE(i) > m_EXO_E_RAW_Threshold)
      EXO_E = fEXO_E(m_EventData, i);
    // std::cout << "test2 " << EXO_E << std::endl;

    if (m_EventData->GetExoEHG(i) > m_EXO_EHG_RAW_Threshold)
      EXO_EHG = fEXO_EHG(m_EventData, i);
    
    if (m_EventData->GetExoTDC(i) > m_EXO_TDC_RAW_Threshold)
      EXO_TDC = fEXO_T(m_EventData, i);
    
    EXO_Outer1 = fEXO_Outer(m_EventData, i, 0);
    EXO_Outer2 = fEXO_Outer(m_EventData, i, 1);
    EXO_Outer3 = fEXO_Outer(m_EventData, i, 2);
    EXO_Outer4 = fEXO_Outer(m_EventData, i, 3);

    if(EXO_E > m_EXO_E_Threshold){
      m_PreTreatedData->SetExo(m_EventData->GetExoCrystal(i), EXO_E*1000,
      EXO_EHG*1000, m_EventData->GetExoTS(i), EXO_TDC, 
      m_EventData->GetExoBGO(i), m_EventData->GetExoCsI(i), EXO_Outer1,
      EXO_Outer2, EXO_Outer3, EXO_Outer4);
    } 
  }
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ResetPreTreatVariable(){
  EXO_E = -1000;
  EXO_EHG = -1000;
  EXO_TDC = -1000;
  EXO_Outer1 = -1000;
  EXO_Outer2 = -1000;
  EXO_Outer3 = -1000;
  EXO_Outer4 = -1000;
}

void TExogamPhysics::BuildPhysicalEvent() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
  PreTreat();

  for(unsigned int i = 0; i < m_PreTreatedData->GetExoMult(); i++){
    
    mean_free_path = ComputeMeanFreePath(m_PreTreatedData->GetExoE(i));
    flange_nbr = MapCrystalFlangeCLover[m_PreTreatedData->GetExoCrystal(i)].first;
    crystal_nbr = MapCrystalFlangeCLover[m_PreTreatedData->GetExoCrystal(i)].second;

    EnergyDoppler = 0; 
    double MaxOuter = 0;
    unsigned int MaxOuter_Id;
    if(m_PreTreatedData->GetExoOuter1(i) > MaxOuter){
      MaxOuter = m_PreTreatedData->GetExoOuter1(i);
      MaxOuter_Id = 0;
    }
    if(m_PreTreatedData->GetExoOuter2(i) > MaxOuter){
      MaxOuter = m_PreTreatedData->GetExoOuter2(i);
      MaxOuter_Id = 1;
    }
    if(m_PreTreatedData->GetExoOuter3(i) > MaxOuter){
      MaxOuter = m_PreTreatedData->GetExoOuter3(i);
      MaxOuter_Id = 2;
    }
    if(m_PreTreatedData->GetExoOuter4(i) > MaxOuter){
      MaxOuter = m_PreTreatedData->GetExoOuter4(i);
      MaxOuter_Id = 3;
    }
    
    if(MaxOuter > 0){
      Exogam_struc = Ask_For_Angles(flange_nbr, mean_free_path);
      double Theta_seg = Exogam_struc.Theta_Crystal_Seg[crystal_nbr][MaxOuter_Id];
      double Phi_seg = Exogam_struc.Phi_Crystal_Seg[crystal_nbr][MaxOuter_Id];
      EnergyDoppler = Doppler_Correction(Theta_seg,Phi_seg,0,0,Beta,m_PreTreatedData->GetExoE(i));
    }
    E_Cal.push_back(m_PreTreatedData->GetExoE(i));
    Flange_N.push_back(flange_nbr);
    Crystal_N.push_back(crystal_nbr);
    E_Doppler.push_back(EnergyDoppler);
  }
/*
  if(PreTreatedData -> GetECCEMult() != PreTreatedData -> GetECCTMult()) cout << PreTreatedData -> GetECCEMult() << " "
  <<  PreTreatedData -> GetECCTMult() << endl;


  for(unsigned int i = 0 ; i < PreTreatedData -> GetECCEMult(); i++) {

    // cout << i << " " << cristal_E << endl;
    // if(PreTreatedData->GetECCTTime(i) > 0)
      {
    ECC_E.push_back(PreTreatedData->GetECCEEnergy(i));
    ECC_T.push_back(PreTreatedData->GetECCTTime(i));
    ECC_CloverNumber.push_back(PreTreatedData->GetECCEClover(i));
    ECC_CristalNumber.push_back(PreTreatedData->GetECCECristal(i));

    //    cout << "BuildPhys " << PreTreatedData->GetECCEClover(i) << " " <<  PreTreatedData->GetECCECristal(i)<< " " <<
  PreTreatedData->GetECCTTime(i) << " " << endl;
      }
  }


  for(unsigned int j = 0 ; j < PreTreatedData -> GetGOCCEEMult(); j++) {
    GOCCE_E.push_back(PreTreatedData->GetGOCCEEEnergy(j));
    GOCCE_CloverNumber.push_back(PreTreatedData->GetGOCCEEClover(j));
    GOCCE_CristalNumber.push_back(PreTreatedData->GetGOCCEECristal(j));
    GOCCE_SegmentNumber.push_back(PreTreatedData->GetGOCCEESegment(j));
  }


  //int NumberOfHitClover = 0;

  int DetectorID = -1;

  for( unsigned short i = 0 ; i < PreTreatedData->GetECCEMult() ; i++ )
    {
      // cout << PreTreatedData->GetECCEClover(i) << endl;
      if( PreTreatedData->GetECCEClover(i) != DetectorID)
  {
    if(i==0)
      {
        NumberOfHitClover++;
      }
    else if(PreTreatedData->GetECCEClover(i)!= PreTreatedData->GetECCEClover(i-1) )
      {
        NumberOfHitClover++;
      }
  }
      if(NumberOfHitClover == 4) break;
      //clover_mult -> Fill(NumberOfHitClover);

    }

  //cout << "NumberOfHitClover " << NumberOfHitClover << endl;

  map<int, vector<int> > MapCristal;
  map<int, vector<int> > MapSegment;

  map<int, vector<int> > :: iterator it;    // iterator used with MapCristal
  map<int, vector<int> > :: iterator at;    // iterator used with MapSegment

  vector<int> PositionOfCristal_Buffer_ECC;
  vector<int> PositionOfSegment_Buffer_GOCCE;


  //Fill map Cristal
  for(int clo = 0; clo < NumberOfClover; clo++)
    {
      for(unsigned int k = 0; k < ECC_CloverNumber.size(); k++)
  {
    if(ECC_CloverNumber.at(k) == clo) // && ECC_CristalNumber.at(k)== cri )
      PositionOfCristal_Buffer_ECC.push_back(k);
  }
      if(PositionOfCristal_Buffer_ECC.size() != 0) MapCristal[clo] = PositionOfCristal_Buffer_ECC;

      PositionOfCristal_Buffer_ECC.clear();

    }


  //Fill map Segment
  for(int clo = 0; clo < NumberOfClover; clo++)
    {
      for(int cri = 0; cri < 4 ; cri++)
  {
    //  for(int seg = 0; seg < 4 ; seg++)
      {
        for(unsigned int m = 0; m < GOCCE_CloverNumber.size(); m++)
    {
      if(GOCCE_CloverNumber.at(m) == clo && GOCCE_CristalNumber.at(m) == cri)// && GOCCE_SegmentNumber.at(m) == seg)
        {
          // PositionOfSegment_Buffer_GOCCE.push_back(4*clo+cri);
          PositionOfSegment_Buffer_GOCCE.push_back(m);
        }
    }
      }
      if(PositionOfSegment_Buffer_GOCCE.size() != 0) MapSegment[4*clo+cri] = PositionOfSegment_Buffer_GOCCE;

      PositionOfSegment_Buffer_GOCCE.clear();
  }
    }

  // Treatment
  for(int clo = 0; clo < NumberOfClover ; clo++)
    {
      double E = 0; double T = 0;
      int mult_cristal = 0;
      int cristal = -1 , segment;

      int cristal_Emax = 0; int cristal_Emin = 0;
      int Emax = 0, Emin = 1000000;
      int Tmin = 0, Tmax = 0;

      //ADD-BACK
      it = MapCristal.find(clo);

      int cristal_cond = 0;

      if(it != MapCristal.end())
  {
    vector<int> PositionOfCristal = it -> second;

    mult_cristal = PositionOfCristal.size();
    //if(mult_cristal!=0) cristal_mult -> Fill(mult_cristal);

    // ADD-BACK
    //cout << "boucle" << endl;

    for(unsigned int k = 0; k < PositionOfCristal.size(); k++)
      {
        int indice = PositionOfCristal.at(k);

        cristal_cond += ECC_CristalNumber.at(indice);
        // cout <<  ECC_CristalNumber.at(k) << " " ECC_E.at(k) << endl;

        if(mult_cristal < 3)
    {
      E+= ECC_E.at(indice);

      if(ECC_E.at(indice) < Emin) {
        cristal_Emin = ECC_CristalNumber.at(indice);
        Emin = ECC_E.at(indice);
        Tmin = ECC_T.at(indice);
      }

      if(ECC_E.at(indice) > Emax) {
        cristal_Emax = ECC_CristalNumber.at(indice);
        Emax = ECC_E.at(indice);
        Tmax = ECC_T.at(indice);
      }
    }

        else // case of multiplicity = 3 or 4
    {
      E = -1; cristal_Emax = -1; cristal_Emin = -1; Tmax = -1; Tmin = -1;
    }

        // cout << ECC_E.at(indice) << " " << Emax << " " << cristal_Emax << " " << Emin << " " << cristal_Emin << endl;

      }

    if( (mult_cristal==1) || (mult_cristal ==2  && cristal_cond %2 == 1) )
      {
        // cout << cristal_cond << endl;

        //cristal = cristal_Emax; T = Tmax;
        //cout << Emax << " " << cristal_Emax << " " << Emin << " " << cristal_Emin << endl;

        if(E > 500) { cristal = cristal_Emax; T = Tmax; }
        else        { cristal = cristal_Emin; T = Tmin; }


        // DOPPLER CORRECTION

        at = MapSegment.find(4*clo+cristal);
        segment = -1;

        if(at != MapSegment.end())
    {
      vector<int> PositionOfSegment = at -> second;     // position of segment k in the vector

      int segment_max = -1, E_temp = -1;

      for(unsigned int m = 0; m < PositionOfSegment.size(); m++)             // loop on hit segments of cristal cri of
  clover clo
        {
          int indice = PositionOfSegment.at(m);

          if(GOCCE_E.at(indice) > 0 && GOCCE_CristalNumber.at(indice) == cristal)
      {
        if( GOCCE_E.at(indice) > E_temp )
          {
            segment_max = GOCCE_SegmentNumber.at(indice) ;
            E_temp = GOCCE_E.at(indice);
          }
      }
        }
      segment = segment_max;
    }

      }


    if(E > 0 && cristal != -1 && segment != -1)
      {
        TotalEnergy_lab.push_back(E);
        Time.push_back(T);
        CloverNumber.push_back(clo);
        CristalNumber.push_back(cristal);
        SegmentNumber.push_back(segment);

        double theta = GetSegmentAngleTheta(clo, cristal, segment);

        Theta.push_back(theta);

        double doppler_E = DopplerCorrection(E, theta);
        DopplerCorrectedEnergy.push_back(doppler_E);

        //  cout << E  << " " << clo << " " << cristal << " " << segment << " " << theta << " " << doppler_E << endl;

      }

  }  // end of condition over CristalMap

    } // loop over NumberOfClover

  CloverMult = GetClover_Mult();

  //cout << "Exogam fine" << endl;
  */
}

double TExogamPhysics::ComputeMeanFreePath(double Energy){
  auto b = Map_PhotonCS.lower_bound(Energy);
  auto a = prev(b);
  if(b == Map_PhotonCS.begin()){
    a = b;
    b++;
  }
  else if(b == Map_PhotonCS.end()){
    b--;
    a = prev(b);
  }
  double coeff = (Energy - a->first)/(b->first - a->first);

  double PhotonCrossSection = a->second + coeff*(b->second - a->second); // mm2/g
  return 1./(GeDensity*PhotonCrossSection);
}

// unsigned int TExogamPhysics::GetFlangeNbr(unsigned int crystal_nbr){

// }

double TExogamPhysics::DopplerCorrection(double E, double Theta) {
  double Pi = 3.141592654;
  TString filename = "configs/beta.txt";
  ifstream file;
  // cout << filename << endl;
  file.open(filename);
  if (!file)
    cout << filename << " was not opened" << endl;

  double E_corr = 0;
  double beta = 0.;
  file >> beta;
  double gamma = 1. / sqrt(1 - beta * beta);

  E_corr = gamma * E * (1. - beta * cos(Theta * Pi / 180.));

  return (E_corr);
}

///////////////////////////////////////////////////////////////////////////

void TExogamPhysics::Clear() {
  // Exogam_struc = {};
  EventMultiplicity = 0;
  ECC_Multiplicity = 0;
  GOCCE_Multiplicity = 0;
  NumberOfHitClover = 0;
  NumberOfHitCristal = 0;

  E_Cal.clear();
  E_Doppler.clear();
  Flange_N.clear();
  Crystal_N.clear();
//
//  ECC_CloverNumber.clear();
//  ECC_CristalNumber.clear();
//  GOCCE_CloverNumber.clear();
//  GOCCE_CristalNumber.clear();
//  GOCCE_SegmentNumber.clear();
//
//  // ECC
//  ECC_E.clear();
//  ECC_T.clear();
//
//  // GOCCE
//  GOCCE_E.clear();
//
//  CristalNumber.clear();
//  SegmentNumber.clear();
//  CloverNumber.clear();
//
//  TotalEnergy_lab.clear();
//  Time.clear();
//  DopplerCorrectedEnergy.clear();
//  Position.clear();
//  Theta.clear();
}
///////////////////////////////////////////////////////////////////////////

////	Innherited from VDetector Class	////

//	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
void TExogamPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Exogam clover found " << endl;

  // FIXME ANGLE FILE??? NOT SURE I GET IT...
  // For Doppler I guess... Something like that should be added later
  // But maybe the more stand R,THETA,PHI or X,Y,Z 
  // vector<string> token = {"ANGLE_FILE"};
  vector<string> token = {"Board, Flange, Channel0, Channel1"};
  // FIXME To be implemented in the future
  // vector<string> token = {"Board, Flange, Channel0, Channel1, R, THETA, PHI"};

  //for (unsigned int i = 0; i < blocks.size(); i++) {
  //  if (blocks[i]->HasTokenList(token)) {
  //    int Board, Flange, Channel0, Channel1; // FIXME!!!! Should come from Data...
  //    AddClover(Board, Flange, Channel0, Channel1);
  //  }
  //  else {
  //    cout << "ERROR: check your input file formatting " << endl;
  //    exit(1);
  //  }
  //}
  ReadAnalysisConfig();
}

void TExogamPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;


  // path to photon cross section
  string CSFilename = "../../Inputs/PhotonCrossSection/CoherentGe.xcom";
  string LineBuffer;

  ifstream CSFile;
  CSFile.open(CSFilename.c_str());
  
  if (!CSFile.is_open()) {
    cout << " No CS file found "
         << CSFilename << endl;
    return;
  }
  while(CSFile.good()){
    double gammaE, CrossSection;
    getline(CSFile, LineBuffer);
    istringstream ss(LineBuffer);
    ss >> gammaE >> CrossSection; // E in MeV, converted to keV, CrossSection in cm2/g  
    gammaE *= 1000.; // Convertion to keV
    CrossSection *= 100.;
    Map_PhotonCS[gammaE] = CrossSection;
  }
  // path to file
  string FileName = "./configs/ConfigExogam.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigExogam.dat found: Default parameters loaded for "
            "Analysis "
         << FileName << endl;
    return;
  }
  
  
  string DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 12, "ConfigExogam") == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {

      whatToDo = "";
      AnalysisConfigFile >> whatToDo;
      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }

      else if (whatToDo == "EXO_Threshold") {
        //AnalysisConfigFile >> DataBuffer;
        //m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        //cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
      }
      else if (whatToDo=="MAP_EXO") {
        unsigned int CrystalNb;
        unsigned int Flange;
        unsigned int CrystalNb2;
        
        AnalysisConfigFile >> DataBuffer;
        CrystalNb = stoi(DataBuffer);
        
        AnalysisConfigFile >> DataBuffer;
        Flange = stoi(DataBuffer);
        
        AnalysisConfigFile >> DataBuffer;
        CrystalNb2 = stoi(DataBuffer);
        MapCrystalFlangeCLover[CrystalNb] = std::make_pair(Flange,CrystalNb2);
        // cout << whatToDo << " " << atoi(DataBuffer.substr(0,1).c_str()) << " " << atoi(DataBuffer.substr(1,1).c_str()) << endl;
      }
      else{
        ReadingStatus = false;
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitSpectra() { m_Spectra = new TExogamSpectra(NumberOfClover); }

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ClearSpectra() {
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TExogamPhysics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}


//////////////////////////////////////////////////////////////////////////
void TExogamPhysics::AddClover(int Board, int Flange, int Channel0, int Channel1) {

}

// FIXME Legacy thing... Might delete later
//////////////////////////////////////////////////////////////////////////
// void TExogamPhysics::AddClover(string AngleFile) {
//   ifstream file;
//   //  TString filename = Form("posBaptiste/angles_exogam_clover%d.txt",NumberOfClover);
//   //  TString filename = Form("posz42_simu50mm/angles_exogam_clover%d.txt",NumberOfClover);
//   //  TString filename = Form("posz42_exp_stat_demiring/angles_exogam_clover%d.txt",NumberOfClover);

//   string path = "configs/";
//   TString filename = path + AngleFile;

//   cout << filename << endl;
//   file.open(filename);
//   if (!file)
//     cout << filename << " was not opened" << endl;

//   vector<double> Angles;
//   vector<vector<double>> Segment_angles;
//   vector<vector<vector<double>>> Cristal_angles;

//   Cristal_angles.clear();

//   double angle;
//   string buffer;

//   for (int i = 0; i < 4; i++) {
//     Segment_angles.clear();

//     for (int j = 0; j < 4; j++) {
//       Angles.clear();

//       for (int k = 0; k < 2; k++) {
//         file >> buffer >> angle;

//         Angles.push_back(angle); // Theta (k = 0)   Phi (k = 1)

//         // cout << angle << endl;
//         if (Angles.size() == 2)
//           cout << "Clover " << NumberOfClover << ": Theta=" << Angles[0] << " Phi=" << Angles[1] << endl;
//       }

//       Segment_angles.push_back(Angles);
//     }

//     Cristal_angles.push_back(Segment_angles);
//   }

//   Clover_Angles_Theta_Phi.push_back(Cristal_angles);

//   file.close();

//   NumberOfClover++;
// }

//	Add Parameter to the CalibrationManger
void TExogamPhysics::AddParameterToCalibrationManager() {

  CalibrationManager* Cal = CalibrationManager::getInstance();

  for (auto it = MapCrystalFlangeCLover.begin(); it != MapCrystalFlangeCLover.end(); it++)
   {  unsigned int i = it->first;
      Cal->AddParameter("EXO", "E" + NPL::itoa(i),
                        "EXO_E" + NPL::itoa(i));
      Cal->AddParameter("EXO", "EHG" + NPL::itoa(i),
                        "EXO_EHG" + NPL::itoa(i));
      // Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T",
                        // "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T");

    for (int j = 0; j < 4; j++) {
      Cal->AddParameter("EXO", "Outer" + NPL::itoa(i) + "_" + NPL::itoa(j),
                        "EXO_Outer" + NPL::itoa(i) + "_" + NPL::itoa(j));
    }
  }
}

//	Activated associated Branches and link it to the private member DetectorData address
//	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
void TExogamPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Exogam", true);
  inputChain->SetBranchStatus("fEXO_*", true);
  inputChain->SetBranchAddress("Exogam", &m_EventData);

  /*
  TList* outputList = RootOutput::getInstance()->GetList();
   clover_mult = new TH1F("clover_mult","clover_mult",20,0,20);
    outputList->Add(clover_mult);
  cristal_mult = new TH1F("cristal_mult","cristal_mult",20,0,20);
  outputList->Add(cristal_mult);
  */
  }
}

/////////////////////////////////////////////////////////////////////
//   Activated associated Branches and link it to the private member DetectorPhysics address
//   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
void TExogamPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("EventMultiplicty", true);
  inputChain->SetBranchStatus("ECC_Multiplicity", true);
  inputChain->SetBranchStatus("GOCCE_Multiplicity", true);
  inputChain->SetBranchStatus("ECC_CloverNumber", true);
  inputChain->SetBranchStatus("ECC_CristalNumber", true);
  inputChain->SetBranchStatus("GOCCE_CloverNumber", true);
  inputChain->SetBranchStatus("GOCCE_CristalNumber", true);
  inputChain->SetBranchStatus("GOCCE_SegmentNumber", true);
  inputChain->SetBranchStatus("ECC_E", true);
  inputChain->SetBranchStatus("ECC_T", true);
  inputChain->SetBranchStatus("GOCCE_E", true);
  inputChain->SetBranchStatus("CristalNumber", true);
  inputChain->SetBranchStatus("SegmentNumber", true);
  inputChain->SetBranchStatus("CloverNumber", true);
  inputChain->SetBranchStatus("CloverMult", true);
  inputChain->SetBranchStatus("TotalEnergy_lab", true);
  inputChain->SetBranchStatus("Time", true);
  inputChain->SetBranchStatus("DopplerCorrectedEnergy", true);
  inputChain->SetBranchStatus("Position", true);
  inputChain->SetBranchStatus("Theta", true);
  inputChain->SetBranchAddress("Exogam", &m_EventPhysics);
  }
}

/////////////////////////////////////////////////////////////////////

//	Create associated branches and associated private member DetectorPhysics address
void TExogamPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Exogam", "TExogamPhysics", &m_EventPhysics);

  // control histograms if needed
  /*
  TList* outputList = RootOutput::getInstance()->GetList();
  controle = new TH1F("controle","histo de controle",20,0,20);
  outputList->Add(controle);
  */
}

void TExogamPhysics::SetTreeReader(TTreeReader* TreeReader) {
   TExogamPhysicsReader::r_SetTreeReader(TreeReader);
 }

/////////////////////////////// DoCalibration Part //////////////////////////:

void TExogamPhysics::InitializeRootHistogramsCalib() {
  std::cout << "Initialize Exogam Histograms" << std::endl;
  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;
  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsE_F(it->first);
    }
  }
  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsEHG_F(it->first);
    }
  }
  //for (it = DoCalibrationT.begin(); it != DoCalibrationT.end(); it++) {
  //  if (it->second) {
  //    InitializeRootHistogramsT_F(it->first);
  //  }
  //}
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
      if (it->second) {
        InitializeRootHistogramsOuter_F(it2->first,it->first);
      }
    }
  }
}

void TExogamPhysics::FillHistogramsCalib() {
  if (NPOptionManager::getInstance()->IsReader())
    m_EventData = &(**r_ReaderEventData);
  
  FillRootHistogramsCalib_F();
}

void TExogamPhysics::InitializeRootHistogramsE_F(unsigned int DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOE = Form("EXO_E%d", DetectorNumber);
  TString htitleEXOE = Form("EXO_E%d", DetectorNumber);
  (*TH1Map)["Exogam"][hnameEXOE] = new TH1F(hnameEXOE, htitleEXOE, 65536, 0, 65536);
}

void TExogamPhysics::InitializeRootHistogramsOuter_F(unsigned int DetectorNumber, unsigned int OuterNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOOuter = Form("EXO_Outer%d_%d", DetectorNumber, OuterNumber);
  TString htitleEXOOuter = Form("EXO_Outer%d_%d", DetectorNumber, OuterNumber);
  (*TH1Map)["Exogam"][hnameEXOOuter] = new TH1F(hnameEXOOuter, htitleEXOOuter, 65536, 0, 65536);
}

void TExogamPhysics::InitializeRootHistogramsEHG_F(unsigned int DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOEHG = Form("EXO_EHG%d", DetectorNumber);
  TString htitleEXOEHG = Form("EXO_EHG%d", DetectorNumber);
  (*TH1Map)["Exogam"][hnameEXOEHG] = new TH1F(hnameEXOEHG, htitleEXOEHG, 65536, 0, 65536);

}
void TExogamPhysics::FillRootHistogramsCalib_F(){
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  TString hname;
  
  for (UShort_t i = 0; i < m_EventData->GetExoMult(); i++) {
    unsigned int DetectorNbr = m_EventData->GetExoCrystal(i);

    if(DoCalibrationE[DetectorNbr] && m_EventData->GetExoE(i) >0){
      hname = Form("EXO_E%d", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoE(i));
    }
    if(DoCalibrationEHG[DetectorNbr] && m_EventData->GetExoEHG(i) >0){
      hname = Form("EXO_EHG%d", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoEHG(i));
    }
    if(DoCalibrationOuter[DetectorNbr][0] && m_EventData->GetExoOuter1(i) >0){
      hname = Form("EXO_Outer%d_0", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter1(i));
    }
    if(DoCalibrationOuter[DetectorNbr][1] && m_EventData->GetExoOuter2(i) >0){
      hname = Form("EXO_Outer%d_1", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter2(i));
    }
    if(DoCalibrationOuter[DetectorNbr][2] && m_EventData->GetExoOuter3(i) >0){
      hname = Form("EXO_Outer%d_2", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter3(i));
    }
    if(DoCalibrationOuter[DetectorNbr][3] && m_EventData->GetExoOuter4(i) >0){
      hname = Form("EXO_Outer%d_3", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter4(i));
    }
  }
}

void TExogamPhysics::DoCalibration() {
  std::cout << "Do Calibration Exogam" << std::endl;
  DefineCalibrationSource(Source_name);
  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;
  
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  std::string make_folder = "mkdir " + Path + OutputName;
  
  MakeInitialCalibFolder(make_folder);
  
  ofstream* calib_file = new ofstream;
  ofstream* dispersion_file = new ofstream;
  
  if(!DoCalibrationE.empty()){
    MakeECalibFolders(make_folder);
    CreateCalibrationEFiles(calib_file, dispersion_file);
  }
  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      DoCalibrationE_F(it->first,"E", calib_file, dispersion_file);
    }
  }
  calib_file->close();
  dispersion_file->close();

  if(!DoCalibrationEHG.empty()){
    MakeEHGCalibFolders(make_folder);
    CreateCalibrationEHGFiles(calib_file, dispersion_file);
  }
  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      DoCalibrationE_F(it->first,"EHG", calib_file, dispersion_file);
    }
  }
  calib_file->close();
  dispersion_file->close();
  
  if(!DoCalibrationOuter.empty()){
    MakeOuterCalibFolders(make_folder);
    CreateCalibrationOuterFiles(calib_file, dispersion_file);
  }
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
      if (it->second) {
        DoCalibrationE_F(it->first,Form("Outer%d_",it2->first), calib_file, dispersion_file);
      }
    }
  }
  calib_file->close();
  dispersion_file->close();
}

void TExogamPhysics::MakeInitialCalibFolder(std::string make_folder) {
  int sys = system(make_folder.c_str());
}

void TExogamPhysics::MakeECalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_E").c_str()); 
}

void TExogamPhysics::MakeEHGCalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_EHG").c_str());
}
void TExogamPhysics::MakeOuterCalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_Outer").c_str());
}

void TExogamPhysics::DoCalibrationE_F(unsigned int DetectorNumber,std::string CalibType, ofstream* calib_file, ofstream* dispersion_file) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();

#if CUBIX
  CubixEnergyCal->Reset();
  std::string hnameEXOE = Form("EXO_%s%d",CalibType.c_str(), DetectorNumber);
  std::string htitleEXOE = Form("EXO_%s%d",CalibType.c_str(), DetectorNumber);
  
  auto hist = ((*TH1Map)["Exogam"][hnameEXOE]);

  CubixEnergyCal->SetDataFromHistTH1(hist,0);
    
  for (auto ie : Source_E)
    CubixEnergyCal->AddPeak(ie);
  
  CubixEnergyCal->SetGain(1.);
  CubixEnergyCal->SetVerbosityLevel(1);
    
  CubixEnergyCal->SetFitPlynomialOrder(FitPolOrder);
  CubixEnergyCal->SetNoOffset(false);
  CubixEnergyCal->UseLeftTail(true);
  CubixEnergyCal->UseRightTail(true);

  CubixEnergyCal->UseFirstDerivativeSearch();
  
  CubixEnergyCal->SetGlobalChannelLimits(hist->GetXaxis()->GetBinLowEdge(1),hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetNbins()));      // limit the search to this range in channels
  CubixEnergyCal->SetGlobalPeaksLimits(15,5);   // default fwhm and minmum amplitude for the peaksearch [15 5]

  CubixEnergyCal->StartCalib();
  vector < Fitted > FitResults = CubixEnergyCal->GetFitResults();
  
  
  
  std:: cout << calib_file << " " << (*calib_file).is_open() << std::endl;
  std:: cout << hnameEXOE << " ";
  (*calib_file) << hnameEXOE << " ";
  if(FitResults.size() > 1)
    {
    for(unsigned int i = 0; i <= FitPolOrder; i++){
      (*calib_file) << scientific << setprecision(6) << setw(14) << CubixEnergyCal->fCalibFunction->GetParameter(i) << " ";
      std::cout << scientific << setprecision(6) << setw(14) << CubixEnergyCal->fCalibFunction->GetParameter(i) << " ";
    }
  }
  else
    {
    for(unsigned int i = 0; i <= FitPolOrder; i++){
      (*calib_file) << scientific << setprecision(6) << setw(14) << 0. << " ";
      std::cout << scientific << setprecision(6) << setw(14) << 0. << " ";
    }
  }
  (*calib_file) << "\n";
  std::cout << "\n";

    
  if(FitResults.size()>1 && CubixEnergyCal->fCalibFunction) {
      auto c = new TCanvas;
      c->SetName("CalibrationResults");
      c->SetTitle("Calibration Results");
      c->Divide(1,2,0.0001,0.0001);
      c->cd(1);
      CubixEnergyCal->fCalibGraph->Draw("ap");
      CubixEnergyCal->fCalibFunction->Draw("same");
      c->cd(2);
      CubixEnergyCal->fResidueGraph->Draw("ape");
      c->Update();
      c->Modified();
  }
  if(FitResults.size() > 1)
  {
  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fCalibGraph->Clone());
  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Energy (MeV)");
  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Calibration_Graph_%s",hnameEXOE.c_str()));

  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fResidueGraph->Clone());
  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetXaxis()->SetTitle("Energy (MeV)");
  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Residue (MeV)");
  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Residue_Graph_%s",hnameEXOE.c_str()));
  }
#else
  std::cout << "Exogam calibration currently not supported without CUBIX. Download CUBIX and set -DCUBIX=1 to use EXOGAM calibration\n";
  exit(1);

#endif

}

void TExogamPhysics::DefineCalibrationSource(std::string source) {
  // 239Pu
  if(source == "60Co"){
    Source_isotope.push_back("$^{60}$Co");
    Source_E.push_back(1.17322);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(99.85);
    
    Source_isotope.push_back("$^{60}$Co");
    Source_E.push_back(1.33249);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(99.98);
  }
  else if(source == "152Eu"){
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.121782);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(28.58);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.344279);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(26.5);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.40801);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(21.0);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.964079);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(14.6);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.11207);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(13.64);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.778904);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(12.94);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.08587);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(10.21);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.244698);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(7.58);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.867378);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(4.25);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.443965);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(2.82);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.411116);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(2.23);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.08974);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(1.73);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.29914);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(1.62);
    
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.21295);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(1.42);
    
  }
  else{
    std::cout << "Please enter a valid source for gamma ray calibration\nCurrently supported sources are 60Co and 152Eu\n";
    exit(1);
  }
}


// FIXME Probably could be done better, currently a but inelegant
void TExogamPhysics::CreateCalibrationEFiles(ofstream* calib_file,
                                                 ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_E";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_E/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_E/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::CreateCalibrationEHGFiles(ofstream* calib_file,
                                                 ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_EHG";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_EHG/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_EHG/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::CreateCalibrationOuterFiles(ofstream* calib_file,
                                                 ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_Outer";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_Outer/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_Outer/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::ReadDoCalibration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");

  vector<string> calibs = {"FirstCr","LastCr","FirstOuter","LastOuter","FitOrder","Source"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(calibs)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Exogam Calibration" << endl;
      unsigned int FirstCr = blocks[i]->GetInt("FirstCr");
      unsigned int LastCr = blocks[i]->GetInt("LastCr");
      unsigned int FirstOuter = blocks[i]->GetInt("FirstOuter");
      unsigned int LastOuter = blocks[i]->GetInt("LastOuter");
      FitPolOrder = blocks[i]->GetInt("FitOrder");
      Source_name = blocks[i]->GetString("Source");
      for(unsigned int k = FirstCr; k <= LastCr; k++){
        DoCalibrationE[k] = true;
        DoCalibrationEHG[k] = true; 
        for(unsigned int p = FirstOuter; p <= LastOuter; p++){
          DoCalibrationOuter[k][p] = true; 
      }
    }
    }
    else {
      cout << "ERROR: Missing token for Exogam DoCalibration blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
}

void TExogamPhysics::WriteHistogramsCalib() {
  std::cout << "Writing Exogam Histograms\n";
  WriteHistogramsE();
  RootHistogramsCalib::getInstance()->GetFile()->Close();
}

void TExogamPhysics::WriteHistogramsE() {
  auto File = RootHistogramsCalib::getInstance()->GetFile();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();

  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;
  std::string hnameEXOE;

  if (!File->GetDirectory("Exogam"))
    File->mkdir("Exogam");
  File->cd("Exogam");

  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      hnameEXOE = Form("EXO_E%d", it->first);
  
      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it->first));
      gDirectory->cd(Form("EXO_Cr%d", it->first));
      
      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      if((*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]->Write();
      if((*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]->Write();
    }
    File->cd("Exogam");
  }
  
  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      hnameEXOE = Form("EXO_EHG%d", it->first);
  
      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it->first));
      gDirectory->cd(Form("EXO_Cr%d", it->first));
      
      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      if((*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_EHG%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_EHG%d",it->first)]->Write();
      if((*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_EHG%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_EHG%d",it->first)]->Write();
    }
    File->cd("Exogam");
  }
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
    if (it->second) {
      hnameEXOE = Form("EXO_Outer%d_%d",it2->first,it->first);
  
      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it2->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it2->first));
      gDirectory->cd(Form("EXO_Cr%d", it2->first));
      
      if (!gDirectory->GetDirectory("Outers"))
        gDirectory->mkdir("Outer");
      gDirectory->cd("Outer");
      
      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      if((*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->Write();
      if((*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->Write();
    }
    File->cd("Exogam");
  }
  }
}
///////////////////////////////////////////////////////////////////////////
namespace EXOGAM_LOCAL {
  //	tranform an integer to a string
  double fEXO_E(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXO/E";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoE(i));
  }
  
  double fEXO_EHG(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXO/EHG";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoEHG(i));
  }
  
  double fEXO_T(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    name += "_TDC";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoTDC(i));
  }
  
  double fEXO_Outer(const TExogamData* m_EventData, const unsigned int& i, const unsigned int OuterNumber) {
    static string name;
    name = "EXO/Outer";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    name += "_";
    name += NPL::itoa(OuterNumber);
    if(OuterNumber == 0)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter1(i));
    else if(OuterNumber == 1)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter2(i));
    else if(OuterNumber == 2)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter3(i));
    else if(OuterNumber == 3)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter4(i));
    else{
      std::cout << "WARNING: Outer number != 0-3, something is wrong\n";
      return 0;
    };
  }
  
  string itoa(int value) {
    std::ostringstream o;

    if (!(o << value))
      return "";

    return o.str();
  }
} // namespace EXOGAM_LOCAL

/////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TExogamPhysics::Construct() { return (NPL::VDetector*)new TExogamPhysics(); }

NPL::VTreeReader* TExogamPhysics::ConstructReader() { return (NPL::VTreeReader*)new TExogamPhysicsReader(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_exogam {
 public:
  proxy_exogam() {
    NPL::DetectorFactory::getInstance()->AddToken("Exogam", "Exogam");
    NPL::DetectorFactory::getInstance()->AddDetector("Exogam", TExogamPhysics::Construct);
    NPL::DetectorFactory::getInstance()->AddDetectorReader("Exogam", TExogamPhysics::ConstructReader);
  }
};

proxy_exogam p;
}
