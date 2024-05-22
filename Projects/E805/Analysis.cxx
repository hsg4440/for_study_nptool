/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  MUST_AND_ZDD analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include<iostream>
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  InitInputBranch();
  InitOutputBranch();
  CATS = (TCATSPhysics*)  m_DetectorManager -> GetDetector("CATSDetector");
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
  ZDD = (TZDDPhysics*)  m_DetectorManager -> GetDetector("ZDD");
  TAC = (TTACPhysics*)  m_DetectorManager -> GetDetector("TAC");
  Exogam = (TExogamPhysics*)  m_DetectorManager -> GetDetector("Exogam");
  reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction->GetBeamEnergy();


  string Path = "../../Inputs/EnergyLoss/";
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string beam=  NPL::ChangeNameToG4Standard(reaction->GetNucleus1()->GetName());
  string heavy_ejectile=  NPL::ChangeNameToG4Standard(reaction->GetNucleus4()->GetName());
  string light=NPL::ChangeNameToG4Standard(reaction->GetNucleus3()->GetName());

  string Reaction_pd_s = "48Cr(p,d)47Cr@1620";
  string Reaction_pt_s = "48Cr(p,t)46Cr@1620";
  string Reaction_p3He_s = "48Cr(p,3He)46V@1620";
  Reaction_pd = new Reaction(Reaction_pd_s);
  Reaction_pt = new Reaction(Reaction_pt_s);
  Reaction_p3He = new Reaction(Reaction_p3He_s);

//
  ProtonSi = NPL::EnergyLoss(Path+ "proton_Si.G4table", "G4Table", 100);
  
  for(unsigned int i = 0; i < ParticleType.size(); i++){
   LightAl[ParticleType[i]] = NPL::EnergyLoss(Path+ParticleType[i]+"_Al.G4table","G4Table",100);
   LightTarget[ParticleType[i]] = NPL::EnergyLoss(Path+ParticleType[i]+"_CH2.G4table","G4Table",100);
  }
  BeamTarget["48Cr"] = NPL::EnergyLoss(Path+"Cr48_CH2.G4table","G4Table",100);

  Reaction_pd->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_pd->GetBeamEnergy(),TargetThickness*0.5,0));
  Reaction_pt->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_pt->GetBeamEnergy(),TargetThickness*0.5,0));
  Reaction_p3He->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_p3He->GetBeamEnergy(),TargetThickness*0.5,0));
  Cal = CalibrationManager::getInstance();
  IsPhysics = NPOptionManager::getInstance()->GetInputPhysicalTreeOption(); 
  TCutGMap = RootHistogramsCalib::getInstance()->GetTCutGMap();
  TFileMap = RootHistogramsCalib::getInstance()->GetTFileMap();

  unsigned int NbCsI = 16;
  unsigned int NbDetectors = 4;
  for(unsigned int k = 0; k < NbDetectors; k++){
    for (unsigned int j = 0; j < NbCsI; j++) {
      string CutsPath = NPOptionManager::getInstance()->GetCutsPath();

      for (unsigned int i = 0; i < ParticleTypeCUT.size(); i++) {
        std::cout << ParticleTypeCUT[i] << "\n";
        TString CutName = Form("%s_hMM%u_CSI%u", ParticleTypeCUT[i].c_str(), k+1, j + 1);
        TString cFileName = CutName + ".root";
        std::cout << CutName << "  " << cFileName << " " << CutsPath + cFileName << "\n";


        if((*TFileMap)["MUST2"][CutName] = new TFile(CutsPath+cFileName))
          (*TCutGMap)["MUST2"][CutName] = (TCutG*)(*TFileMap)["MUST2"][CutName]->FindObjectAny(CutName);
      }
    }
  }
}
  ///////////////////////////// Initialize some important parameters //////////////////////////////////


bool Analysis::UnallocateBeforeBuild(){
  GATCONFMASTER.clear();
  GATCONFMASTERTS.clear();

  GATCONFMASTER = **GATCONFMASTER_;
  return (GATCONFMASTER.size() == 1 && GATCONFMASTER[0] > 0); 
  //return true;
}

bool Analysis::UnallocateBeforeTreat(){
  GATCONFMASTERTS = **GATCONFMASTERTS_;
  return true;
}

bool Analysis::FillOutputCondition(){
  // return bCATS;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

    ReInit();
    
    TreatCATS();
    if(bCATS){
      TreatMUST2();
      TreatEXO();
    }
    
    
    
}

void Analysis::TreatCATS(){
  if(CATS->PositionOnTargetX > -1000 && CATS->PositionOnTargetY > -1000){
    BeamImpact = TVector3(CATS->PositionOnTargetX,CATS->PositionOnTargetY,0); 
    BeamDirection = TVector3(CATS->PositionX[0] - CATS->PositionX[1],CATS->PositionY[0] - CATS->PositionY[1],-(CATS->PositionZ[0] - CATS->PositionZ[1]));
    bCATS = true;
  }
  else bCATS = false;
  
}

void Analysis::TreatZDD(){
  
}

void Analysis::TreatTAC(){
  
}

void Analysis::TreatMUST2(){

  int M2_size = M2->Si_E.size();
  for(unsigned int countMust2 = 0 ; countMust2 < M2_size ; countMust2++){
    M2_TelescopeM++;
    int TelescopeNumber = M2->TelescopeNumber[countMust2];

    // Part 1 : Impact Angle
    ThetaM2Surface = 0;
    ThetaNormalTarget = 0;
      
    //BeamImpact = TVector3(0,0,0);
    //BeamDirection = TVector3(0,0,1);
    TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact;
    M2_ThetaLab.push_back(HitDirection.Angle( BeamDirection ));
    //std::cout << BeamImpact.X() << " " << BeamImpact.Y() << " "  << BeamImpact.Z() << std::endl;
    //std::cout << BeamDirection.X() << " " << BeamDirection.Y() << " "  << BeamDirection.Z() << std::endl << std::endl;;

    M2_X.push_back(M2 -> GetPositionOfInteraction(countMust2).X());
    M2_Y.push_back(M2 -> GetPositionOfInteraction(countMust2).Y());
    M2_Z.push_back(M2 -> GetPositionOfInteraction(countMust2).Z());

    ThetaM2Surface = HitDirection.Angle(- M2 -> GetTelescopeNormal(countMust2) );
    ThetaNormalTarget = HitDirection.Angle( TVector3(0,0,1) ) ;

    // Part 2 : Impact Energy
    int CristalNb = M2->CsI_N[countMust2];

    Si_E_M2 = M2->Si_E[countMust2];
    CsI_E_M2= M2->CsI_E[countMust2];
    double Si_X = M2->Si_X[countMust2];
    unsigned int Pixel = M2->Pixel[countMust2];
    unsigned int PID= 0;
    
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(TelescopeNumber);
    name += "_Si_X_O";
    name += NPL::itoa(Si_X);
    name += "_E";
    double StripXEnergy_O = Cal->ApplyCalibration(name,M2->Si_E_Raw[countMust2]);
    // std::cout << Si_E_M2 << " " << StripXEnergy_O << std::endl;

    for(unsigned int i = 0; i < ParticleType.size(); i++){
      Energy[ParticleType[i]] = 0;
      CsI_Energy[ParticleType[i]] = 0;

    if(M2->CsI_E_Raw[countMust2] > 8192){
      // The energy in CsI is calculate form dE/dx Table because
      std::string name;
      if(M2->GetCalPixel() && std::find(ParticleTypePixel.begin(),ParticleTypePixel.end(),ParticleType[i])!= ParticleTypePixel.end()){
        name = "MUST2/"+ParticleType[i]+"_T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_Pixel"+NPL::itoa(Pixel)+"_E";
        // Checking if Pixel calib worked, else using global calib
        // std::cout << name << " " << Cal->GetValue(name,0) << std::endl;
        if(Cal->GetValue(name,0) == 0)
          name = "MUST2/"+ParticleType[i]+"_T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_E";
      }
      else{
        name = "MUST2/"+ParticleType[i]+"_T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_E";

      }
      CsI_Energy[ParticleType[i]] =  Cal->ApplyCalibration(name,M2->CsI_E_Raw[countMust2]);
      Energy[ParticleType[i]] = CsI_Energy[ParticleType[i]];
      
      }
    for(unsigned int i = 0; i < ParticleTypeCUT.size(); i++){
      TString CutName = Form("%s_hMM%u_CSI%u", ParticleTypeCUT[i].c_str(), TelescopeNumber, CristalNb);
      if ((*TCutGMap)["MUST2"][CutName] != 0 &&
          (*TCutGMap)["MUST2"][CutName]->IsInside(M2->CsI_E_Raw[countMust2], StripXEnergy_O)) {
            PID = i+1;
          }
    }

    Energy[ParticleType[i]] += Si_E_M2;
    
    if(Energy[ParticleType[i]] > 0)
      Energy[ParticleType[i]] = LightTarget[ParticleType[i]].EvaluateInitialEnergy(Energy[ParticleType[i]] ,TargetThickness*0.5, ThetaNormalTarget);
    else
      Energy[ParticleType[i]] = -1000;


    }
    // Target Correction
    //M2_ELab.push_back(LightTarget.EvaluateInitialEnergy(Energy ,TargetThickness*0.5, ThetaNormalTarget));
    M2_ELab.push_back(Energy["proton"]);
    M2_dE.push_back(Si_E_M2);

    // Part 3 : Excitation Energy Calculation_npreader_
    M2_Ex_p.push_back(reaction->ReconstructRelativistic( Energy["proton"] , M2_ThetaLab[countMust2] ));
    M2_Ex_d.push_back(Reaction_pd->ReconstructRelativistic( Energy["deuteron"] , M2_ThetaLab[countMust2] ));
    M2_Ex_t.push_back(Reaction_pt->ReconstructRelativistic( Energy["triton"] , M2_ThetaLab[countMust2] ));
    M2_Ex_a.push_back(Reaction_p3He->ReconstructRelativistic( Energy["alpha"] , M2_ThetaLab[countMust2] ));
    
    TLorentzVector PHeavy_pd = Reaction_pd->LorentzAfterReaction(Energy["deuteron"] , M2_ThetaLab[countMust2]);
    TLorentzVector PHeavy_pt = Reaction_pt->LorentzAfterReaction(Energy["triton"] , M2_ThetaLab[countMust2]);
    TLorentzVector PHeavy_p3He = Reaction_p3He->LorentzAfterReaction(Energy["alpha"] , M2_ThetaLab[countMust2]);
    Beta_pd.push_back(PHeavy_pd.Beta());
    Beta_pt.push_back(PHeavy_pt.Beta());
    Beta_p3He.push_back(PHeavy_p3He.Beta());


    M2_CsI_E_p.push_back(CsI_Energy["proton"]);
    M2_CsI_E_d.push_back(CsI_Energy["deuteron"]);
    M2_CsI_E_t.push_back(CsI_Energy["triton"]);
    M2_CsI_E_a.push_back(CsI_Energy["alpha"]);
    Si_E_O.push_back(StripXEnergy_O);
    PID_M2.push_back(PID);
    
    M2_ThetaLab[countMust2]=M2_ThetaLab[countMust2]/deg;

    // Part 4 : Theta CM Calculation
    M2_ThetaCM.push_back(reaction->EnergyLabToThetaCM( M2_ELab[countMust2] , M2_ThetaLab[countMust2])/deg);

  }

}

void Analysis::TreatEXO(){
  int EXO_AB_size = Exogam->E_AB.size();
  for(unsigned int countExo = 0 ; countExo < EXO_AB_size; countExo++){
  // Doing Doppler correction only if one reaction occurs
    if(Beta_pd.size() == 1){
      EXO_Doppler_pd.push_back(Doppler_Correction(Exogam->Theta_D[countExo], Exogam->Phi_D[countExo], 0,0,Beta_pd[0],Exogam->E_AB[countExo]));
      EXO_Doppler_pt.push_back(Doppler_Correction(Exogam->Theta_D[countExo], Exogam->Phi_D[countExo], 0,0,Beta_pt[0],Exogam->E_AB[countExo]));
      EXO_Doppler_p3He.push_back(Doppler_Correction(Exogam->Theta_D[countExo], Exogam->Phi_D[countExo], 0,0,Beta_p3He[0],Exogam->E_AB[countExo]));
    }
  }
}


void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("GATCONF",&GATCONFMASTER);
  RootOutput::getInstance()->GetTree()->Branch("GATCONFTS",&GATCONFMASTERTS);
  
  // Brnch with old alpha calibration: only to make cuts before calib
  RootOutput::getInstance()->GetTree()->Branch("Si_E_O",&Si_E_O);
  RootOutput::getInstance()->GetTree()->Branch("PID_M2",&PID_M2);

  
  RootOutput::getInstance()->GetTree()->Branch("M2_TelescopeM",&M2_TelescopeM,"M2_TelescopeM/s");
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_p",&M2_CsI_E_p);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_d",&M2_CsI_E_d);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_t",&M2_CsI_E_t);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_a",&M2_CsI_E_a);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_p",&M2_Ex_p);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_d",&M2_Ex_d);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_t",&M2_Ex_t);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_a",&M2_Ex_a);
  
  RootOutput::getInstance()->GetTree()->Branch("M2_ELab",&M2_ELab);
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaLab",&M2_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaCM",&M2_ThetaCM);
  RootOutput::getInstance()->GetTree()->Branch("M2_X",&M2_X);
  RootOutput::getInstance()->GetTree()->Branch("M2_Y",&M2_Y);
  RootOutput::getInstance()->GetTree()->Branch("M2_Z",&M2_Z);
  
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_pd",&EXO_Doppler_pd);
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_pt",&EXO_Doppler_pt);
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_p3He",&EXO_Doppler_p3He);
  
  RootOutput::getInstance()->GetTree()->Branch("Beta_pd",&Beta_pd);
  RootOutput::getInstance()->GetTree()->Branch("Beta_p3He",&Beta_p3He);
  RootOutput::getInstance()->GetTree()->Branch("Beta_pt",&Beta_pt);
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){

  if(!NPOptionManager::getInstance()->GetInputPhysicalTreeOption()){
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    GATCONFMASTER_ = new TTreeReaderValue<vector<unsigned int>>(*inputTreeReader,"GATCONF");
    GATCONFMASTERTS_ = new TTreeReaderValue<vector<unsigned long long>>(*inputTreeReader,"GATCONFTS");

  }
}
////////////////////////////////////////////////////////////////////////////////

void Analysis::ReInit(){

  Energy.clear();
  CsI_Energy.clear();
  Si_E_M2 = -1000;
  CsI_E_M2 = -1000;
  Theta_seg = -1000;
  Phi_seg = -1000;
  M2_TelescopeM= 0; 
  M2_Ex_p.clear();
  M2_Ex_d.clear();
  M2_Ex_t.clear();
  M2_Ex_a.clear();

  PID_M2.clear(); 
  Si_E_O.clear(); 
  M2_CsI_E_p.clear();
  M2_CsI_E_d.clear();
  M2_CsI_E_t.clear();
  M2_CsI_E_a.clear();
  
  // M2_ECsI_from_deltaE.clear();
  //ExNoBeam=ExNoProto.clear();
  //EDC.clear();
  M2_ELab.clear();
  //BeamEnergy .clear();
  M2_ThetaLab.clear();
  M2_ThetaCM.clear();
  M2_X.clear();
  M2_Y.clear();
  M2_Z.clear();
  M2_dE.clear();

  EXO_Doppler_p3He.clear();
  EXO_Doppler_pt.clear();
  EXO_Doppler_pd.clear();
  
  Beta_p3He.clear();
  Beta_pt.clear();
  Beta_pd.clear();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy{
  public:
    proxy(){
      NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy p;
}

