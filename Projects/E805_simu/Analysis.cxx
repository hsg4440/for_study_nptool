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
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
  
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
  
  ProtonSi = NPL::EnergyLoss(Path+ "proton_Si.G4table", "G4Table", 100);
  
  for(unsigned int i = 0; i < ParticleType.size(); i++){
    LightAl[ParticleType[i]] = NPL::EnergyLoss(Path+ParticleType[i]+"_Al.G4table","G4Table",100);
    LightTarget[ParticleType[i]] = NPL::EnergyLoss(Path+ParticleType[i]+"_CH2.G4table","G4Table",100);
  }
  BeamTarget["48Cr"] = NPL::EnergyLoss(Path+"Cr48_CH2.G4table","G4Table",100);

  Reaction_pd->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_pd->GetBeamEnergy(),TargetThickness*0.5,0));
  Reaction_pt->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_pt->GetBeamEnergy(),TargetThickness*0.5,0));
  Reaction_p3He->SetBeamEnergy(BeamTarget["48Cr"].Slow(Reaction_p3He->GetBeamEnergy(),TargetThickness*0.5,0));
  reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction->GetBeamEnergy();
  

  Cal = CalibrationManager::getInstance();  


//
}
  ///////////////////////////// Initialize some important parameters //////////////////////////////////


bool Analysis::UnallocateBeforeBuild(){
    return true;
}

bool Analysis::UnallocateBeforeTreat(){
  return true;
}

bool Analysis::FillOutputCondition(){
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){


    ReInit();
    //////////////////// MUST2 Part ////////////////////
    int M2_size = M2->Si_E.size();
    for(unsigned int countMust2 = 0 ; countMust2 < M2_size ; countMust2++){
      M2_TelescopeM++;
      // MUST2
      int TelescopeNumber = M2->TelescopeNumber[countMust2];

      // Part 1 : Impact Angle
      ThetaM2Surface = 0;
      ThetaNormalTarget = 0;
      

      BeamImpact = TVector3(0,0,0); 
      
      BeamDirection = TVector3(0,0,1);
      
      TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
      M2_ThetaLab.push_back(HitDirection.Angle( BeamDirection ));

      M2_X.push_back(M2 -> GetPositionOfInteraction(countMust2).X());
      M2_Y.push_back(M2 -> GetPositionOfInteraction(countMust2).Y());
      M2_Z.push_back(M2 -> GetPositionOfInteraction(countMust2).Z());

      ThetaM2Surface = HitDirection.Angle(- M2 -> GetTelescopeNormal(countMust2) );
      ThetaNormalTarget = HitDirection.Angle( TVector3(0,0,1) ) ;

      // Part 2 : Impact Energy
      int CristalNb = M2->CsI_N[countMust2];

      Si_E_M2 = M2->Si_E[countMust2];
      CsI_E_M2= M2->CsI_E[countMust2];
      // if(CsI_E_M2 > 0)
      // std::cout << "Analysis " << CsI_E_M2 << " " << M2->CsI_E.size() << " " << CristalNb << "\n \n";
      Energy = 0;
      if(CsI_E_M2 > 0 && Si_E_M2 > 0)
        {
        Energy += CsI_E_M2;
        Energy += Si_E_M2;
        }
      else
        Energy = -1000;
      if(Energy > 0)
        Energy = LightTarget["deuteron"].EvaluateInitialEnergy(Energy,TargetThickness*0.5, ThetaNormalTarget);


      // Evaluate energy using the thickness
      //Energy = LightAl.EvaluateInitialEnergy( Energy ,0.4*micrometer , ThetaM2Surface);
      // Target Correction
      //M2_ELab.push_back(LightTarget.EvaluateInitialEnergy(Energy ,TargetThickness*0.5, ThetaNormalTarget));
      M2_ELab.push_back(Energy);
      M2_dE.push_back(Si_E_M2);

      // Part 3 : Excitation Energy Calculation
      M2_Ex.push_back(reaction->ReconstructRelativistic( Energy, M2_ThetaLab[countMust2] ));
       
      M2_ThetaLab[countMust2]=M2_ThetaLab[countMust2]/deg;

      // Part 4 : Theta CM Calculation
      M2_ThetaCM.push_back(reaction->EnergyLabToThetaCM( M2_ELab[countMust2] , M2_ThetaLab[countMust2])/deg);

      //if(proton_cut[TelescopeNumber]){
      //  std::cout << "Test :" << proton_cut[TelescopeNumber] << " \n";
      //  
    }//end loop MUST2

  
  //for(unsigned int countMust2 = 0 ; countMust2 < M2->Si_E.size() ; countMust2++){
  //Si_E_M2 = M2->Si_E[countMust2];
  //CsI_E_M2 = M2->CsI_E[countMust2];
  //ThetaM2Surface = 0;
  //  if(Si_E_M2 > 0 && CsI_E_M2 > 8192){
  //    //double EfromdeltaE = ProtonSi.EvaluateEnergyFromDeltaE(
  //    //  Si_E_M2, 300*um, ThetaM2Surface, 6.0 * MeV, 300.0 * MeV,
  //    //  0.001 * MeV, 10000);
  //    double EfromdeltaE = (CsI_E_M2-8192)*0.1;
  //    M2_ECsI_from_deltaE.push_back(EfromdeltaE);
  //    if(EfromdeltaE > 0){
  //      Beta_light = sqrt(1./(1.+1./(pow(EfromdeltaE/911. + 1,2)-1)));
  //      Beta_from_deltaE.push_back(Beta_light);
  //      if(Beta_light>0){
  //        double Beth = log(2*511.*Beta_light*Beta_light/(0.174*(1-Beta_light*Beta_light))) - Beta_light*Beta_light; 
  //        Beth_from_deltaE.push_back(Beth);
  //      }
  //    }
  //  }
  //}
}



void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("M2_TelescopeM",&M2_TelescopeM,"M2_TelescopeM/s");
  // RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_p",&M2_CsI_E_p);
  // RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_d",&M2_CsI_E_d);
  // RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_t",&M2_CsI_E_t);
  // RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_a",&M2_CsI_E_a);
  // RootOutput::getInstance()->GetTree()->Branch("M2_Ex_p",&M2_Ex_p);
   RootOutput::getInstance()->GetTree()->Branch("M2_Ex",&M2_Ex);
  // RootOutput::getInstance()->GetTree()->Branch("M2_Ex_d",&M2_Ex_d);
  // RootOutput::getInstance()->GetTree()->Branch("M2_Ex_t",&M2_Ex_t);
  // RootOutput::getInstance()->GetTree()->Branch("M2_Ex_a",&M2_Ex_a);
  //RootOutput::getInstance()->GetTree()->Branch("ExNoBeam",&ExNoBeam,"ExNoBeam/D");
  //RootOutput::getInstance()->GetTree()->Branch("ExNoProton",&ExNoProton,"ExNoProton/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_ELab",&M2_ELab);
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaLab",&M2_ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaCM",&M2_ThetaCM);
  RootOutput::getInstance()->GetTree()->Branch("M2_X",&M2_X);
  RootOutput::getInstance()->GetTree()->Branch("M2_Y",&M2_Y);
  RootOutput::getInstance()->GetTree()->Branch("M2_Z",&M2_Z);
  RootOutput::getInstance()->GetTree()->Branch("M2_dE",&M2_dE);
  RootOutput::getInstance()->GetTree()->Branch("CsI_E_M2",&CsI_E_M2);
  // RootOutput::getInstance()->GetTree()->Branch("M2_ECsI_from_deltaE",&M2_ECsI_from_deltaE);
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){
  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
}
////////////////////////////////////////////////////////////////////////////////

void Analysis::ReInit(){

  CsI_Energy.clear();
  Si_E_M2 = -1000;
  CsI_E_M2 = -1000;
  Theta_seg = -1000;
  Phi_seg = -1000;
  M2_TelescopeM= 0; 
  M2_Ex.clear();
  
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

