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
  InitialConditions = new TInitialConditions;
  ReactionConditions = new TReactionConditions; 
  InitOutputBranch();
  InitInputBranch();
  ZDD= (TZDDPhysics*) m_DetectorManager->GetDetector("ZDD");
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
  reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction->GetBeamEnergy();


  string Path = "../../Inputs/EnergyLoss/";
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string beam=  NPL::ChangeNameToG4Standard(reaction->GetNucleus1()->GetName());
  string heavy_ejectile=  NPL::ChangeNameToG4Standard(reaction->GetNucleus4()->GetName());
  string light=NPL::ChangeNameToG4Standard(reaction->GetNucleus3()->GetName());

  
  Beam_Target = NPL::EnergyLoss(beam+"_"+TargetMaterial+".G4table","G4Table",100);
  Heavy_Target = NPL::EnergyLoss(heavy_ejectile+"_"+TargetMaterial+".G4table","G4Table",100);
  LightTarget = NPL::EnergyLoss(light+"_"+TargetMaterial+".G4table","G4Table",100 );
  LightAl = NPL::EnergyLoss(light+"_Al.G4table","G4Table",100);
  LightSi = NPL::EnergyLoss(light+"_Si.G4table","G4Table",100);
  
  ///////////////////////// Creating EnerlyLoss Objects for ZDD ////////////////////////////////////
  for(int i = 0; i < ZDD->Get_AC_Material().size(); i++){
    Heavy_IC_Mylar.push_back(NPL::EnergyLoss(heavy_ejectile+"_"+ZDD->Get_AC_Material()[i]+".G4table","G4Table",100));
  }
  for(int i = 0; i < ZDD->Get_Entry_Exit_Material().size(); i++){
    Heavy_IC_Windows.push_back(NPL::EnergyLoss(heavy_ejectile+"_"+ZDD->Get_Entry_Exit_Material()[i]+".G4table","G4Table",100));
  }
  for(int i = 0; i < ZDD->Get_m_Gas_Gap_Gas().size(); i++){
    double pressure = ZDD->Get_m_Gas_Gap_Pressure()[i]/bar;
    string pres = to_string(pressure);
    pres.erase(pres.find_last_not_of('0')+1,pres.size());
    if(pressure >= 1){
      pres.erase(pres.find_last_not_of('.')+1,pres.size());
    }
    int Temp = ZDD->Get_m_Gas_Gap_Temperature()[i];
    Heavy_IC_Gas.push_back(NPL::EnergyLoss(heavy_ejectile+"_"+ZDD->Get_m_Gas_Gap_Gas()[i]+"_"+pres+"bar_"+to_string(Temp)+"K.G4table","G4Table",100));

  }
  for(int i = 0; i < ZDD->Get_m_Drift_Chamber_Gas().size(); i++){
    double pressure = ZDD->Get_m_Drift_Chamber_Pressure()[i]/bar;
    string pres = to_string(pressure);
    pres.erase(pres.find_last_not_of('0')+1,pres.size());
    if(pressure >= 1){
      pres.erase(pres.find_last_not_of('.')+1,pres.size());
    }
    int Temp = ZDD->Get_m_Drift_Chamber_Temperature()[i];
    Heavy_DC_Gas.push_back(NPL::EnergyLoss(heavy_ejectile+"_"+ZDD->Get_m_Drift_Chamber_Gas()[i]+"_"+pres+"bar_"+to_string(Temp)+"K.G4table","G4Table",100));

  }
  ///////////////////////// Creating EnergyLoss Objects for MUST2 ////////////////////////////////////

  //////////////////////////// Treating Drift Chamber (for this part, we assume the drift speed, in further analysis we shall use CATS)
  //////////////////////////// (information to retrieve drift coefficients)


  ///////////////////////////// Initialize some important parameters //////////////////////////////////

  Drift_Speed =  1 * cm / microsecond;
  ZDir.SetXYZ(0.,0.,1.);


  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  M2_X = 0;
  M2_Y = 0;
  M2_Z = 0;
  M2_dE = 0;
  BeamDirection = TVector3(0,0,1);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
  if(CheckGoodEvent()){
    ReInit();
    SetParticles();
    
    TVector3 BeamDir = InitialConditions->GetBeamDirection();
    IncidentTheta = BeamDir.Angle(ZDir);
    

    std::vector<double> dE_IC(ZDD->GetICcounter(),0);

    ///////////////////////////// Using Initial Conditions as CATS ///////////////
    // This part seems to be working fine
    double xinit = InitialConditions->GetIncidentPositionX();
    double yinit = InitialConditions->GetIncidentPositionY();
    double zinit = InitialConditions->GetIncidentPositionZ();
    double ztarget = m_DetectorManager->GetTargetZ();
    double coefficient = (ztarget - zinit)/BeamDir.Z(); // coeff such that k*zdir = (ztarget - zinit)
    xtarget = xinit + coefficient*BeamDir.X();
    ytarget = yinit + coefficient*BeamDir.y();
    BeamImpact = TVector3(xtarget,ytarget,ztarget); 
    BeamDirection = TVector3(xtarget - xinit, ytarget - yinit, ztarget - zinit);
    
    FinalBeamEnergy = Beam_Target.Slow(OriginalBeamEnergy,TargetThickness*0.5,BeamDirection.Angle(ZDir)); 
    reaction->SetBeamEnergy(FinalBeamEnergy);


    //////////////////////////// Treating DC Position ///////////////////////////
    
    for(int i = 0; i < ZDD->DC_DetectorNumber.size(); i++){
      if(ZDD->DC_DetectorNumber[i] == 1){
        ZDD_DC_X = Drift_Speed*(ZDD->DC_DriftTime[i]*microsecond) - Drift_Chamber_Length/2;
      }
      else if(ZDD->DC_DetectorNumber[i] == 2){
        ZDD_DC_Y = Drift_Speed*(ZDD->DC_DriftTime[i]*microsecond) - Drift_Chamber_Width/2;
      }
    }

    // Faut encore réfléchir aux angles, y'a un truc pas net, l'angle après la cible est bcp trop large
    TVector3 AfterTargetDir_X(ZDD_DC_X - xtarget, 0,ZDD_R + ZDD->Get_m_Drift_Chamber_Z()[0] + ZDD->Get_m_Drift_Chamber_Thickness()[0]*0.5 - ztarget);
    ZDD_ThetaAfterTarget_X = AfterTargetDir_X.Angle(ZDir);
    if(ZDD_DC_X < xtarget){
      ZDD_ThetaAfterTarget_X *= -1;
    }
    TVector3 AfterTargetDir((ZDD_R + ZDD->Get_m_Drift_Chamber_Z()[1] + ZDD->Get_m_Drift_Chamber_Thickness()[1]*0.5 -ztarget)*tan(ZDD_ThetaAfterTarget_X) - xtarget,
    ZDD_DC_Y - ytarget,
    ZDD_R + ZDD->Get_m_Drift_Chamber_Z()[1] + ZDD->Get_m_Drift_Chamber_Thickness()[1]*0.5 - ztarget);
    ZDD_ThetaAfterTarget = AfterTargetDir.Angle(ZDir);
    /*if(tan(ZDD_ThetaAfterTarget_X) < 0){
      ZDD_ThetaAfterTarget *= -1;
    }*/
    // Faire gaffe: se poser la question qui est X et qui est Y;
    
    /////////////////////// Treating Energy ///////////////////////////////
    ZDD_dE_tot = 0;
    ZDD_E_tot = 0;
    ZDD_E_Plastic = 0;
    for(int i = 0; i < ZDD->Plastic_DetectorNumber.size(); i++){
      ZDD_E_Plastic+= ZDD->Plastic_Energy[i];
    }

    for(int i = 0; i < ZDD->IC_DetectorNumber.size(); i++){
      dE_IC[ZDD->IC_DetectorNumber[i]]+= ZDD->IC_Energy[i];
    }
    
    ////////////////////////// Note on how to compute full energy recovery: Plastic, then 1 plan of Kapton, then gas, then Mylar,
    ////////////////////////// then all IC + Mylar, then again gas and Kapton, and Mylar at the entry and exit of each DC

    //ZDD_ThetaAfterTarget = 0;
    ZDD_E_tot = ZDD_E_Plastic;
    ZDD_E_tot = Heavy_IC_Windows[1].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_Entry_Exit_Thickness()[1], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Gas[1].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_m_Gas_Gap_Thickness()[1], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Mylar[ZDD->Get_AC_Z().size()-1].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[ZDD->Get_AC_Z().size()-1], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
    // Ok ça c'est bon
    // 
    for(int i = ZDD->IC_Energy.size() -1; i >= 0; i--){
      //if(ZDD->IC_Energy[i] > 5){
      ZDD_E_tot+= ZDD->IC_Energy[i];
      ZDD_dE_tot+= ZDD->IC_Energy[i];
      ZDD_E_tot = Heavy_IC_Mylar[i+Nb_Mylar_Before_IC].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[i+Nb_Mylar_Before_IC], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
      //}
    }
    ZDD_E_tot = Heavy_IC_Gas[0].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_m_Gas_Gap_Thickness()[0], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Windows[0].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_Entry_Exit_Thickness()[0], ZDD_ThetaIC + ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Mylar[3].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[3], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_DC_Gas[1].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_m_Drift_Chamber_Thickness()[1], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Mylar[2].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[2], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Mylar[1].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[1], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_DC_Gas[0].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_m_Drift_Chamber_Thickness()[0], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_IC_Mylar[0].EvaluateInitialEnergy(ZDD_E_tot, ZDD->Get_AC_Thickness()[0], ZDD_ThetaAfterTarget);
    ZDD_E_tot = Heavy_Target.EvaluateInitialEnergy(ZDD_E_tot, TargetThickness/2 , ZDD_ThetaAfterTarget);

    //ZDD_E_tot = Beam_Target.EvaluateInitialEnergy(ZDD_E_tot, TargetThickness/2 , ZDD_ThetaAfterTarget);

    // We reproduce the succession of materials crossed

    //////////////////// MUST2 Part ////////////////////

    for(unsigned int countMust2 = 0 ; countMust2 < M2->Si_E.size() ; countMust2++){
      /************************************************/
      //Part 0 : Get the usefull Data
      // MUST2
      int TelescopeNumber = M2->TelescopeNumber[countMust2];

      /************************************************/
      // Part 1 : Impact Angle
      ThetaM2Surface = 0;
      ThetaNormalTarget = 0;
      TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
      M2_ThetaLab = HitDirection.Angle( BeamDirection );

      M2_X = M2 -> GetPositionOfInteraction(countMust2).X();
      M2_Y = M2 -> GetPositionOfInteraction(countMust2).Y();
      M2_Z = M2 -> GetPositionOfInteraction(countMust2).Z();

      ThetaM2Surface = HitDirection.Angle(- M2 -> GetTelescopeNormal(countMust2) );
      ThetaNormalTarget = HitDirection.Angle( TVector3(0,0,1) ) ;

      /************************************************/

      /************************************************/
      // Part 2 : Impact Energy
      Energy = M2_ELab = 0;
      Si_E_M2 = M2->Si_E[countMust2];
      CsI_E_M2= M2->CsI_E[countMust2];

      // if CsI
      if(CsI_E_M2>0 ){
        // The energy in CsI is calculate form dE/dx Table because
        Energy = CsI_E_M2;
        Energy = LightAl.EvaluateInitialEnergy( Energy ,0.4*micrometer , ThetaM2Surface);
        Energy+=Si_E_M2;
      }

      else
        Energy = Si_E_M2;


      // Evaluate energy using the thickness
      M2_ELab = LightAl.EvaluateInitialEnergy( Energy ,0.4*micrometer , ThetaM2Surface);
      // Target Correction
      M2_ELab   = LightTarget.EvaluateInitialEnergy( M2_ELab ,TargetThickness*0.5, ThetaNormalTarget);
      M2_dE = Si_E_M2;

      /************************************************/

      /************************************************/
      // Part 3 : Excitation Energy Calculation
      //std::cout << reaction->GetBeamEnergy() << std::endl;
      M2_Ex = reaction->ReconstructRelativistic( M2_ELab , M2_ThetaLab );
      reaction->SetBeamEnergy(InitialConditions->GetIncidentFinalKineticEnergy());
      M2_ExNoBeam=reaction->ReconstructRelativistic( M2_ELab , M2_ThetaLab );
      reaction->SetBeamEnergy(FinalBeamEnergy);
      M2_ExNoProton = reaction->ReconstructRelativistic( ReactionConditions->GetKineticEnergy(0) , ReactionConditions->GetParticleDirection(0).Angle(TVector3(0,0,1)) );
      M2_ThetaLab=M2_ThetaLab/deg;

      /************************************************/

      /************************************************/
      // Part 4 : Theta CM Calculation
      M2_ThetaCM  = reaction->EnergyLabToThetaCM( M2_ELab , M2_ThetaLab)/deg;
      /************************************************/
    }//end loop MUST2



  }
  //ZDD->Clear();
}

void Analysis::SetParticles(){
  BeamPart->SetUp(BeamName);
  BeamPart->SetKineticEnergy(InitialConditions->GetIncidentInitialKineticEnergy());
  Beta = BeamPart->GetBeta();
  Gamma = BeamPart->GetGamma();
  Velocity = BeamPart->GetVelocity();
}

bool Analysis::CheckGoodEvent(){

  return //CheckIC() && CheckPlastics() && CheckDC();  
  true;
}

bool Analysis::CheckIC(){

  std::vector<int> IC_vec(ZDD->GetICcounter(), 0);
  bool result = true; 
  for(int i = 0; i < ZDD->IC_DetectorNumber.size(); i++){
    if(ZDD->IC_Energy[i] > 0){
      IC_vec[ZDD->IC_DetectorNumber[i]]++; 
    }

  }
  for(int i = 0; i < ZDD->GetICcounter(); i++){
    if(IC_vec[i] != 1){
      result = false; 
    }
  }
  return result;
}
bool Analysis::CheckPlastics(){
  int counter = 0;
  for(int i = 0; i < ZDD->Plastic_DetectorNumber.size(); i++){
    if(ZDD->Plastic_Energy[i] > 0){
      counter++;
    }
  }
  if(counter == 1){
    return true;
  }
  else{
    return false;
  }  
}

bool Analysis::CheckDC(){
  int DC1_counter = 0, DC2_counter = 0;
  for(int i = 0; i < ZDD->DC_DetectorNumber.size(); i++){
    if(ZDD->DC_DetectorNumber[i] == 1){
      DC1_counter++;
    }
    else if(ZDD->DC_DetectorNumber[i] == 2){
      DC2_counter++;
    }
  }
  return DC2_counter == 1 && DC1_counter == 1;
}

void Analysis::InitOutputBranch() {
  /////////////////////// ZDD related branch ////////////////////////////////
  RootOutput::getInstance()->GetTree()->Branch("ZDD_E_tot",&ZDD_E_tot,"ZDD_E_tot/D");
  //RootOutput::getInstance()->GetTree()->Branch("ZDD_IC_Energy",&ZDD_IC_Energy,"ZDD_IC_Energy[10]/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_dE_tot",&ZDD_dE_tot,"ZDD_dE_tot/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_E_Plastic",&ZDD_E_Plastic,"ZDD_E_Plastic/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_DC_X",&ZDD_DC_X,"ZDD_DC_X/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_DC_Y",&ZDD_DC_Y,"ZDD_DC_Y/D");
  RootOutput::getInstance()->GetTree()->Branch("xtarget",&xtarget,"xtarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ytarget",&ytarget,"ytarget/D");
  RootOutput::getInstance()->GetTree()->Branch("Beta",&Beta,"Beta/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_ThetaAfterTarget",&ZDD_ThetaAfterTarget,"ZDD_ThetaAfterTarget/D");
  RootOutput::getInstance()->GetTree()->Branch("ZDD_ThetaAfterTarget_X",&ZDD_ThetaAfterTarget_X,"ZDD_ThetaAfterTarget_X/D");
  RootOutput:: getInstance()->GetTree()->Branch("InitialConditions",&InitialConditions);
  
  //////////////////////// MUST related branch ////////////////////////////////
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex",&M2_Ex,"M2_Ex/D");
  //RootOutput::getInstance()->GetTree()->Branch("ExNoBeam",&ExNoBeam,"ExNoBeam/D");
  //RootOutput::getInstance()->GetTree()->Branch("ExNoProton",&ExNoProton,"ExNoProton/D");
  //RootOutput::getInstance()->GetTree()->Branch("EDC",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_ELab",&M2_ELab,"M2_ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaLab",&M2_ThetaLab,"M2_ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_ThetaCM",&M2_ThetaCM,"M2_ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_X",&M2_X,"M2_X/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_Y",&M2_Y,"M2_Y/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_Z",&M2_Z,"M2_Z/D");
  RootOutput::getInstance()->GetTree()->Branch("M2_dE",&M2_dE,"M2_dE/D");
}

void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&InitialConditions);
  RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
  RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",&ReactionConditions);
}
////////////////////////////////////////////////////////////////////////////////

void Analysis::ReInit(){
  ZDD_E_tot = -1000;
  ZDD_E_Plastic = -1000;
  ZDD_dE_tot= -1000;
  ZDD_DC_X = -1000;
  ZDD_DC_Y = -1000;
  xtarget = -1000;
  ytarget = -1000;
  Beta = -1000;
  
  M2_Ex = -1000 ;
  //ExNoBeam=ExNoProton=-1000;
  //EDC= -1000;
  M2_ELab = -1000;
  //BeamEnergy = -1000;
  M2_ThetaLab = -1000;
  M2_ThetaCM = -1000;
  M2_X = -1000;
  M2_Y = -1000;
  M2_Z = -1000;
  M2_dE= -1000;
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

