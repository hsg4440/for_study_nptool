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
   
  reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction->GetBeamEnergy();


  string Path = "../../Inputs/EnergyLoss/";
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string beam=  NPL::ChangeNameToG4Standard(reaction->GetNucleus1()->GetName());
  string heavy_ejectile=  NPL::ChangeNameToG4Standard(reaction->GetNucleus4()->GetName());
  string light=NPL::ChangeNameToG4Standard(reaction->GetNucleus3()->GetName());



//
  ProtonSi = NPL::EnergyLoss(Path+ "proton_Si.G4table", "G4Table", 100);

  Cal = CalibrationManager::getInstance();  
}
  ///////////////////////////// Initialize some important parameters //////////////////////////////////


bool Analysis::UnallocateBeforeBuild(){
  // std::cout << "test unallocate" << std::endl;
  GATCONFMASTER = **GATCONFMASTER_;
  return (GATCONFMASTER > 0); 
  //return true;
}

bool Analysis::UnallocateBeforeTreat(){
  for(int i = 0; i < 10; i++){
  PlasticRaw[i] = (*PlasticRaw_   )[i];
  PlasticRawTS[i] = (*PlasticRaw_TS_)[i];
  }
  
  for(int i = 0; i < 6; i++){
  IC_ZDDRaw[i] =  (*IC_ZDDRaw_)[i];
  IC_ZDDRawTS[i] =  (*IC_ZDDRaw_TS_)[i];
  }
  
  TAC_CATS_PL = **TAC_CATS_PL_;
  TAC_CATS_PLTS = **TAC_CATS_PL_TS_;
  
  TAC_CATS_HF = **TAC_CATS_HF_;
  TAC_CATS_HFTS = **TAC_CATS_HF_TS_;
  
  TAC_CATS_EXOGAM = **TAC_CATS_EXOGAM_;
  TAC_CATS_EXOGAMTS = **TAC_CATS_EXOGAM_TS_;
  
  TAC_MMG_CATS2 = **TAC_MMG_CATS2_;
  TAC_MMG_CATS2TS = **TAC_MMG_CATS2_TS_;
  
  TAC_MMG_CATS1 = **TAC_MMG_CATS1_;
  TAC_MMG_CATS1TS = **TAC_MMG_CATS1_TS_;
  
  TAC_MMG_EXOGAM = **TAC_MMG_EXOGAM_;
  TAC_MMG_EXOGAMTS = **TAC_MMG_EXOGAM_TS_;
  
  TAC_CATS1_CATS2 = **TAC_CATS1_CATS2_;
  TAC_CATS1_CATS2TS = **TAC_CATS1_CATS2_TS_;
  
  TAC_D4_CATS1 =**TAC_D4_CATS1_;
  TAC_D4_CATS1 =**TAC_D4_CATS1_TS_;
  
  TAC_PL_1 = **TAC_PL_1_;
  TAC_PL_1TS = **TAC_PL_1_TS_;
  TAC_PL_2 = **TAC_PL_2_;
  TAC_PL_2TS = **TAC_PL_2_TS_;
  TAC_PL_3 = **TAC_PL_3_;
  TAC_PL_3TS = **TAC_PL_3_TS_;
  TAC_PL_4 = **TAC_PL_4_;
  TAC_PL_4TS = **TAC_PL_4_TS_;
  TAC_PL_5 = **TAC_PL_5_;
  TAC_PL_5TS = **TAC_PL_5_TS_;
  return true;
}

bool Analysis::FillOutputCondition(){
  return true;
  //return (CATS->MapX.size() > 0 && CATS->MapY.size()>0);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

    //  if(M2->CsI_E.size() > 0) 
    //  std::cout << "Analysis test " << M2->CsI_E[0] << " " << M2->CsI_E.size() << " " << "\n \n";

    ReInit();
    // std::cout << CATS->PositionX.size() << std::endl;
    //////////////////// MUST2 Part ////////////////////
    int M2_size = M2->Si_E.size();
    if(CATS->PositionX.size() == 2 && CATS->PositionY.size() == 2){
    for(unsigned int countMust2 = 0 ; countMust2 < M2_size ; countMust2++){
      M2_TelescopeM++;
      // MUST2
      int TelescopeNumber = M2->TelescopeNumber[countMust2];

      // Part 1 : Impact Angle
      ThetaM2Surface = 0;
      ThetaNormalTarget = 0;
      
      BeamImpact = TVector3(CATS->PositionOnTargetX,CATS->PositionOnTargetY,0); 
      // BeamImpact = TVector3(0,0,0); 
      
      BeamDirection = TVector3(CATS->PositionX[1] - CATS->PositionX[0],CATS->PositionY[1] - CATS->PositionY[0],CATS->PositionZ[1] - CATS->PositionZ[0]);
      // std::cout << CATS->PositionX[0] - CATS->PositionX[1] << " " << CATS->PositionY[0] - CATS->PositionY[1] << " " << CATS->PositionZ[0] - CATS->PositionZ[1] << std::endl;
      // BeamDirection = TVector3(0,0,1);
      
      TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
      M2_ThetaLab.push_back(HitDirection.Angle( BeamDirection ));

      //std::cout <<  M2 -> GetPositionOfInteraction(countMust2).X() << "  " << M2 -> GetPositionOfInteraction(countMust2).Y() << "  "<< M2 -> GetPositionOfInteraction(countMust2).Z() <<  std::endl;
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
      for(unsigned int i = 0; i < ParticleType.size(); i++){
        Energy[ParticleType[i]] = 0;
        CsI_Energy[ParticleType[i]] = 0;

      if(CsI_E_M2>8192 ){
        // The energy in CsI is calculate form dE/dx Table because
        CsI_Energy[ParticleType[i]] =  Cal->ApplyCalibration("MUST2/"+ParticleType[i]+"_T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_E",CsI_E_M2);
        Energy[ParticleType[i]] = CsI_Energy[ParticleType[i]];
        //std::cout << ParticleType[i]+"MUST2/T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_E" << " " <<  Energy[ParticleType[i]] << "\n";
        //Energy = LightAl.EvaluateInitialEnergy( Energy ,0.4*micrometer , ThetaM2Surface);
        //Energy+=Si_E_M2;
      }

//      else
      Energy[ParticleType[i]] += Si_E_M2;


      }
      // Evaluate energy using the thickness
      //Energy = LightAl.EvaluateInitialEnergy( Energy ,0.4*micrometer , ThetaM2Surface);
      // Target Correction
      //M2_ELab.push_back(LightTarget.EvaluateInitialEnergy(Energy ,TargetThickness*0.5, ThetaNormalTarget));
      M2_ELab.push_back(Energy["proton"]);
      M2_dE.push_back(Si_E_M2);

      // Part 3 : Excitation Energy Calculation
      M2_Ex_p.push_back(reaction->ReconstructRelativistic( Energy["proton"] , M2_ThetaLab[countMust2] ));
      M2_Ex_d.push_back(reaction->ReconstructRelativistic( Energy["deuteron"] , M2_ThetaLab[countMust2] ));
      M2_Ex_t.push_back(reaction->ReconstructRelativistic( Energy["triton"] , M2_ThetaLab[countMust2] ));
      M2_Ex_a.push_back(reaction->ReconstructRelativistic( Energy["alpha"] , M2_ThetaLab[countMust2] ));
      
      M2_CsI_E_p.push_back(CsI_Energy["proton"]);
      M2_CsI_E_d.push_back(CsI_Energy["deuteron"]);
      M2_CsI_E_t.push_back(CsI_Energy["triton"]);
      M2_CsI_E_a.push_back(CsI_Energy["alpha"]);
      
      M2_ThetaLab[countMust2]=M2_ThetaLab[countMust2]/deg;

      // Part 4 : Theta CM Calculation
      M2_ThetaCM.push_back(reaction->EnergyLabToThetaCM( M2_ELab[countMust2] , M2_ThetaLab[countMust2])/deg);

      //if(proton_cut[TelescopeNumber]){
      //  std::cout << "Test :" << proton_cut[TelescopeNumber] << " \n";
      //  
    }//end loop MUST2
    }

  
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
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_p",&M2_CsI_E_p);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_d",&M2_CsI_E_d);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_t",&M2_CsI_E_t);
  RootOutput::getInstance()->GetTree()->Branch("M2_CsI_E_a",&M2_CsI_E_a);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_p",&M2_Ex_p);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_d",&M2_Ex_d);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_t",&M2_Ex_t);
  RootOutput::getInstance()->GetTree()->Branch("M2_Ex_a",&M2_Ex_a);
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
  RootOutput::getInstance()->GetTree()->Branch("M2_ECsI_from_deltaE",&M2_ECsI_from_deltaE);
  RootOutput::getInstance()->GetTree()->Branch("GATCONF",&GATCONFMASTER);
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_PL",&TAC_CATS_PL,"TAC_CATS_PL/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_PLTS",&TAC_CATS_PLTS,"TAC_CATS_PLTS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_HF",&TAC_CATS_HF,"TAC_CATS_HF/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_HFTS",&TAC_CATS_HFTS,"TAC_CATS_HFTS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_EXOGAM",&TAC_CATS_EXOGAM,"TAC_CATS_EXOGAM/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS_EXOGAMTS",&TAC_CATS_EXOGAMTS,"TAC_CATS_EXOGAMTS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_CATS2",&TAC_MMG_CATS2,"TAC_MMG_CATS2/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_CATS2TS",&TAC_MMG_CATS2TS,"TAC_MMG_CATS2TS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_CATS1",&TAC_MMG_CATS1,"TAC_MMG_CATS1/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_CATS1TS",&TAC_MMG_CATS1TS,"TAC_MMG_CATS1TS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_EXOGAM",&TAC_MMG_EXOGAM,"TAC_MMG_EXOGAM/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_MMG_EXOGAMTS",&TAC_MMG_EXOGAMTS,"TAC_MMG_EXOGAMTS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS1_CATS2",&TAC_CATS1_CATS2,"TAC_CATS1_CATS2/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_CATS1_CATS2TS",&TAC_CATS1_CATS2TS,"TAC_CATS1_CATS2TS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_D4_CATS1",&TAC_D4_CATS1,"TAC_D4_CATS1/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_D4_CATS1TS",&TAC_D4_CATS1TS,"TAC_D4_CATS1TS/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_1",&TAC_PL_1,"TAC_PL_1/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_1TS",&TAC_PL_1TS,"TAC_PL_1TS/l");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_2",&TAC_PL_2,"TAC_PL_2/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_2TS",&TAC_PL_2TS,"TAC_PL_2TS/l");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_3",&TAC_PL_3,"TAC_PL_3/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_3TS",&TAC_PL_3TS,"TAC_PL_3TS/l");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_4",&TAC_PL_4,"TAC_PL_4/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_4TS",&TAC_PL_4TS,"TAC_PL_4TS/l");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_5",&TAC_PL_5,"TAC_PL_5/s");
  RootOutput:: getInstance()->GetTree()->Branch("TAC_PL_5TS",&TAC_PL_5TS,"TAC_PL_5TS/l");
  RootOutput:: getInstance()->GetTree()->Branch("PlasticRaw",PlasticRaw,"PlasticRaw[10]/s");
  RootOutput:: getInstance()->GetTree()->Branch("PlasticRawTS",PlasticRawTS,"PlasticRawTS[10]/l");
  
  RootOutput:: getInstance()->GetTree()->Branch("IC_ZDDRaw",IC_ZDDRaw,"IC_ZDDRaw[6]/s");
  RootOutput:: getInstance()->GetTree()->Branch("IC_ZDDRawTS",IC_ZDDRawTS,"IC_ZDDRawTS[6]/l");
  
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){
  TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
   GATCONFMASTER_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"GATCONFMASTER");
  //DATATRIG_CATS_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"DATATRIG_CATS");
  PlasticRaw_   = new TTreeReaderArray<UShort_t>(*inputTreeReader,"PlasticRaw");
  PlasticRaw_TS_ = new TTreeReaderArray<ULong64_t>(*inputTreeReader,"PlasticRawTS");
  
  IC_ZDDRaw_ = new TTreeReaderArray<UShort_t>(*inputTreeReader,"IC_ZDDRaw");
  IC_ZDDRaw_TS_= new TTreeReaderArray<ULong64_t>(*inputTreeReader,"IC_ZDDRawTS");
  
  TAC_CATS_PL_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_CATS_PL");
  TAC_CATS_PL_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_CATS_PLTS");
  
  TAC_CATS_HF_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_CATS_HF");
  TAC_CATS_HF_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_CATS_HFTS");
  
  TAC_CATS_EXOGAM_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_CATS_EXOGAM");
  TAC_CATS_EXOGAM_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_CATS_EXOGAMTS");
  
  TAC_MMG_CATS2_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_MMG_CATS2");
  TAC_MMG_CATS2_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_MMG_CATS2TS");
  
  TAC_MMG_CATS1_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_MMG_CATS1");
  TAC_MMG_CATS1_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_MMG_CATS1TS");
  
  TAC_MMG_EXOGAM_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_MMG_EXOGAM");
  TAC_MMG_EXOGAM_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_MMG_EXOGAMTS");
  
  TAC_CATS1_CATS2_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_CATS2_CATS1");
  TAC_CATS1_CATS2_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_CATS2_CATS1TS");
  
  TAC_D4_CATS1_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_D4_CATS1");
  TAC_D4_CATS1_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_D4_CATS1TS");
  
  TAC_PL_1_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_PL_1");
  TAC_PL_1_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_PL_1TS");
  TAC_PL_2_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_PL_2");
  TAC_PL_2_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_PL_2TS");
  TAC_PL_3_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_PL_3");
  TAC_PL_3_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_PL_3TS");
  TAC_PL_4_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_PL_4");
  TAC_PL_4_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_PL_4TS");
  TAC_PL_5_= new TTreeReaderValue<UShort_t>(*inputTreeReader,"TAC_PL_5");
  TAC_PL_5_TS_= new TTreeReaderValue<ULong64_t>(*inputTreeReader,"TAC_PL_5TS");
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
  
  M2_CsI_E_p.clear();
  M2_CsI_E_d.clear();
  M2_CsI_E_t.clear();
  M2_CsI_E_a.clear();
  
  M2_ECsI_from_deltaE.clear();
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

