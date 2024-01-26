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
}
  ///////////////////////////// Initialize some important parameters //////////////////////////////////


bool Analysis::UnallocateBeforeBuild(){
  GATCONFMASTER.clear();
  GATCONFMASTERTS.clear();

  GATCONFMASTER = **GATCONFMASTER_;
  return (GATCONFMASTER.size() == 1 && GATCONFMASTER[0] > 0); 
  // return true;
}

bool Analysis::UnallocateBeforeTreat(){
  GATCONFMASTERTS = **GATCONFMASTERTS_;
  return true;
}

bool Analysis::FillOutputCondition(){
  return bCATS;
  //return true;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

    ReInit();
    // std::cout << "TEST " << GATCONFMASTER.size() << std::endl;
    
    //////////////////// MUST2 Part ////////////////////
    TreatCATS();
    if(bCATS){
      TreatMUST2();
      TreatEXO();
    }
    //if(bCATS){
    //  TreatZDD();
    //  TreatTAC();
    //}
  /*//for(unsigned int countMust2 = 0 ; countMust2 < M2->Si_E.size() ; countMust2++){
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
*/
}

void Analysis::TreatCATS(){
  if(CATS->PositionOnTargetX > -1000 && CATS->PositionOnTargetY > -1000){
    BeamImpact = TVector3(CATS->PositionOnTargetX,CATS->PositionOnTargetY,0); 
    BeamDirection = TVector3(CATS->PositionX[0] - CATS->PositionX[1],CATS->PositionY[0] - CATS->PositionY[1],CATS->PositionZ[0] - CATS->PositionZ[1]);
    bCATS = true;
  }
  else bCATS = false;
  
  // BeamImpact = TVector3(0,0,0); 
  // std::cout << "Position On target : " << CATS->PositionOnTargetX << " " << CATS->PositionOnTargetY << std::endl; 
  // BeamDirection = TVector3(0,0,1);
  // std::cout << "Position XY " <<  CATS->PositionX[1] - CATS->PositionX[0] << " " << CATS->PositionY[1] - CATS->PositionY[0] << " " << CATS->PositionZ[1] - CATS->PositionZ[0] << std::endl;
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
    
    for(unsigned int i = 0; i < ParticleType.size(); i++){
      Energy[ParticleType[i]] = 0;
      CsI_Energy[ParticleType[i]] = 0;

    if(M2->CsI_E_Raw[countMust2] > 8192){
      // The energy in CsI is calculate form dE/dx Table because
      std::string name = "MUST2/"+ParticleType[i]+"_T"+NPL::itoa(TelescopeNumber)+"_CsI"+NPL::itoa(CristalNb)+"_E";
      CsI_Energy[ParticleType[i]] =  Cal->ApplyCalibration(name,M2->CsI_E_Raw[countMust2]);
      Energy[ParticleType[i]] = CsI_Energy[ParticleType[i]];
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

    // Part 3 : Excitation Energy Calculation
    M2_Ex_p.push_back(reaction->ReconstructRelativistic( Energy["proton"] , M2_ThetaLab[countMust2] ));
    M2_Ex_d.push_back(Reaction_pd->ReconstructRelativistic( Energy["deuteron"] , M2_ThetaLab[countMust2] ));
    M2_Ex_t.push_back(Reaction_pt->ReconstructRelativistic( Energy["triton"] , M2_ThetaLab[countMust2] ));
    M2_Ex_a.push_back(Reaction_p3He->ReconstructRelativistic( Energy["alpha"] , M2_ThetaLab[countMust2] ));
    
    TLorentzVector PHeavy_pd = Reaction_pd->LorentzAfterReaction(Energy["deuteron"] , M2_ThetaLab[countMust2]);
    TLorentzVector PHeavy_pt = Reaction_pd->LorentzAfterReaction(Energy["triton"] , M2_ThetaLab[countMust2]);
    TLorentzVector PHeavy_p3He = Reaction_pd->LorentzAfterReaction(Energy["alpha"] , M2_ThetaLab[countMust2]);
    Beta_pd.push_back(PHeavy_pd.Beta());
    Beta_pt.push_back(PHeavy_pt.Beta());
    Beta_p3He.push_back(PHeavy_p3He.Beta());


    M2_CsI_E_p.push_back(CsI_Energy["proton"]);
    M2_CsI_E_d.push_back(CsI_Energy["deuteron"]);
    M2_CsI_E_t.push_back(CsI_Energy["triton"]);
    M2_CsI_E_a.push_back(CsI_Energy["alpha"]);
    
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
  // }
  /*
  if(M2->Si_E.size() > 0){
    if(Inner6MVM > 0){
      
      int ExoMultMax = 2;
      if(OutersVM < 2){
        ExoMultMax == OutersVM;
      }
      for(int ExoMult = 0; ExoMult < ExoMultMax;ExoMult++){
        ExogamDetNb[ExoMult] = (OutersVN[ExoMult] - 32)/16;
        CristalNb[ExoMult] = (OutersVN[ExoMult] - (2+ ExogamDetNb[ExoMult])*16)/4;
        SegmentNb[ExoMult] = (OutersVN[ExoMult] - 16*(ExogamDetNb[ExoMult] + 2) - 4*(CristalNb[ExoMult]));
      }
      
      if(Inner6MVM == 1){
        if(OutersVM > 0 && BGOV[0] < 20){
            Theta_seg = Exogam_Clovers_struc[ExogamDetNb[0]].Theta_Crystal_Seg[CristalNb[0]][SegmentNb[0]];
            Phi_seg = Exogam_Clovers_struc[ExogamDetNb[0]].Phi_Crystal_Seg[CristalNb[0]][SegmentNb[0]];
            EnergyDoppler = Doppler_Correction(Theta_seg,Phi_seg,0,0,Beta,Inner6MV[0]);
            EnergyAddBackDoppler = EnergyDoppler;
          }
      }
      if(Inner6MVM == 2 && OutersVM > 1 && BGOV[0] < 20 && BGOV[1] < 20){
        if(Inner6MVN[0]/4 == Inner6MVN[1]/4){
          EnergyAddBack = Inner6MV[0] + Inner6MV[1];
        }
      }
      if(EnergyAddBack > -1000){
        for(int i = 0;i < ExoMultMax; i++){
          if(OutersV[i] > OutersV[highest_E]){
            highest_E = i;
          }
        }
        Theta_seg = Exogam_Clovers_struc[ExogamDetNb[highest_E]].Theta_Crystal_Seg[CristalNb[highest_E]][SegmentNb[highest_E]];
        Phi_seg = Exogam_Clovers_struc[ExogamDetNb[highest_E]].Phi_Crystal_Seg[CristalNb[highest_E]][SegmentNb[highest_E]];
        EnergyAddBackDoppler = Doppler_Correction(Theta_seg,Phi_seg,0,0,Beta,EnergyAddBack);
      }
    } 
  }
  
  
*/ 
}


void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("GATCONF",&GATCONFMASTER);
  RootOutput::getInstance()->GetTree()->Branch("GATCONFTS",&GATCONFMASTERTS);
  
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
  // RootOutput::getInstance()->GetTree()->Branch("M2_dE",&M2_dE);
  // RootOutput::getInstance()->GetTree()->Branch("CsI_E_M2",&CsI_E_M2);
  
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_pd",&EXO_Doppler_pd);
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_pt",&EXO_Doppler_pt);
  RootOutput::getInstance()->GetTree()->Branch("EXO_Doppler_p3He",&EXO_Doppler_p3He);
  
  RootOutput::getInstance()->GetTree()->Branch("Beta_pd",&Beta_pd);
  RootOutput::getInstance()->GetTree()->Branch("Beta_p3He",&Beta_p3He);
  RootOutput::getInstance()->GetTree()->Branch("Beta_pt",&Beta_pt);
  /*
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
  
  RootOutput:: getInstance()->GetTree()->Branch("EnergyDoppler",&EnergyDoppler,"EnergyDoppler/F");
  RootOutput:: getInstance()->GetTree()->Branch("EnergyAddBack",&EnergyAddBack,"EnergyAddBack/F");
  RootOutput:: getInstance()->GetTree()->Branch("EnergyAddBackDoppler",&EnergyAddBackDoppler,"EnergyAddBackDoppler/F");
  
  RootOutput:: getInstance()->GetTree()->Branch("Inner6MVM",&Inner6MVM,"Inner6MVM/I");
  RootOutput:: getInstance()->GetTree()->Branch("Inner6MV",Inner6MV,"Inner6MV[Inner6MVM]/F");
  RootOutput:: getInstance()->GetTree()->Branch("Inner6MVN",Inner6MVN,"Inner6MVN[Inner6MVM]/s");
  RootOutput:: getInstance()->GetTree()->Branch("Inner6MVTS",Inner6MVTS,"Inner6MVTS[Inner6MVM]/l");
  RootOutput:: getInstance()->GetTree()->Branch("BGOVM",&BGOVM,"BGOVM/I");
  RootOutput:: getInstance()->GetTree()->Branch("BGOV",BGOV,"BGOV[BGOVM]/F");
  RootOutput:: getInstance()->GetTree()->Branch("BGOVN",BGOVN,"BGOVN[BGOVM]/s"); 
  RootOutput:: getInstance()->GetTree()->Branch("DeltaTVM",&DeltaTVM,"DeltaTVM/I");
  RootOutput:: getInstance()->GetTree()->Branch("DeltaTV",DeltaTV,"DeltaTV[DeltaTVM]/F");
  RootOutput:: getInstance()->GetTree()->Branch("DeltaTVN",DeltaTVN,"DeltaTVN[DeltaTVM]/s");
  RootOutput:: getInstance()->GetTree()->Branch("DeltaTVTS",DeltaTVTS,"DeltaTVTS[DeltaTVM]/l");
*/
}

void Analysis::UnallocateVariables(){
}

void Analysis::InitInputBranch(){

  if(!NPOptionManager::getInstance()->GetInputPhysicalTreeOption()){
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    GATCONFMASTER_ = new TTreeReaderValue<vector<unsigned int>>(*inputTreeReader,"GATCONF");
    GATCONFMASTERTS_ = new TTreeReaderValue<vector<unsigned long long>>(*inputTreeReader,"GATCONFTS");

  }
 // else{
 //   RootInput::getInstance()->GetChain()->SetBranchAddress("GATCONF",&GATCONFMASTER);
 //   RootInput::getInstance()->GetChain()->SetBranchAddress("GATCONFTS",&GATCONFMASTERTS);
 // }

  //DATATRIG_CATS_ = new TTreeReaderValue<unsigned short>(*inputTreeReader,"DATATRIG_CATS");
  /*PlasticRaw_   = new TTreeReaderArray<UShort_t>(*inputTreeReader,"PlasticRaw");
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
  
  Inner6MVM_ = new TTreeReaderValue<int>(*inputTreeReader,"Inner6MRawM");
  Inner6MV_ = new TTreeReaderArray<float>(*inputTreeReader,"Inner6MRaw");
  Inner6MVN_ = new TTreeReaderArray<unsigned short>(*inputTreeReader,"Inner6MRawNr");
  //RootInput:: getInstance()->GetChain()->SetBranchAddress("Inner6MVTS",Inner6MVTS);

  OutersVM_ = new TTreeReaderValue<int>(*inputTreeReader,"OutersRawM");
  OutersV_ = new TTreeReaderArray<float>(*inputTreeReader,"OutersRaw");
  OutersVN_ = new TTreeReaderArray<unsigned short>(*inputTreeReader,"OutersRawNr");
  
  BGOVM_ = new TTreeReaderValue<int>(*inputTreeReader,"BGORawM");
  BGOV_ = new TTreeReaderArray<float>(*inputTreeReader,"BGORaw");
  BGOVN_ = new TTreeReaderArray<unsigned short>(*inputTreeReader,"BGORawNr");
*/
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

