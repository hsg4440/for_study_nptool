/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Edited for e793s                                                          *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    : May 2021                                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Class describing the property of an Analysis object                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include<iostream>
#include<string>

#include<fstream>

using namespace std;

#include "Analysis.h"
#include "NPFunction.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPOptionManager.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
 ///////////////////////////////////////////////////////////////////////////////  
  
//  if(NPOptionManager::getInstance()->HasDefinition("simulation")){
    cout << " == == == == SIMULATION == == == ==" << endl;
    isSim=true;
//  } else {
//    cout << " == == == == EXPERIMENT == == == ==" << endl;
//    isSim=false;
//  }

  agata_zShift=51*mm;
  //BrhoRef=0.65;

  if(isSim){
    Initial = new TInitialConditions();
    ReactionConditions = new TReactionConditions(); 
  }

  // initialize input and output branches
  InitOutputBranch();
  InitInputBranch();
  
  M2 = (TMust2Physics*)  m_DetectorManager -> GetDetector("M2Telescope");
  MG = (TMugastPhysics*) m_DetectorManager -> GetDetector("Mugast");
  ML = (TModularLeafPhysics*) m_DetectorManager -> GetDetector("ModularLeaf"); 
//  CATS = (TCATSPhysics*) m_DetectorManager->GetDetector("CATSDetector");

  // get reaction information
  reaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction.GetBeamEnergy();
  reaction.Print(); //TESTING PARSER

  // get beam position from .reaction file
  Beam = (NPL::Beam*) reaction.GetParticle1(); 
  XBeam = Beam->GetMeanX();
  YBeam = Beam->GetMeanY();

  cout << " ---------------- Beam Position ---------------- " << endl;
  cout << " \tX = " << XBeam << " mm\tY = " << YBeam  << " mm" << endl;
  cout << " ----------------------------------------------- " << endl;

  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  
  cout << " ---------------- Target Values ---------------- " << endl;
  cout << " \tZ = " << m_DetectorManager->GetTargetZ() << " mm\tT = " << TargetThickness  << " mm" << endl;
  cout << " ----------------------------------------------- " << endl;
  
  // energy losses
  string light = NPL::ChangeNameToG4Standard(reaction.GetParticle3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(reaction.GetParticle1()->GetName());
  LightTarget = NPL::EnergyLoss(light+"_"+TargetMaterial+".G4table","G4Table",100 );
  LightAl = NPL::EnergyLoss(light+"_Al.G4table" ,"G4Table",100);
  LightSi = NPL::EnergyLoss(light+"_Si.G4table" ,"G4Table",100);
  BeamTargetELoss = NPL::EnergyLoss(beam+"_"+TargetMaterial+".G4table","G4Table",100);

  FinalBeamEnergy = BeamTargetELoss.Slow(OriginalBeamEnergy, 0.5*TargetThickness, 0);
  reaction.SetBeamEnergy(FinalBeamEnergy); 

  cout << "Beam energy at mid-target: " << FinalBeamEnergy << endl;

  // initialize various parameters
  Rand = TRandom3();
  ParticleMult = 0;
  GammaMult = 0;
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  elab_tmp = 0;
  thetalab_tmp = 0;
  philab_tmp = 0;
  ThetaGDSurface = 0;
  X.clear();
  Y.clear();
  Z.clear();
  dE = 0;
  BeamDirection = TVector3(0,0,1);
  //nbTrack=0;
  //nbHits=0;
  //count=0;
  AHeavy=reaction.GetParticle4()->GetA();
  ALight=reaction.GetParticle3()->GetA(); 
  MHeavy=reaction.GetParticle4()->Mass();
  MLight=reaction.GetParticle3()->Mass();
  bool writetoscreen=true;

  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    GATCONF_Counter[i] = 0 ; 
  }

  ThetaCM_detected->Sumw2();
  ThetaLab_detected->Sumw2();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){

  // Reinitiate calculated variable
  ReInitValue();

  if(isSim){
    //ThetaCM_emmitted->Fill(Initial->GetThetaCM(0));
    ThetaLab_emmitted->Fill(Initial->GetThetaLab_WorldFrame(0));
    
    TVector3 labDir = Initial->GetParticleDirection(0);
    TVector3 comDir = labDir;

    comDir.SetZ( cos(0.5*Initial->GetThetaLab_WorldFrame(0)*deg) );
    comDir.SetX( labDir.X() / (2*comDir.Z()) );
    comDir.SetY( labDir.Y() / (2*comDir.Z()) );

    /*
    cout << Initial->GetThetaLab_WorldFrame(0) << "  " 
	 << labDir.Theta()/deg << " "
	 << comDir.Theta()/deg << " "
	 << endl;
    */

    ThetaCM_emmitted->Fill(comDir.Theta()/deg);
  }

  GATCONF_MASTER=ML->GetCalibratedValue("GATCONF_MASTER");
  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    if(GATCONF_MASTER & (unsigned int)pow(2,i)){ // test if ith bit is on
      GATCONF_Counter[i]++ ; // increment the array
      //      cout <<  std::bitset<16>(GATCONF_MASTER) << " " << i+1 << endl;
    }
  }

//  if(isSim){
//    OriginalELab = ReactionConditions->GetKineticEnergy(0);
//    OriginalThetaLab = ReactionConditions->GetTheta(0);
//    BeamEnergy = ReactionConditions->GetBeamEnergy();
//  }

  TVector3 BeamDirection(0.,0.,1.);
  BeamImpact = TVector3(XBeam,YBeam,m_DetectorManager->GetTargetZ()); 

  ParticleMult=M2->Si_E.size()+MG->DSSD_E.size();
  //ParticleMult=M2->Si_E.size();
  //  FinalBeamEnergy=BeamCD2.Slow(OriginalBeamEnergy,0.5*TargetThickness,BeamDirection.Angle(TVector(0,0,1)));
  //reaction.SetBeamEnergy(FinalBeamEnergy);


  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// LOOP on MUST2  ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  unsigned int sizeM2 = M2->Si_E.size();
  for(unsigned int countMust2 = 0; countMust2 < sizeM2; countMust2++){

    /************************************************/
    // Part 0 : Get the useful Data
    // MUST2
    int TelescopeNumber = M2->TelescopeNumber[countMust2];

    if(isSim){
      if(TelescopeNumber==5){
        //ThetaCM_detected->Fill(Initial->GetThetaCM(0));
	ThetaLab_detected->Fill(Initial->GetThetaLab_WorldFrame(0));

        TVector3 labDir = Initial->GetParticleDirection(0);
        TVector3 comDir = labDir;

        comDir.SetZ( cos(0.5*Initial->GetThetaLab_WorldFrame(0)*deg) );
        comDir.SetX( labDir.X() / (2*comDir.Z()) );
        comDir.SetY( labDir.Y() / (2*comDir.Z()) );

        /*
	// Check the angles agree
        cout << Initial->GetThetaLab_WorldFrame(0) << "  " 
	     << labDir.Theta()/deg << " "
	     << comDir.Theta()/deg << " "
	     << endl;
        */

        ThetaCM_detected->Fill(comDir.Theta()/deg);

      }
    }

    /************************************************/
    // Part 1 : Impact Angle
    ThetaM2Surface = 0;
    ThetaNormalTarget = 0;
    thetalab_tmp = 0;
    philab_tmp = 0;
    TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
    thetalab_tmp = HitDirection.Angle(BeamDirection);
    
    /******************************************************************************/
    /*** THIS SECTION IS FOR USE WHEN CONVERTING ISOTROPIC LAB TO ISOTROPIC CoM ***/
//      if(writetoscreen){
//        cout << "====================================" << endl;
//        cout << "   CONVERTING ISO LAB TO ISO CoM!   " << endl;
//        cout << "====================================" << endl;
//	writetoscreen=false;
//      }
      /* Generate a normalized vector/unit vector for use in theta calculation*/
//        TVector3 NormHitDir = HitDirection.Unit();

      /* thetaCM = 0.5 thetaLab = 0.5 acos(zLab) */
//        thetalab_tmp = 0.5 * acos(NormHitDir.Z());
        /* Now, thetalab_tmp is in CoM frame */

      /* zCM = cos(thetaCM) */
//        HitDirection.SetZ(cos(thetalab_tmp));
        /* Now, (xLab, yLab,  zCM) */

      /* xCM = xLab/(2 zCM) */
//        HitDirection.SetX(HitDirection.X()/(2*HitDirection.Z()));
        /* Now, (xCM,  yLab,  zCM) */

      /* yCM = yLab/(2 zCM) */
//        HitDirection.SetY(HitDirection.Y()/(2*HitDirection.Z()));
        /* Now, (xCM,  yCM,  zCM) */
    /******************************************************************************/
    /******************************************************************************/
    
    philab_tmp = HitDirection.Phi();

    X.push_back(M2->GetPositionOfInteraction(countMust2).X());
    Y.push_back(M2->GetPositionOfInteraction(countMust2).Y());
    Z.push_back(M2->GetPositionOfInteraction(countMust2).Z());

    ThetaM2Surface = HitDirection.Angle(- M2->GetTelescopeNormal(countMust2) );
    ThetaNormalTarget = HitDirection.Angle( TVector3(0.0,0.0,1) ) ;

    /************************************************/
    // Part 2 : Impact Energy
    Energy = 0;
    elab_tmp = 0;
    Si_E_M2 = M2->Si_E[countMust2];
    CsI_E_M2= M2->CsI_E[countMust2];

    // if CsI
    if(CsI_E_M2>0 ){
      // The energy in CsI is calculate form dE/dx Table because
      Energy = CsI_E_M2;
      if(!isSim){
        Energy = LightAl.EvaluateInitialEnergy(Energy, 0.4*micrometer, ThetaM2Surface);
      }
      Energy+=Si_E_M2;
    }
    else
      Energy = Si_E_M2;

    RawEnergy.push_back(Energy);

    /*
    cout << M2->TelescopeNumber[0] << "   "
         << M2->GetTelescopeNormal(countMust2).X() << "   "
         << M2->GetTelescopeNormal(countMust2).Y() << "   "
         << M2->GetTelescopeNormal(countMust2).Z() << endl;
    */

    if(!isSim){
      // Evaluate energy using the thickness
      elab_tmp = LightAl.EvaluateInitialEnergy(Energy, 0.4*micrometer, ThetaM2Surface);
      // Target Correction
      elab_tmp = LightTarget.EvaluateInitialEnergy(elab_tmp, 0.5*TargetThickness, ThetaNormalTarget);
    } else {elab_tmp = Energy;}

    ELab.push_back(elab_tmp);

    /************************************************/

    /************************************************/
    // Part 3 : Excitation Energy Calculation
    Ex.push_back(reaction.ReconstructRelativistic(elab_tmp,thetalab_tmp));
    Ecm.push_back(Energy*(AHeavy+ALight)/(4*AHeavy*cos(thetalab_tmp)*cos(thetalab_tmp)));
    /************************************************/

    /************************************************/
    // Part 4 : Theta CM Calculation
    
    ThetaCM.push_back(reaction.EnergyLabToThetaCM(elab_tmp, thetalab_tmp)/deg);
    /************************************************/

    ThetaLab.push_back(thetalab_tmp/deg);
    PhiLab.push_back(philab_tmp/deg);
  }

/*
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// LOOP on MUGAST ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  unsigned int sizeMG = MG->DSSD_E.size();
  for(unsigned int countMugast = 0; countMugast<sizeMG; countMugast++){
    // Part 1 : Impact Angle
    ThetaMGSurface = 0;
    ThetaNormalTarget = 0;
    thetalab_tmp = 0;
    philab_tmp = 0;
    TVector3 HitDirection = MG->GetPositionOfInteraction(countMugast) - BeamImpact;
    thetalab_tmp = HitDirection.Angle(BeamDirection);
    philab_tmp = HitDirection.Phi();

    X.push_back( MG -> GetPositionOfInteraction(countMugast).X());
    Y.push_back( MG -> GetPositionOfInteraction(countMugast).Y());
    Z.push_back( MG -> GetPositionOfInteraction(countMugast).Z());

    ThetaMGSurface = HitDirection.Angle( MG -> GetTelescopeNormal(countMugast) );
    ThetaNormalTarget = HitDirection.Angle( TVector3(0.0,0.0,1.0) ) ;

//cout << MG->TelescopeNumber[0] << "\t" 
//     << "\t" << MG->GetTelescopeNormal(countMugast).X()
//     << "\t" << MG->GetTelescopeNormal(countMugast).Y()
//     << "\t" << MG->GetTelescopeNormal(countMugast).Z() << endl;

    // Part 2 : Impact Energy
    Energy = elab_tmp = 0;
    Energy = MG->GetEnergyDeposit(countMugast);
    RawEnergy.push_back(Energy);

    if(!isSim){
      elab_tmp = LightAl.EvaluateInitialEnergy(
		    Energy,              //particle energy after Al
		    0.4*micrometer,      //thickness of Al
		    ThetaMGSurface);     //angle of impingement
      elab_tmp = LightTarget.EvaluateInitialEnergy(
		    elab_tmp,            //particle energy after leaving target
		    TargetThickness*0.5, //distance passed through target
		    ThetaNormalTarget);  //angle of exit from target
    } else {elab_tmp = Energy;}

    ELab.push_back(elab_tmp);

    // Part 3 : Excitation Energy Calculation
    if(!isSim){
      Ex.push_back(reaction.ReconstructRelativistic(elab_tmp,thetalab_tmp));
      Ecm.push_back(elab_tmp*(AHeavy+ALight)/(4*AHeavy*cos(thetalab_tmp)*cos(thetalab_tmp)));
    }

    // Part 4 : Theta CM Calculation
    ThetaLab.push_back(thetalab_tmp/deg);
    PhiLab.push_back(philab_tmp/deg);
    ThetaCM.push_back(reaction.EnergyLabToThetaCM(elab_tmp, thetalab_tmp)/deg);

    if(sizeMG==1){
      MG_T = MG->DSSD_T[0];
      MG_E = MG->DSSD_E[0];
      MG_X = MG->DSSD_X[0];
      MG_Y = MG->DSSD_Y[0];
      MG_D = MG->TelescopeNumber[0];
    }

  }//end loop Mugast

*/

  ////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////// LOOP on AGATA ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // Agata by track
  /*
  for(int j=0; j<nbTrack; j++){ // all multiplicity
    TLorentzVector GammaLV;
    // Measured E
    double Egamma=trackE[j]/1000.; // From keV to MeV 
    // Gamma detection position
    // TrackZ1 to be corrected there is a shift of +51mm
    TVector3 GammaHit(trackX1[j],trackY1[j],trackZ1[j]+agata_zShift); 
    // TVector3 GammaHit(trackX1[0],trackY1[0],trackZ1[0]); 
    // Gamma Direction 
    TVector3 GammaDirection = GammaHit-BeamImpact;
    GammaDirection = GammaDirection.Unit();
    // Beta from Two body kinematic
    //TVector3 beta = reaction.GetEnergyImpulsionLab_4().BoostVector();
    // Beta from the Beam mid target 
    reaction.GetKinematicLine4();
    TVector3 beta(0,0,-reaction.GetNucleus4()->GetBeta());

    double ThetaGamma = GammaDirection.Angle(BeamDirection)/deg;
    // Construct LV in lab frame
    GammaLV.SetPx(Egamma*GammaDirection.X());
    GammaLV.SetPy(Egamma*GammaDirection.Y());
    GammaLV.SetPz(Egamma*GammaDirection.Z());
    GammaLV.SetE(Egamma);
    // Boost back in CM
    GammaLV.Boost(beta);
    // Get EDC
    EDC.push_back(GammaLV.Energy());
  }
  */

/********
  // Agata add back is not always multiplicity 1 ?? NO, not necessarily!
  for(int i=0; i<nbAdd; i++){
  //if(nbAdd==1){
    TLorentzVector GammaLV;
    // Measured E
    double Egamma = AddE[i]/1000.; // From keV to MeV 
    // Gamma detection position
    // TrackZ1 to be corrected there is a shift of +51mm
    TVector3 GammaHit(AddX[i],AddY[i],AddZ[i]+agata_zShift); 
    // TVector3 GammaHit(trackX1[0],trackY1[0],trackZ1[0]); 
    // Gamma Direction 
    TVector3 GammaDirection = GammaHit-BeamImpact;
    GammaDirection = GammaDirection.Unit();
    // Beta from Two body kinematic
    //TVector3 beta = reaction.GetEnergyImpulsionLab_4().BoostVector();
    // Beta from the Beam mid target 
    reaction.GetKinematicLine4();
    // TVector3 beta(0,0,-reaction.GetNucleus4()->GetBeta());
    TVector3 beta(0,0,-0.1257);

    double ThetaGamma = GammaDirection.Angle(BeamDirection)/deg;
    // Construct LV in lab frame
    GammaLV.SetPx(Egamma*GammaDirection.X());
    GammaLV.SetPy(Egamma*GammaDirection.Y());
    GammaLV.SetPz(Egamma*GammaDirection.Z());
    GammaLV.SetE(Egamma);
    // Boost back in CM
    GammaLV.Boost(beta);
    // Get EDC
    AddBack_EDC.push_back(GammaLV.Energy());
    AddBack_EDC2.push_back(GammaLV.Energy());

  }

  float offset[21] = {0.,27.4,46.7,41.1,71.2,58.7,32.1,40.4,-46.3,15.1,9.9,-39.1,-39.8,-5.3,18.5,0.7,-28.5,-7,-22,-8.5,-8.5};
  float slope[21] = {1.,0.96,0.93,0.93,0.89,0.91,0.95,0.93,1.06,0.97,0.99,1.05,1.04,1.01,0.97,1.,1.03,0.99,1.,1.02,1.02};

  for(unsigned int i=0; i<MW_Nr; i++){
    MWT[i] = -1500;
    for(unsigned short j=0; j<21; j++)
      if(MW_N[i]==j)
        MWT[i] = (MW_T[i]-offset[j])/slope[j];
  }

********/
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
  cout << endl << "\e[1;32m GATCONF Statistics " << endl ;
  for(int i=0;i<GATCONF_SIZE;i++){ // loop over the bits
    cout << GATCONF_Label[i]  << "\t\t" << GATCONF_Counter[i] << endl ; 
  }
  cout << endl ;

  if(isSim){

    TObjArray HistList(0);

    auto Efficiency_CM = new TH1F(*ThetaCM_detected);
    //auto Efficiency_CM = new TH1F(*ThetaCM_emmitted);
    Efficiency_CM->SetName("Efficiency_CM");
    Efficiency_CM->SetTitle("Efficiency_CM");
    Efficiency_CM->Sumw2();
    auto Efficiency_Lab = new TH1F(*ThetaLab_detected);
    //auto Efficiency_Lab = new TH1F(*ThetaLab_emmitted);
    Efficiency_Lab->SetName("Efficiency_Lab");
    Efficiency_Lab->SetTitle("Efficiency_Lab");
    Efficiency_Lab->Sumw2();
     
    auto SolidAngle_CM = new TH1F(*ThetaCM_detected);
    //auto SolidAngle_CM = new TH1F(*ThetaCM_emmitted);
    SolidAngle_CM->SetName("SolidAngle_CM");
    SolidAngle_CM->SetTitle("SolidAngle_CM");
    auto SolidAngle_Lab = new TH1F(*ThetaLab_detected);
    //auto SolidAngle_Lab = new TH1F(*ThetaLab_emmitted);
    SolidAngle_Lab->SetName("SolidAngle_Lab");
    SolidAngle_Lab->SetTitle("SolidAngle_Lab");

    Efficiency_CM->Divide(ThetaCM_emmitted);
    //Efficiency_CM->Divide(ThetaCM_detected);
    Efficiency_Lab->Divide(ThetaLab_emmitted);
    //Efficiency_Lab->Divide(ThetaLab_detected);

    double dt = 180./Efficiency_Lab->GetNbinsX();
    cout << "Angular infinitesimal = " << dt << "deg " << endl;
    auto Cline = new TF1("Cline",Form("1./(2*%f*sin(x*%f/180.)*%f*%f/180.)",M_PI,M_PI,dt,M_PI),0,180);

    SolidAngle_CM->Divide(ThetaCM_emmitted);
    //SolidAngle_CM->Divide(ThetaCM_detected);
    SolidAngle_CM->Divide(Cline,1);
    SolidAngle_Lab->Divide(ThetaLab_emmitted);
    //SolidAngle_Lab->Divide(ThetaLab_detected);
    SolidAngle_Lab->Divide(Cline,1);

    HistList.Add(ThetaCM_emmitted);
    HistList.Add(ThetaLab_emmitted);
    HistList.Add(ThetaCM_detected);
    HistList.Add(ThetaLab_detected);
    HistList.Add(Efficiency_CM);
    HistList.Add(Efficiency_Lab);
    HistList.Add(SolidAngle_CM);
    HistList.Add(SolidAngle_Lab);
    HistList.Add(Cline);

    //HistoFile->Write();
    auto HistoFile = new TFile("SolidAngle_HistoFile.root","RECREATE");
    HistList.Write();
    HistoFile->Close();
  }
}

  ////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch(){
  //RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex);
  //RootOutput::getInstance()->GetTree()->Branch("EDC",&EDC,"EDC/D");
  RootOutput::getInstance()->GetTree()->Branch("EDC",&EDC);
  RootOutput::getInstance()->GetTree()->Branch("AddBack_EDC",&AddBack_EDC);
  RootOutput::getInstance()->GetTree()->Branch("AddBack_EDC2",&AddBack_EDC2);
  RootOutput::getInstance()->GetTree()->Branch("EAgata",&EAgata,"EAgata/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab);
  RootOutput::getInstance()->GetTree()->Branch("Ecm",&Ecm);
  RootOutput::getInstance()->GetTree()->Branch("RawEnergy",&RawEnergy); // CPx ADDITION
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("PhiLab",&PhiLab);
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM);
  RootOutput::getInstance()->GetTree()->Branch("Run",&Run,"Run/I");
  RootOutput::getInstance()->GetTree()->Branch("X",&X);
  RootOutput::getInstance()->GetTree()->Branch("Y",&Y);
  RootOutput::getInstance()->GetTree()->Branch("Z",&Z);
  RootOutput::getInstance()->GetTree()->Branch("dE",&dE,"dE/D");
  RootOutput::getInstance()->GetTree()->Branch("MG_T",MG_T);
  RootOutput::getInstance()->GetTree()->Branch("MG_E",MG_E);
  RootOutput::getInstance()->GetTree()->Branch("MG_X",MG_X);
  RootOutput::getInstance()->GetTree()->Branch("MG_Y",MG_Y);
  RootOutput::getInstance()->GetTree()->Branch("MG_D",MG_D);

  // Vamos 
  RootOutput::getInstance()->GetTree()->Branch("LTS",&LTS,"LTS/l");
/*    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS1",&T_FPMW_CATS1,"T_FPMW_CATS1/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS2_C",&T_FPMW_CATS2_C,"T_FPMW_CATS2_C/F");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF",&T_FPMW_HF,"T_FPMW_CATS1/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF_C",&T_FPMW_HF_C,"T_FPMW_CATS1/F");
    RootOutput::getInstance()->GetTree()->Branch ("IC",IC, "IC[12]/F");
    RootOutput::getInstance()->GetTree()->Branch ("ICRawPUMult01",ICRawPUMult01, "ICRawPUMult01[7]/s");

    RootOutput::getInstance()->GetTree()->Branch("T_CATS2_HF",&T_CATS2_HF,"T_CATS2_HF/s");
    RootOutput::getInstance()->GetTree()->Branch("T_MUGAST_FPMW",&T_MUGAST_FPMW,"T_MUGAST_FPMW/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_CATS2",&T_FPMW_CATS2,"T_FPMW_CATS2/s");
    RootOutput::getInstance()->GetTree()->Branch("T_FPMW_HF",&T_FPMW_HF,"T_FPMW_HF/s");

    // DC0
    RootOutput::getInstance()->GetTree()->Branch ("DC0_X", &DC0_X, "DC0_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC0_Y", &DC0_Y, "DC0_Y/F");
    // DC1
    RootOutput::getInstance()->GetTree()->Branch ("DC1_X", &DC1_X, "DC1_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC1_Y", &DC1_Y, "DC1_Y/F");
    // DC2
    RootOutput::getInstance()->GetTree()->Branch ("DC2_X", &DC2_X, "DC2_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC2_Y", &DC2_Y, "DC2_Y/F");
    // DC3
    RootOutput::getInstance()->GetTree()->Branch ("DC3_X", &DC3_X, "DC3_X/F");
    RootOutput::getInstance()->GetTree()->Branch ("DC3_Y", &DC3_Y, "DC3_Y/F");

    RootOutput::getInstance()->GetTree()->Branch ("Xf", &Xf, "Xf/F");
    RootOutput::getInstance()->GetTree()->Branch ("Tf", &Tf, "Tf/F");

    RootOutput::getInstance()->GetTree()->Branch ("Brho", &Brho, "Brho/F");
    RootOutput::getInstance()->GetTree()->Branch ("Theta", &Theta, "Theta/F");
    RootOutput::getInstance()->GetTree()->Branch ("Phi", &Phi, "Phi/F");
    RootOutput::getInstance()->GetTree()->Branch ("Path", &Path, "Path/F");

    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_1_1", &EWIRE_1_1, "EWIRE_1_1/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_1_2", &EWIRE_1_2, "EWIRE_1_2/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_2_1", &EWIRE_2_1, "EWIRE_2_1/s");
    RootOutput::getInstance()->GetTree()->Branch ("EWIRE_2_2", &EWIRE_2_2, "EWIRE_2_2/s");
*/
  RootOutput::getInstance()->GetTree()->Branch ("MW_Nr", &MW_Nr, "MW_Nr/I");
  RootOutput::getInstance()->GetTree()->Branch ("MW_T", MW_T, "MW_T[MW_Nr]/F");
  RootOutput::getInstance()->GetTree()->Branch ("MW_N", MW_N, "MW_N[MW_Nr]/s");
  RootOutput::getInstance()->GetTree()->Branch ("MWT", MWT, "MWT[MW_Nr]/F");

  RootOutput::getInstance()->GetTree()->Branch ("AGAVA_VAMOSTS", &AGAVA_VAMOSTS, "AGAVA_VAMOSTS/l");
/*    RootOutput::getInstance()->GetTree()->Branch("mT",&mT,"mT/D");
    RootOutput::getInstance()->GetTree()->Branch("mV",&mV,"mV/D");
    RootOutput::getInstance()->GetTree()->Branch("mD",&mD,"mD/D");
    RootOutput::getInstance()->GetTree()->Branch("mBeta",&mBeta,"mBeta/D");
    RootOutput::getInstance()->GetTree()->Branch("mGamma",&mGamma,"mGamma/D");
    RootOutput::getInstance()->GetTree()->Branch("mM_Q",&mM_Q,"mM_Q/D");
    RootOutput::getInstance()->GetTree()->Branch("mM",&mM,"mM/D");
    RootOutput::getInstance()->GetTree()->Branch("mE",&mE,"mE/D");
    RootOutput::getInstance()->GetTree()->Branch("mdE",&mdE,"mdE/D");
*/

    // DC0->3_X/_Y Xf, Tf, Brho, Path, EWIRE_1_1 all, MW_Nr, MW_T, AGAVA_VAMOSTS, 
    // Agata
    // Time stamp of the agata trigger
    /*
    RootOutput::getInstance()->GetTree()->Branch("TStrack",&TStrack,"TStrack/l");

    // Array of reconstructed tracks
    RootOutput::getInstance()->GetTree()->Branch("nbTrack",&nbTrack,"nbTrack/I");
    RootOutput::getInstance()->GetTree()->Branch("trackE",trackE,"trackE[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackX1",trackX1,"trackX1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackY1",trackY1,"trackY1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackZ1",trackZ1,"trackZ1[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackT",trackT,"trackT[nbTrack]/F");
    RootOutput::getInstance()->GetTree()->Branch("trackCrystalID",trackCrystalID,"trackCrystalID[nbTrack]/I");
    */
    
  //Agata Addback
  RootOutput::getInstance()->GetTree()->Branch("nbAdd",&nbAdd,"nbAdd/I");
  RootOutput::getInstance()->GetTree()->Branch("TSHit",&TSHit,"TSHit/l");
  RootOutput::getInstance()->GetTree()->Branch("AddE",AddE,"AddE[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddX",AddX,"AddX[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddY",AddY,"AddY[nbAdd]/F");
  RootOutput::getInstance()->GetTree()->Branch("AddZ",AddZ,"AddZ[nbAdd]/F");

  if(isSim){
    RootOutput::getInstance()->GetTree()->Branch("OriginalELab",
	      &OriginalELab,"OriginalELab/D");
    RootOutput::getInstance()->GetTree()->Branch("OriginalThetaLab",
	      &OriginalThetaLab,"OriginalThetaLab/D");
    RootOutput::getInstance()->GetTree()->Branch("BeamEnergy",
	      &BeamEnergy,"BeamEnergy/D");
  }
}


////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  SetBranchStatus();
  // RootInput:: getInstance()->GetChain()->SetBranchAddress("GATCONF",&vGATCONF);
  //
  // Vamos
  RootInput::getInstance()->GetChain()->SetBranchAddress("LTS",&LTS);
/*    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS1",&T_FPMW_CATS1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS2_C",&T_FPMW_CATS2_C);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_HF_C",&T_FPMW_HF_C);

    RootInput::getInstance()->GetChain()->SetBranchAddress("T_CATS2_HF",&T_CATS2_HF);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_MUGAST_FPMW",&T_MUGAST_FPMW);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_CATS2",&T_FPMW_CATS2);
    RootInput::getInstance()->GetChain()->SetBranchAddress("T_FPMW_HF",&T_FPMW_HF);

    RootInput::getInstance()->GetChain()->SetBranchAddress("Xf", &Xf);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Tf", &Tf);

    RootInput::getInstance()->GetChain()->SetBranchAddress("Brho", &Brho);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Theta", &Theta);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Phi", &Phi);
    RootInput::getInstance()->GetChain()->SetBranchAddress("Path", &Path);

    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_1_1", &EWIRE_1_1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_1_2", &EWIRE_1_2);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_2_1", &EWIRE_2_1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("EWIRE_2_2", &EWIRE_2_2);
*/
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTVM", &MW_Nr);
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTV", MW_T);
  RootInput::getInstance()->GetChain()->SetBranchAddress("MWTVN", MW_N);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AGAVA_VAMOSTS", &AGAVA_VAMOSTS);

  // Agata
     /*
    RootInput::getInstance()->GetChain()->SetBranchAddress("TStrack",&TStrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("nbTrack",&nbTrack);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackE",trackE);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackX1",trackX1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackY1",trackY1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackZ1",trackZ1);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackT",trackT);
    RootInput::getInstance()->GetChain()->SetBranchAddress("trackCrystalID",trackCrystalID);
    */
  RootInput::getInstance()->GetChain()->SetBranchAddress("TSHit",&TSHit);
  RootInput::getInstance()->GetChain()->SetBranchAddress("nbAdd",&nbAdd);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddE",AddE);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddX",AddX);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddY",AddY);
  RootInput::getInstance()->GetChain()->SetBranchAddress("AddZ",AddZ);

  if(isSim){
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&Initial);
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
    //RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchAddress("ReactionConditions",
      &ReactionConditions);
  }
}

void Analysis::SetBranchStatus(){
  // Set Branch status 
  RootInput::getInstance()->GetChain()->SetBranchStatus("LTS",true);
 /*   RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS1",true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS2_C",true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_HF",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_HF_C",true);
    
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_CATS2_HF",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_MUGAST_FPMW",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("T_FPMW_CATS2",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("Xf",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("Tf",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("Brho",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("Path",true );

    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_1_1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_1_2",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_2_1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("EWIRE_2_2",true );
*/
  RootInput::getInstance()->GetChain()->SetBranchStatus("MW_Nr",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MW_T",true  );

  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTV",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTVN",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("MWTVM",true  );

  RootInput::getInstance()->GetChain()->SetBranchStatus("AGAVA_VAMOSTS",true  );
  // Agata
    /*
    RootInput::getInstance()->GetChain()->SetBranchStatus("TStrack",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("nbTrack",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackE",true  );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackX1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackY1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackZ1",true );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackT",true  );
    RootInput::getInstance()->GetChain()->SetBranchStatus("trackCrystalID",true);
    */
  RootInput::getInstance()->GetChain()->SetBranchStatus("TSHit",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("nbAdd",true );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddE",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddX",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddY",true  );
  RootInput::getInstance()->GetChain()->SetBranchStatus("AddZ",true  );

  if(isSim){
    RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("ReactionConditions",true );
    RootInput:: getInstance()->GetChain()->SetBranchStatus("fRC_*",true );
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue(){
  Ex.clear();
  Ecm.clear();
  AddBack_EDC.clear();
  AddBack_EDC2.clear();
  AddBack_EDC2.push_back(-1.0); //offset by 1
  EAgata = -1000;
  ELab.clear();
  RawEnergy.clear(); 
  ThetaLab.clear();
  PhiLab.clear();
  ThetaCM.clear();
  //relative_angle = -1000;
  //sum_angle = -1000;
  X.clear();
  Y.clear();
  Z.clear();
  dE= -1000;
  ParticleMult = 0;
  GammaMult = 0;
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  BeamDirection = TVector3(0,0,1);
  elab_tmp = 0;
  thetalab_tmp = 0;
  philab_tmp = 0;

  MG_T=-1000;
  MG_E=-1000;
  MG_X=-1000;
  MG_Y=-1000;
  MG_D=-1000;

  //EDC= -1000;
  EDC.clear();
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct(){
  return (NPL::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_analysis{
    public:
      proxy_analysis(){
        NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
      }
  };

  proxy_analysis p_analysis;
}
