/*****************************************************************************
 * Copyright (C) 2009-2014    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : march 2025                                               *
 * Last update    :                                                          *
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
using namespace std;
#include"Analysis.h"
#include"NPAnalysisFactory.h"
#include"NPDetectorManager.h"
#include"NPOptionManager.h"
#include"RootOutput.h"
#include"RootInput.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis(){
}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init(){
  M2= (TMust2Physics*) m_DetectorManager->GetDetector("MUST2Array");
  SSSD= (TSSSDPhysics*) m_DetectorManager->GetDetector("SSSD");
  InitOutputBranch();
  InitInputBranch();
  Rand = TRandom3();
  He10Reaction= new NPL::Reaction();
  He10Reaction->ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  DetectorNumber = 0 ;
  ThetaNormalTarget = 0 ;
  ThetaM2Surface = 0; 
  X_M2 = 0 ;
  Y_M2 = 0 ;
  Z_M2 = 0 ;
  Si_E_M2 = 0 ;
  CsI_E_M2 = 0 ; 
  E_SSSD = 0 ;
  Energy = 0;
  E_M2 = 0;
  Si_X_M2 = 0;
  Si_Y_M2 = 0;
  ZTarget = 0;
  TargetThickness = m_DetectorManager->GetTargetThickness();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent(){
    // Reinitiate calculated variable
    ReInitValue();

    // Get the Init information on beam position and energy
    // and apply by hand the experimental resolution
    // This is because the beam diagnosis are not simulated

    // PPAC position resolution on target is assumed to be 1mm
    double XTarget = Rand.Gaus(Initial->GetIncidentPositionX(),1);
    double YTarget = Rand.Gaus(Initial->GetIncidentPositionY(),1);
    TVector3 BeamDirection = Initial->GetBeamDirection();

    // Beam energy is measured using F3 and F2 plastic TOF
    double BeamEnergy = Rand.Gaus(Initial->GetIncidentInitialKineticEnergy(),4.5);
    BeamEnergy = Li11CD2.Slow(BeamEnergy,TargetThickness/2.,0);
    He10Reaction->SetBeamEnergy(BeamEnergy);

    //////////////////////////// LOOP on MUST2 + SSSD Hit //////////////////
    for(unsigned int countSSSD = 0 ; countSSSD < SSSD->Energy.size() ; countSSSD++){
      for(unsigned int countMust2 = 0 ; countMust2 < M2->Si_E.size() ; countMust2++){
        /************************************************/
        //Part 0 : Get the usefull Data
        // MUST2
        int X = M2->Si_X[countMust2];
        int Y = M2->Si_Y[countMust2];
        int TelescopeNumber = M2->TelescopeNumber[countMust2];
        Si_X_M2 = X ;
        Si_Y_M2 = Y ;
        //SSSD
        int SiNumber = SSSD->DetectorNumber[countSSSD];

        /************************************************/
        // Matching between Thin Si and MUST2, and Forward Telescope Only
        if(TelescopeNumber==SiNumber && TelescopeNumber<5){
          DetectorNumber = TelescopeNumber ;
          /************************************************/
          // Part 1 : Impact Angle
          ThetaM2Surface = 0;
          ThetaNormalTarget = 0;
          if(XTarget>-1000 && YTarget>-1000){
            TVector3 BeamImpact(XTarget,YTarget,0);
            TVector3 HitDirection = M2 -> GetPositionOfInteraction(countMust2) - BeamImpact ;
            ThetaLab = HitDirection.Angle( BeamDirection );

            ThetaM2Surface = HitDirection.Angle(- M2 -> GetTelescopeNormal(countMust2) );
            ThetaNormalTarget = HitDirection.Angle( TVector3(0,0,1) ) ;
            X_M2 = M2 -> GetPositionOfInteraction(countMust2).X() ;
            Y_M2 = M2 -> GetPositionOfInteraction(countMust2).Y() ;
            Z_M2 = M2 -> GetPositionOfInteraction(countMust2).Z() ;
          }

          else{
            BeamDirection = TVector3(-1000,-1000,-1000);
            ThetaM2Surface    = -1000  ;
            ThetaNormalTarget = -1000  ;
          }

          /************************************************/

          /************************************************/

          // Part 2 : Impact Energy
          Energy = ELab = 0;
          Si_E_M2 = M2->Si_E[countMust2];
          CsI_E_M2= M2->CsI_E[countMust2];
          E_SSSD = SSSD->Energy[countSSSD];

          // if CsI
          if(CsI_E_M2>0 ){
            // The energy in CsI is calculate form dE/dx Table because 
            // 20um resolution is poor
            Energy = 
              He3Si.EvaluateEnergyFromDeltaE(Si_E_M2,300*micrometer,
                  ThetaM2Surface, 0.01*MeV, 
                  450.*MeV,0.001*MeV ,1000);
            E_M2=CsI_E_M2;
          }

          else
            Energy = Si_E_M2;

          E_M2 += Si_E_M2;

          // Evaluate energy using the thickness 
          ELab = He3Al.EvaluateInitialEnergy( Energy ,2*0.4*micrometer , ThetaM2Surface); 
          ELab = He3Si.EvaluateInitialEnergy( ELab ,20*micrometer , ThetaM2Surface);
          ELab = He3Al.EvaluateInitialEnergy( ELab ,0.4*micrometer , ThetaM2Surface);
          // Target Correction
          ELab   = He3CD2.EvaluateInitialEnergy( ELab ,TargetThickness/2., ThetaNormalTarget);
          /************************************************/

          /***********888888888888888888888****************/
          // Part 3 : Excitation Energy Calculation
          Ex = He10Reaction -> ReconstructRelativistic( ELab , ThetaLab );
          ThetaLab=ThetaLab/deg;
          /************************************************/

          /************************************************/
          // Part 4 : Theta CM Calculation
          ThetaCM  = He10Reaction -> EnergyLabToThetaCM( ELab , 0)/deg;
          /************************************************/
        }
      } //end loop SSSD
    }//end loop MUST2
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End(){
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Ex",&Ex,"Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab",&ELab,"ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab",&ThetaLab,"ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM",&ThetaCM,"ThetaCM/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch(){
  RootInput:: getInstance()->GetChain()->SetBranchAddress("InitialConditions",&Initial);
  RootInput:: getInstance()->GetChain()->SetBranchStatus("InitialConditions",true );
  RootInput:: getInstance()->GetChain()->SetBranchStatus("fIC_*",true );
}

////////////////////////////////////////////////////////////////////////////////     
void Analysis::ReInitValue(){
  Ex = -1000 ;
  ELab = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPA::VAnalysis* Analysis::Construct(){
  return (NPA::VAnalysis*) new Analysis();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy{
  public:
    proxy(){
      NPA::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct);
    }
};

proxy p;
}

