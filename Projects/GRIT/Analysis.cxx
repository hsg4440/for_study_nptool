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
#include <iostream>
using namespace std;

#include "Analysis.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPFunction.h"
#include "NPOptionManager.h"
////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis() {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
  if (NPOptionManager::getInstance()->HasDefinition("simulation")) {
    cout << "Considering input data as simulation" << endl;
    simulation = true;
  }
  else {
    cout << "Considering input data as real" << endl;
    simulation = false;
  }
  agata_zShift = 51 * mm;

  // initialize input and output branches
  if (simulation) {
    Initial = new TInitialConditions();
    ReactionConditions = new TReactionConditions();
  }

  InitOutputBranch();
  InitInputBranch();
  // get MUST2 and Gaspard objects
  GPDTrack = (GaspardTracker*)m_DetectorManager->GetDetector("GaspardTracker");

  // get reaction information
  reaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction.GetBeamEnergy();
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  // Cryo target case
  WindowsThickness = 0;        // m_DetectorManager->GetWindowsThickness();
  string WindowsMaterial = ""; // m_DetectorManager->GetWindowsMaterial();

  // energy losses
  string light = NPL::ChangeNameToG4Standard(reaction.GetNucleus3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(reaction.GetNucleus1()->GetName());
  LightTarget = NPL::EnergyLoss(light + "_" + TargetMaterial + ".G4table", "G4Table", 100);
  //  LightAl = NPL::EnergyLoss(light + "_Al.G4table", "G4Table", 100);
  LightSi = NPL::EnergyLoss(light + "_Si.G4table", "G4Table", 100);
  BeamCD2 = NPL::EnergyLoss(beam + "_" + TargetMaterial + ".G4table", "G4Table", 100);

  FinalBeamEnergy = BeamCD2.Slow(OriginalBeamEnergy, TargetThickness * 0.5, 0);
  // FinalBeamEnergy = OriginalBeamEnergy;
  cout << "Original beam energy: " << OriginalBeamEnergy << " MeV      Mid-target beam energy: " << FinalBeamEnergy
       << "MeV " << endl;
  reaction.SetBeamEnergy(FinalBeamEnergy);

  if (WindowsThickness) {
    cout << "Cryogenic target with windows" << endl;
    BeamWindow = new NPL::EnergyLoss(beam + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
    LightWindow = new NPL::EnergyLoss(light + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
  }

  else {
    BeamWindow = NULL;
    LightWindow = NULL;
  }

  // initialize various parameters
  Rand = TRandom3();
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
  ThetaMGSurface = 0;
  Si_E_M2 = 0;
  CsI_E_M2 = 0;
  Energy = 0;
  ThetaGDSurface = 0;
  X = 0;
  Y = 0;
  Z = 0;
  dE = 0;
  BeamDirection = TVector3(0, 0, 1);
  nbTrack = 0;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
  // Reinitiate calculated variable
  ReInitValue();
  double XTarget, YTarget;
  TVector3 BeamDirection;
  XTarget = 0;
  YTarget = 0;
  BeamDirection = TVector3(0, 0, 1);
  OriginalELab = ReactionConditions->GetKineticEnergy(0);
  OriginalThetaLab = ReactionConditions->GetTheta(0);
  BeamEnergy = ReactionConditions->GetBeamEnergy();
  BeamImpact = TVector3(XTarget, YTarget, 0);
  // determine beam energy for a randomized interaction point in target
  // 1% FWHM randominastion (E/100)/2.35
  // reaction.SetBeamEnergy(Rand.Gaus(ReactionConditions->GetIncidentFinalKineticEnergy(),ReactionConditions->GetIncidentFinalKineticEnergy()/235));

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// GASPARD  //////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  if ((ELab = GPDTrack->GetEnergyDeposit()) > 0) {
    TVector3 HitDirection = GPDTrack->GetPositionOfInteraction() - BeamImpact;
    ELab = GPDTrack->GetEnergyDeposit();
    ThetaLab = HitDirection.Angle(BeamDirection);
    ELab = LightTarget.EvaluateInitialEnergy(ELab, TargetThickness * 0.5, 0);

    // Part 3 : Excitation Energy Calculation
    Ex = reaction.ReconstructRelativistic(ELab, ThetaLab);

    // Part 4 : Theta CM Calculation
    ThetaCM = reaction.EnergyLabToThetaCM(ELab, ThetaLab) / deg;
    ThetaLab = ThetaLab / deg;
  }
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("Ex", &Ex, "Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab", &ELab, "ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab", &ThetaLab, "ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM", &ThetaCM, "ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Run", &Run, "Run/I");
  RootOutput::getInstance()->GetTree()->Branch("X", &X, "X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y", &Y, "Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z", &Z, "Z/D");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch() {
  // RootInput:: getInstance()->GetChain()->SetBranchAddress("GATCONF",&vGATCONF);
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions", true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*", true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions", &Initial);
  RootInput::getInstance()->GetChain()->SetBranchStatus("ReactionConditions", true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fRC_*", true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions", &ReactionConditions);
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue() {
  Ex = -1000;
  ELab = -1000;
  BeamEnergy = -1000;
  ThetaLab = -1000;
  ThetaCM = -1000;
  X = -1000;
  Y = -1000;
  Z = -1000;
  dE = -1000;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the AnalysisFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct() { return (NPL::VAnalysis*)new Analysis(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_analysis {
 public:
  proxy_analysis() { NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct); }
};

proxy_analysis p_analysis;
}

