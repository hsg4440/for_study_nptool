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

  simulation = true;

  // initialize input and output branches
  if (simulation) {
    Initial = new TInitialConditions();
    ReactionConditions = new TReactionConditions();
  }
  // ExogamData = new TExogamData();
  InitOutputBranch();
  InitInputBranch();
  // get MUST2 and Gaspard objects
  M2 = (TMust2Physics*)m_DetectorManager->GetDetector("M2Telescope");
  // if (!simulation)
  // CATS = (TCATSPhysics*)m_DetectorManager->GetDetector("CATSDetector");

  // get reaction information
  reaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = reaction.GetBeamEnergy();
  std::cout << OriginalBeamEnergy << std::endl;
  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();

  // energy losses
  string light = NPL::ChangeNameToG4Standard(reaction.GetNucleus3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(reaction.GetNucleus1()->GetName());
  cout << light << " " << beam << " " << TargetMaterial << " " << TargetThickness << endl;
  LightTarget = NPL::EnergyLoss(light + "_" + TargetMaterial + ".G4table", "G4Table", 100);
  LightMyl = NPL::EnergyLoss(light + "_Mylar.G4table", "G4Table", 100);
  LightAl = NPL::EnergyLoss(light + "_Al.G4table", "G4Table", 100);
  LightSi = NPL::EnergyLoss(light + "_Si.G4table", "G4Table", 100);
  BeamTarget = NPL::EnergyLoss(beam + "_" + TargetMaterial + ".G4table", "G4Table", 100);

  FinalBeamEnergy = BeamTarget.Slow(OriginalBeamEnergy, TargetThickness * 0.5, 0);
  // FinalBeamEnergy = OriginalBeamEnergy;
  cout << "Original beam energy: " << OriginalBeamEnergy << " MeV      Mid-target beam energy: " << FinalBeamEnergy
       << "MeV " << endl;
  reaction.SetBeamEnergy(FinalBeamEnergy);

  // initialize various parameters
  Rand = TRandom3();
  DetectorNumber = 0;
  ThetaNormalTarget = 0;
  ThetaM2Surface = 0;
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

  // For simulation
  XTarget = 0;
  YTarget = 0;
  BeamDirection = TVector3(0, 0, 1);

  BeamImpact = TVector3(XTarget, YTarget, 0);
  // determine beam energy for a randomized interaction point in target
  // 1% FWHM randominastion (E/100)/2.35
  // reaction.SetBeamEnergy(Rand.Gaus(BeamEnergy, BeamEnergy * 1. / 100. / 2.35));
  // reaction.SetBeamEnergy(BeamEnergy);

  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////////// LOOP on MUST2  ////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  int M2_size = M2->Si_E.size();
  for (unsigned int countMust2 = 0; countMust2 < M2_size; countMust2++) {
    /************************************************/
    // Part 0 : Get the usefull Data
    //  MUST2
    int TelescopeNumber = M2->TelescopeNumber[countMust2];

    /************************************************/
    // Part 1 : Impact Angle
    ThetaM2Surface = 0;
    ThetaNormalTarget = 0;
    TVector3 HitDirection = M2->GetPositionOfInteraction(countMust2) - BeamImpact;
    double Theta = HitDirection.Angle(BeamDirection);

    X = M2->GetPositionOfInteraction(countMust2).X();
    Y = M2->GetPositionOfInteraction(countMust2).Y();
    Z = M2->GetPositionOfInteraction(countMust2).Z();

    ThetaM2Surface = HitDirection.Angle(-M2->GetTelescopeNormal(countMust2));
    ThetaNormalTarget = HitDirection.Angle(TVector3(0, 0, 1));

    /************************************************/

    /************************************************/
    // Part 2 : Impact Energy
    Energy = 0;
    Si_E_M2 = M2->Si_E[countMust2];
    CsI_E_M2 = M2->CsI_E[countMust2];

    // if CsI
    if (CsI_E_M2 > 0) {
      // The energy in CsI is calculate form dE/dx Table because
      Energy = CsI_E_M2;
      Energy = LightMyl.EvaluateInitialEnergy(Energy, 3 * micrometer, ThetaM2Surface);
      Energy = LightAl.EvaluateInitialEnergy(Energy, 0.4 * micrometer, ThetaM2Surface);
      Energy += Si_E_M2;
    }
    else {
      Energy = Si_E_M2;
    }

    Energy = LightAl.EvaluateInitialEnergy(Energy, 0.4 * micrometer, ThetaM2Surface);
    // Evaluate energy using the thickness
    // Target Correction
    Energy = LightTarget.EvaluateInitialEnergy(Energy, TargetThickness * 0.5, Theta);

    // What is written in the tree
    ThetaLab.push_back(Theta / deg);
    Ex.push_back(reaction.ReconstructRelativistic(Energy, Theta));
    ELab.push_back(Energy);
    /************************************************/

  } // end loop MUST2
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  // RootOutput::getInstance()->GetTree()->Branch("Exogam", &ExogamData);
  RootOutput::getInstance()->GetTree()->Branch("Ex", &Ex);
  RootOutput::getInstance()->GetTree()->Branch("ELab", &ELab);
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab", &ThetaLab);
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM", &ThetaCM, "ThetaCM/D");
  RootOutput::getInstance()->GetTree()->Branch("Run", &Run, "Run/I");
  RootOutput::getInstance()->GetTree()->Branch("X", &X, "X/D");
  RootOutput::getInstance()->GetTree()->Branch("Y", &Y, "Y/D");
  RootOutput::getInstance()->GetTree()->Branch("Z", &Z, "Z/D");
  RootOutput::getInstance()->GetTree()->Branch("dE", &dE, "dE/D");
  if (!simulation) {
  }
  else {
    RootOutput::getInstance()->GetTree()->Branch("OriginalELab", &OriginalELab, "OriginalELab/D");
    RootOutput::getInstance()->GetTree()->Branch("OriginalThetaLab", &OriginalThetaLab, "OriginalThetaLab/D");
    RootOutput::getInstance()->GetTree()->Branch("BeamEnergy", &BeamEnergy, "BeamEnergy/D");
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch() {
  // RootInput:: getInstance()->GetChain()->SetBranchAddress("GATCONF",&vGATCONF);
  if (!simulation) {
  }
  else {
    RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions", true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("fIC_*", true);
    RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions", &Initial);
    RootInput::getInstance()->GetChain()->SetBranchStatus("ReactionConditions", true);
    RootInput::getInstance()->GetChain()->SetBranchStatus("fRC_*", true);
    RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions", &ReactionConditions);
  }
  // RootInput::getInstance()->GetChain()->SetBranchStatus("Exogam", true);
  RootInput::getInstance()->GetChain()->SetBranchStatus("fExo*", true);
  // RootInput::getInstance()->GetChain()->SetBranchAddress("Exogam", &ExogamData);
  // RootInput::getInstance()->GetChain()->SetBranchStatus("ZDD", true);
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue() {
  ExNoBeam = ExNoProton = -1000;
  EDC = -1000;
  BeamEnergy = -1000;
  ThetaCM = -1000;
  X = -1000;
  Y = -1000;
  Z = -1000;
  dE = -1000;
  ELab.clear();
  ThetaLab.clear();
  Ex.clear();
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
