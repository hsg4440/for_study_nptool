/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto, Omar Nasr                                  *
 * contact address: pilottoelia@gmail.com, omar.nasr@etu.unicaen.fr          *
 *                                                                           *
 * Creation Date  : September 2021                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 * Use to kill the beam track and replace it with the reaction product       *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

//C++ libraries
#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>//FIXME
#include <Randomize.hh>//FIXME
//G4 libraries
#include "G4Electron.hh"//FIXME
#include "G4Gamma.hh"//FIXME
#include "G4EmCalculator.hh"//FIXME
#include "G4VPhysicalVolume.hh"//FIXME
#include "G4IonTable.hh"//FIXME
#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
//nptool libraries
#include "NPFunction.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSFunction.hh"
//other
#include "SamuraiFieldPropagation.hh"
#include "SamuraiFieldMap.h"
#include "RootOutput.h"
#include "TLorentzVector.h"

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::SamuraiFieldPropagation(G4String modelName, G4Region* envelope)
  : G4VFastSimulationModel(modelName, envelope) {
  //ReadConfiguration();

    m_Map = NULL;
    m_Initialized = false;

    //ABLA = new G4AblaInterface();
  }

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::SamuraiFieldPropagation(G4String modelName)
  : G4VFastSimulationModel(modelName) {}

////////////////////////////////////////////////////////////////////////////////
NPS::SamuraiFieldPropagation::~SamuraiFieldPropagation() {}

////////////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::AttachReactionConditions() {
  // Reasssigned the branch address
  /*
  if (RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress(
		  "ReactionConditions", &m_ReactionConditions);*/
}

////////////////////////////////////////////////////////////////////////////////
/*
void NPS::SamuraiFieldPropagation::ReadConfiguration() {
  NPL::InputParser input(NPOptionManager::getInstance()->GetReactionFile());
  
  if(input.GetAllBlocksWithToken("TwoBodyReaction").size()>0) m_ReactionType =TwoBody;
  else if(input.GetAllBlocksWithToken("QFSReaction").size()>0) m_ReactionType =QFS;
  else if(input.GetAllBlocksWithToken("FusionReaction").size()>0) m_ReactionType =Fusion;

  // Two body
  if (m_ReactionType==TwoBody ) {
    m_Reaction.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_Reaction.GetParticle1()->GetName());
    if(m_Reaction.GetParticle3()->GetName() != ""){
      m_active = true;
      m_ReactionConditions = new TReactionConditions();
      AttachReactionConditions();
      if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
        RootOutput::getInstance()->GetTree()->Branch(
            "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
    }
  } 

  // QFS
  else  if (m_ReactionType==QFS) {
    m_QFS.ReadConfigurationFile(input);
    m_BeamName = NPL::ChangeNameToG4Standard(m_QFS.GetParticleA()->GetName());
    m_active = true;
    m_ReactionConditions = new TReactionConditions();
    AttachReactionConditions();
    if (!RootOutput::getInstance()->GetTree()->FindBranch("ReactionConditions"))
      RootOutput::getInstance()->GetTree()->Branch(
          "ReactionConditions", "TReactionConditions", &m_ReactionConditions);
  }

  // Fusion
  else  if (m_ReactionType==Fusion) {
    vector<InputBlock*> blocks=input.GetAllBlocksWithToken("FusionReaction");
    m_BeamName = NPL::ChangeNameToG4Standard(blocks[0]->GetString("Beam"));
    m_TargetNuclei = blocks[0]->GetString("Target");
    m_FusionProduct = blocks[0]->GetString("Product");
    m_FusionExcitation = blocks[0]->GetDouble("ExcitationEnergy","MeV");
    m_active = true;
    // not used
    m_ReactionConditions = new TReactionConditions();
  }
  else {
    m_active = false;
  }
  
}
*/
////////////////////////////////////////////////////////////////////////////////
G4bool NPS::SamuraiFieldPropagation::IsApplicable(const G4ParticleDefinition& particleType) {
  if (particleType.GetPDGCharge() == 0) return false;
  return true;
}

////////////////////////////////////////////////////////////////////////////////
G4bool NPS::SamuraiFieldPropagation::ModelTrigger(const G4FastTrack& fastTrack) {
  return true;
}

////////////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::DoIt(const G4FastTrack& fastTrack,
    G4FastStep&        fastStep) {
  //std::cout << "\nDOIT" << std::endl;
  //std::cout << "FIELDMAP " << m_FieldMap << std::endl;
  
  if(!m_Initialized){
    m_Map = new SamuraiFieldMap();
    //cout << "before load map\n";
    m_Map->LoadMap(0, m_FieldMap, 10);//FIXME

    //Needed for the RungeKutta method
    m_Map->SetStepLimit( 10.*m / m_StepSize );
    m_Map->SetRmax( m_Rmax );

    m_Initialized = true;
  }
  
  switch (m_Meth){
    case RungeKutta:
      RungeKuttaPropagation(fastTrack, fastStep);
      break;
    case EliaOmar:
      EliaOmarPropagation(fastTrack, fastStep);
      break;
    default:
      cout << endl;
      cout << "//////WARNING///////" << endl;
      cout << "In SamuraiFieldPropagation.cc:\n";
      cout << "Propagation method not defined - THIS MESSAGE SHOULD NEVER APPEAR\n";
      cout << endl;

  }
  
  return;
}

/////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::EliaOmarPropagation (const G4FastTrack& fastTrack,
							G4FastStep& fastStep){
  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  
  static G4VSolid* solid;//Stores the current solid 
  solid = PrimaryTrack->GetTouchable()->GetSolid();
  
  double speed = PrimaryTrack->GetVelocity();
  
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();
  G4ThreeVector localVel = localDir * speed;
  
  G4ThreeVector newPosition;
  G4ThreeVector newDir;
  G4ThreeVector newMomentum;  
  
  //static int count = 0;
  // count++;

  G4ThreeVector B = VtoG4( m_Map->InterpolateB(G4toV(localPosition/mm)) );
  double magB = B.mag() / tesla;

    
  if(magB != 0){

    //unit conversion factors 
    double ConF_p = 1 / (joule * c_light / (m/s)) ; // MeV/c to kg*m/s (SI units)
    double ConF_m = 1 / ( joule * (c_light/(m/s)) * (c_light/(m/s)) ) ; // MeV/c^2 to kg (SI units)
  
    double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
    double mass = PrimaryTrack->GetParticleDefinition()->GetPDGMass() * ConF_m;

    //Calculate curvature radius and pulsation
    G4ThreeVector rho = -(localMomentum * ConF_p).cross(B / tesla);
    rho = rho / (charge * magB * magB) * meter;
    G4ThreeVector omega = ( -(charge/mass) * B / tesla ) * hertz;

    //Distance the particle will travel along B and perpendicular to B
    double L_B = B.unit().dot(localDir) * m_StepSize;
    double L_perp = sqrt(m_StepSize*m_StepSize - L_B*L_B);
    
    //Calculate new position
    double angle = L_perp / rho.mag() * rad;
    G4ThreeVector rho2 = rho;
    rho2.rotate(angle, omega);

    //motion along the B direction
    G4ThreeVector B_motion = L_B * B.unit();

    //Calculating new kinematic properties
    newPosition = localPosition - rho + rho2 + B_motion;
    newDir = localDir;
    newDir.rotate(angle, omega);
    newMomentum = localMomentum.mag() * newDir;

  }
  else{
    
    //Calculating new kinematic properties with no magnetic field
    newPosition = localPosition + localDir * m_StepSize;
    newDir = localDir;
    newMomentum = localMomentum;
  }
  
  if (solid->Inside(newPosition) != kInside){
      G4ThreeVector toOut = solid->DistanceToOut(localPosition, localDir) * localDir;
      newPosition = localPosition + toOut;
  }
  
  double time  = PrimaryTrack->GetGlobalTime()+(newPosition - localPosition).mag()/speed;

  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  return;
  
}

/////////////////////////////////////////////////////////////////////////
void NPS::SamuraiFieldPropagation::RungeKuttaPropagation (const G4FastTrack& fastTrack,
							  G4FastStep& fastStep){
  static int counter = 0;//debugging purposes
  static bool inside = false;//true if previous step is inside the volume
  static vector<TVector3> trajectory;
  static int count;//keeps track of the step reached inside trajectory

  // Get the track info
  const G4Track* PrimaryTrack = fastTrack.GetPrimaryTrack();
  
  G4ThreeVector localDir = fastTrack.GetPrimaryTrackLocalDirection();
  G4ThreeVector localPosition = fastTrack.GetPrimaryTrackLocalPosition();
  G4ThreeVector localMomentum = fastTrack.GetPrimaryTrackLocalMomentum();

  double speed = PrimaryTrack->GetVelocity();
  double charge = PrimaryTrack->GetParticleDefinition()->GetPDGCharge() / coulomb;
  double ConF_p = 1 / (joule * c_light / (m/s)) ; // MeV/c to kg*m/s (SI units)

  //Initially inside is false
  if (!inside){
    count = 2;//skip first two positions as they are the same as the current position
    double Brho = localMomentum.mag() * ConF_p / charge;
    TVector3 pos (localPosition.x(), localPosition.y(), localPosition.z());
    TVector3 dir (localDir.x(), localDir.y(), localDir.z());

    m_Map->SetTimeIntervalSize( m_StepSize / speed );

    trajectory.clear();
    trajectory = m_Map->Propagate(Brho, pos, dir);

    inside = true;
  }

  G4ThreeVector newPosition (trajectory[count].x(), trajectory[count].y(), trajectory[count].z());

  //Check if newPosition is outside
  if (inside){
    G4VSolid* solid = PrimaryTrack->GetTouchable()->GetSolid();

    if (solid->Inside(newPosition) != kInside){
      inside = false;
      G4ThreeVector toOut = solid->DistanceToOut(localPosition, localDir) * localDir;
      newPosition = localPosition + toOut;
      counter++;
    }
  }

  G4ThreeVector newDir = (newPosition - localPosition).unit();
  G4ThreeVector newMomentum = newDir * localMomentum.mag(); 
  double time  = PrimaryTrack->GetGlobalTime()+(newPosition - localPosition).mag()/speed;

  fastStep.ProposePrimaryTrackFinalPosition( newPosition );
  fastStep.SetPrimaryTrackFinalMomentum ( newMomentum );//FIXME
  fastStep.ProposePrimaryTrackFinalTime( time );

  count++;
  
  return;
}


