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

#ifndef SamuraiFieldPropagation_h
#define SamuraiFieldPropagation_h

#include "G4VFastSimulationModel.hh"
#include "G4Abla.hh"
#include "G4AblaInterface.hh"
#include "G4Fragment.hh"
#include "PhysicsList.hh"
#include "NPReaction.h"
#include "NPQFS.h"
#include "TReactionConditions.h"

//FIXME
#include <string>
#include <vector>

class G4VPhysicalVolume;
class SamuraiFieldMap;
namespace NPS{
  enum PropagationMethod{
    RungeKutta,
    EliaOmar
    };
  
  class SamuraiFieldPropagation : public G4VFastSimulationModel{
    public:
      SamuraiFieldPropagation (G4String, G4Region*);
      SamuraiFieldPropagation (G4String);
      ~SamuraiFieldPropagation ();

    public:
    //void ReadConfiguration();
      G4bool IsApplicable(const G4ParticleDefinition&);
      G4bool ModelTrigger(const G4FastTrack &);
      void DoIt(const G4FastTrack&, G4FastStep&);
 
    private:
    //G4AblaInterface* ABLA;

    //Utility functions
    vector <double> G4toV(G4ThreeVector V){
      return vector <double> { V.x(),V.y(),V.z() }; 
    }
    G4ThreeVector VtoG4(vector <double> V){
      return G4ThreeVector( V.at(0), V.at(1), V.at(2) ) ;
    }
    
    string Prt (G4ThreeVector V){//FIXME
      return to_string(V.getR()) + " " + to_string(V.getPhi()) + " " + to_string(V.getTheta());
      }

    string Cart (G4ThreeVector V){//FIXME
      return to_string(V.getX()) + " " + to_string(V.getY()) + " " + to_string(V.getZ());
      }

    bool m_Initialized; //Map Loaded
    double m_StepSize;
    double m_Angle; //Angle of rotation of the whole magnet - FIXME
    string m_FieldMap;//FIXME
    PropagationMethod m_Meth;// "I AM THE ONE WHO KNOCKS!!!" cit. Heisenberg
    SamuraiFieldMap* m_Map;

    G4ThreeVector MagField (G4ThreeVector pos);//FIXME

    void EliaOmarPropagation (const G4FastTrack& fastTrack, G4FastStep& fastStep);
    void RungeKuttaPropagation (const G4FastTrack& fastTrack, G4FastStep& fastStep);
    
    
      
    private:
    //TReactionConditions* m_ReactionConditions;
 
    public:
    void AttachReactionConditions();//FIXME
    
    void SetStepSize(double step){m_StepSize=step;};
    void SetAngle(double angle){m_Angle=angle;};//FIXME
    void SetFieldMap(string fieldmap){m_FieldMap=fieldmap;};//FIXME
    void SetMethod(PropagationMethod method){m_Meth=method;};
  };
}


#endif 
