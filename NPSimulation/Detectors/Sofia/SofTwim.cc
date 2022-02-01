/******************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project          *
 *                                                                            *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                 *
 * For the list of contributors see $NPTOOL/Licence/Contributors              *
 ******************************************************************************/

/******************************************************************************
 * Original Author: Pierre Morfouace contact address: pierre.morfouace2@cea.fr*
 *                                                                            *
 * Creation Date  : November 2020                                             *
 * Last update    :                                                           *
 *----------------------------------------------------------------------------*
 * Decription:                                                                *
 *  This class describe a simple SofTwim setup for simulation                   *
 *                                                                            *
 *----------------------------------------------------------------------------*
 * Comment:                                                                   *
 *                                                                            *
 ******************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

// NPTool header
#include "SofTwim.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SofTwim_NS{
  // Energy and time Resolution
  const double EnergyThreshold       = 0.1*MeV;
  const double ResoTime              = 0.007*ns;
  const double ResoEnergy            = 1.0*MeV;

  const double twim_anode_width = 10.*cm;
  const double twim_anode_height = 10.*cm;
  const double twim_anode_thickness = 3.1*cm;

  const double twim_cathode_width = 30*um;
  const double twim_cathode_height = 20*cm;
  const double twim_cathode_length = 50*cm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SofTwim Specific Method
SofTwim::SofTwim(){
  m_Event = new TSofTwimData() ;
  m_TwimScorer = 0;

  m_TwinMusic= 0;
  m_TwimGas= "P10_1atm";
  m_Pressure= 1*bar;

  // RGB Color + Transparency
  m_VisCathode = new G4VisAttributes(G4Colour(0.7, 0.4, 0.1, 1));   
  m_VisSquare = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 1));   

}

SofTwim::~SofTwim(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTwim::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTwim::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTwim::BuildTwinMusic(){

  double twim_width_out = 22*cm;
  double twim_height_out = 22*cm;
  double twim_thickness_out = 55*cm;
  double twim_width_in = 21*cm;
  double twim_height_in = 21*cm;
  double twim_thickness_in = 55.1*cm;


  double sector_width = 10*cm;
  double sector_height = 10*cm;
  double sector_thickness = 50*cm;

  if(!m_TwinMusic){
    m_TwinMusic = new G4AssemblyVolume();
    
    // Full Twin volume
    G4Box* Twimbox1 = new G4Box("Twim_Box1", twim_width_in*0.5, twim_height_in*0.5, twim_thickness_in*0.5);
    G4Box* Twimbox2 = new G4Box("Twim_Box2", twim_width_out*0.5, twim_height_out*0.5, twim_thickness_out*0.5);
    G4VSolid* Twimbox = new G4SubtractionSolid("Twim_Box",Twimbox2,Twimbox1);
    G4Material* TwimMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4LogicalVolume* LogicTwimBox = new G4LogicalVolume(Twimbox, TwimMaterial, "logic_twim", 0,0,0);

    //m_VisSquare->SetForceWireframe(1);
    LogicTwimBox->SetVisAttributes(m_VisSquare);

    // Sector Twin volume
    G4Box* Sectorbox = new G4Box("Sector_Box", sector_width*0.5, sector_height*0.5, sector_thickness*0.5);
    G4LogicalVolume* LogicalSector = new G4LogicalVolume(Sectorbox, TwimMaterial, "logic_twim", 0,0,0);
    LogicalSector->SetVisAttributes(G4VisAttributes::GetInvisible());
    // Drift Anode Area //
    G4Box* Anodebox = new G4Box("Anode_Box", SofTwim_NS::twim_anode_width*0.5, SofTwim_NS::twim_anode_height*0.5, SofTwim_NS::twim_anode_thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_TwimGas);
    G4LogicalVolume* m_AnodeDriftArea = new G4LogicalVolume(Anodebox, DetectorMaterial, "logic_twim_anode", 0, 0, 0);
    G4VisAttributes* m_VisTwimAnode = new G4VisAttributes(G4Colour(0.3,0.4,0.5,0.5));
    m_AnodeDriftArea->SetVisAttributes(m_VisTwimAnode);
    m_AnodeDriftArea->SetSensitiveDetector(m_TwimScorer);

    // Cathode plane in the middle //
    G4Box* cathode_box = new G4Box("Cathode_box", SofTwim_NS::twim_cathode_width*0.5, SofTwim_NS::twim_cathode_height*0.5, SofTwim_NS::twim_cathode_length*0.5);
    G4Material* CathodeMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4LogicalVolume* LogicalCathode = new G4LogicalVolume(cathode_box, CathodeMaterial, "logic_cathode", 0,0,0);
    LogicalCathode->SetVisAttributes(m_VisCathode);

    G4RotationMatrix* Rv = new G4RotationMatrix(0,0,0);
    Rv->rotateZ(90*deg);
    G4ThreeVector Tv;
    Tv.setX(0);
    Tv.setY(0);
    Tv.setZ(0);

    m_TwinMusic->AddPlacedVolume(LogicTwimBox,Tv,0);
    m_TwinMusic->AddPlacedVolume(LogicalCathode,Tv,Rv);

    for(unsigned int i=0; i<4; i++){
      if(i==2){
        Tv.setX(0.5*SofTwim_NS::twim_anode_width+SofTwim_NS::twim_cathode_width);
        Tv.setY(0.5*SofTwim_NS::twim_anode_height);
      }
      if(i==3){
        Tv.setX(-0.5*SofTwim_NS::twim_anode_width-SofTwim_NS::twim_cathode_width);
        Tv.setY(0.5*SofTwim_NS::twim_anode_height);
      }
      if(i==0){
        Tv.setX(-0.5*SofTwim_NS::twim_anode_width-SofTwim_NS::twim_cathode_width);
        Tv.setY(-0.5*SofTwim_NS::twim_anode_height);
      }
      if(i==1){
        Tv.setX(0.5*SofTwim_NS::twim_anode_width+SofTwim_NS::twim_cathode_width);
        Tv.setY(-0.5*SofTwim_NS::twim_anode_height);
      }

      Tv.setZ(0);
      int anode_nbr= 0;
      for(unsigned int j=0; j<16; j++){
        anode_nbr++;
        Tv.setZ(j*SofTwim_NS::twim_anode_thickness -7.5*SofTwim_NS::twim_anode_thickness);
        m_TwinMusic->AddPlacedVolume(m_AnodeDriftArea,Tv,0);
      }
    }
  }

  return m_TwinMusic;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void SofTwim::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTwim");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","TwimGas","Pressure"};
  vector<string> sphe = {"R","Theta","Phi","TwimGas","Pressure"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTwim " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_TwimGas = blocks[i]->GetString("TwimGas");
      m_Pressure = blocks[i]->GetDouble("Pressure","bar");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTwim " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_TwimGas = blocks[i]->GetString("TwimGas");
      m_Pressure = blocks[i]->GetDouble("Pressure","bar");

      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void SofTwim::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*SofTwim_NS::twim_cathode_length*0.5;
    // Building Detector reference frame
    G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double kk = -sin(m_Theta[i]);

    G4ThreeVector Y(ii,jj,kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v = v.unit();
    u = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
    BuildTwinMusic()->MakeImprint(world,Det_pos,Rot);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SofTwim::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SofTwim")){
    pTree->Branch("SofTwim", "TSofTwimData", &m_Event) ;
  }
  pTree->SetBranchAddress("SofTwim", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void SofTwim::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // TWIM scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_TwimScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SofTwim_NS::ResoEnergy);
    if(Energy>SofTwim_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),SofTwim_NS::ResoTime);
      int SectionNbr;
      int AnodeNbr = level[1]-2;
      if(AnodeNbr<17){
        SectionNbr = 1;
        AnodeNbr = AnodeNbr;
      }
      else if(AnodeNbr>16 && AnodeNbr<33){
        SectionNbr = 2;
        AnodeNbr = AnodeNbr-(SectionNbr-1)*16;
      }
      else if(AnodeNbr>32 && AnodeNbr<49){  
        SectionNbr = 3;
        AnodeNbr = AnodeNbr-(SectionNbr-1)*16;
      }
      else if(AnodeNbr>48){
        SectionNbr = 4;
        AnodeNbr = AnodeNbr-(SectionNbr-1)*16;
      }
      m_Event->SetSectionNbr(SectionNbr);
      m_Event->SetAnodeNbr(AnodeNbr);
      m_Event->SetEnergy(Energy);
      m_Event->SetDriftTime(Time); 
    }
  }
  Scorer->clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void SofTwim::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_TwimScorer = CheckScorer("TwimScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; 
  level.push_back(1);
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  //G4VPrimitiveScorer* TwinCalorimeter= new CalorimeterScorers::PS_Calorimeter("TwinCalorimeter",level, 0) ;
  //G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_TwimScorer->RegisterPrimitive(Calorimeter);
  //m_TwimScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_TwimScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SofTwim::Construct(){
  return  (NPS::VDetector*) new SofTwim();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_SofTwim{
    public:
      proxy_nps_SofTwim(){
        NPS::DetectorFactory::getInstance()->AddToken("SofTwim","SofTwim");
        NPS::DetectorFactory::getInstance()->AddDetector("SofTwim",SofTwim::Construct);
      }
  };

  proxy_nps_SofTwim p_nps_SofTwim;
}
