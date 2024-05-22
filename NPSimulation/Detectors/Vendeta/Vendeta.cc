/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Vendeta simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "Vendeta.hh"
#include "ProcessScorers.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"
// CAD Mesj
#include "CADMesh.hh"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Vendeta_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.01*MeV;
  const double ResoTime = 0.454*ns ; // reso Cf
  /* const double ResoTime = 0.61*ns ; // reso U8 */
  /* const double ResoTime = 0.*ns ; // reso U8 */
  const double ResoEnergyLG = 0.43*MeV ;
  const double ResoEnergyHG = 0.255*MeV ;
  const double Thickness = 51.*mm ;
  const double Radius = 127./2*mm ;
  
  // Lead shield
  const double Lead_Radius = 9*cm;
  const double Lead_Thickness = 9*mm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Vendeta Specific Method
Vendeta::Vendeta(){
  m_Event = new TVendetaData() ;
  m_VendetaScorer = 0;
  m_VendetaDetector = 0;
  m_SensitiveCell = 0;
  m_MecanicalStructure = 0;  
  m_Build_MecanicalStructure = 0;
  m_LeadShield = 0;
  m_BuildLeadShield = 1;

  // RGB Color + Transparency
  m_VisAl      = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));   
  m_VisEJ309   = new G4VisAttributes(G4Colour(0.2, 0.85, 0.85, 1));   
  m_VisMuMetal = new G4VisAttributes(G4Colour(0.55, 0.5, 0.5, 0.7));   
  m_VisPyrex   = new G4VisAttributes(G4Colour(0.1, 0.5, 0.7, 1));   
  m_VisEJ560   = new G4VisAttributes(G4Colour(0.6, 0.6, 0.2, 1));   
  m_VisInox    = new G4VisAttributes(G4Colour(0.6, 0.5, 0.6, 1));   
  m_VisLeadShield    = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 1));   

  // Material definition
  m_Vacuum  = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_Al      = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_Inox    = MaterialManager::getInstance()->GetMaterialFromLibrary("Inox");
  m_EJ309   = MaterialManager::getInstance()->GetMaterialFromLibrary("EJ309");
  m_EJ560   = MaterialManager::getInstance()->GetMaterialFromLibrary("EJ560");
  m_Pyrex   = MaterialManager::getInstance()->GetMaterialFromLibrary("Pyrex");
  m_MuMetal = MaterialManager::getInstance()->GetMaterialFromLibrary("mumetal");
}

Vendeta::~Vendeta(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Vendeta::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Vendeta::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Vendeta::BuildSensitiveCell(){

  if(!m_SensitiveCell){
    G4Tubs* scin_cell = new G4Tubs("Vendeta_scin", 0, Vendeta_NS::Radius, Vendeta_NS::Thickness*0.5, 0., 360*deg);
    m_SensitiveCell = new G4LogicalVolume(scin_cell,m_EJ309,"logic_Vendeta_scin",0,0,0);
    m_SensitiveCell->SetVisAttributes(m_VisEJ309);
    m_SensitiveCell->SetSensitiveDetector(m_VendetaScorer);
  }

  return m_SensitiveCell;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Vendeta::BuildVendetaDetector(){
  string basepath = getenv("NPTOOL");
  string path = basepath + "/NPSimulation/Detectors/Vendeta/stl/";
  string stl_file;

  G4ThreeVector Pos = G4ThreeVector(0,0,0);
  Pos.setX(-160);
  Pos.setY(117.5);
  G4RotationMatrix* Rot = new G4RotationMatrix();
  Rot->rotateX(90*deg);

  if(!m_VendetaDetector){

    m_VendetaDetector = new G4AssemblyVolume();

    // *** Sensitive Volume *** //
    G4ThreeVector Pos_cell = G4ThreeVector(0,0,0.5*Vendeta_NS::Thickness + 2.*mm);
    m_VendetaDetector->AddPlacedVolume(BuildSensitiveCell(),Pos_cell,0);

    
    G4LogicalVolume* LogicVolume;
    // *** 1 *** //
    stl_file = path + "cell_body.stl";
    auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) stl_file.c_str());
    mesh->SetScale(mm);
    
    auto cad_solid = mesh->GetSolid();
    LogicVolume = new G4LogicalVolume(cad_solid,m_Al,"cell_body",0,0,0);
    LogicVolume->SetVisAttributes(m_VisAl);
    m_VendetaDetector->AddPlacedVolume(LogicVolume,Pos,Rot);

    
    // *** 2 *** //
    stl_file = path + "Vendeta_Al.stl";
    mesh = CADMesh::TessellatedMesh::FromSTL((char*) stl_file.c_str());
    mesh->SetScale(mm);
    
    cad_solid = mesh->GetSolid();
    LogicVolume = new G4LogicalVolume(cad_solid,m_Al,"Vendeta_Al",0,0,0);
    LogicVolume->SetVisAttributes(m_VisAl);
    m_VendetaDetector->AddPlacedVolume(LogicVolume,Pos,Rot);

    // *** 3 *** //
    stl_file = path + "mumetal_shield.stl";
    mesh = CADMesh::TessellatedMesh::FromSTL((char*) stl_file.c_str());
    mesh->SetScale(mm);
    
    cad_solid = mesh->GetSolid();
    LogicVolume = new G4LogicalVolume(cad_solid,m_MuMetal,"mumetal",0,0,0);
    LogicVolume->SetVisAttributes(m_VisMuMetal);
    m_VendetaDetector->AddPlacedVolume(LogicVolume,Pos,Rot);
 
    // *** 4 *** //
    stl_file = path + "EJ560.stl";
    mesh = CADMesh::TessellatedMesh::FromSTL((char*) stl_file.c_str());
    mesh->SetScale(mm);
    
    cad_solid = mesh->GetSolid();
    LogicVolume = new G4LogicalVolume(cad_solid,m_EJ560,"EJ560",0,0,0);
    LogicVolume->SetVisAttributes(m_VisEJ560);
    m_VendetaDetector->AddPlacedVolume(LogicVolume,Pos,Rot);
 
    // *** 5 *** //
    stl_file = path + "quartz_window.stl";
    mesh = CADMesh::TessellatedMesh::FromSTL((char*) stl_file.c_str());
    mesh->SetScale(mm);
    
    cad_solid = mesh->GetSolid();
    LogicVolume = new G4LogicalVolume(cad_solid,m_Pyrex,"quartz_window",0,0,0);
    LogicVolume->SetVisAttributes(m_VisPyrex);
    m_VendetaDetector->AddPlacedVolume(LogicVolume,Pos,Rot);
 

    if(m_BuildLeadShield){
      G4Tubs* lead = new G4Tubs("lead_shield", 0, Vendeta_NS::Lead_Radius, Vendeta_NS::Lead_Thickness*0.5, 0, 360*deg);
      G4Material* LeadMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
      m_LeadShield = new G4LogicalVolume(lead, LeadMaterial, "logic_lead_shield",0,0,0);
      m_LeadShield->SetVisAttributes(m_VisLeadShield);
      Pos.setX(0);
      Pos.setY(0);
      Pos.setZ( Vendeta_NS::Lead_Thickness*0.5-5*mm);
      m_VendetaDetector->AddPlacedVolume(m_LeadShield, Pos, 0);
    }

  }
  return m_VendetaDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Vendeta::BuildMecanicalStructure()
{
  if(!m_MecanicalStructure){
    m_MecanicalStructure = new G4AssemblyVolume();

    G4RotationMatrix* Rot = new G4RotationMatrix();
    G4ThreeVector Pos = G4ThreeVector(0,0,0);

    string basepath = getenv("NPTOOL");
    // *** Steel part of the strucutre *** //
    string path = basepath + "/NPSimulation/Detectors/Vendeta/Structure_meca_stl/Structure_Acier.stl";

    auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
    mesh->SetScale(mm);

    auto cad_solid = mesh->GetSolid();
    m_MecanicalStructure_Steel = new G4LogicalVolume(cad_solid,m_Inox,"Structure_Steel",0,0,0);
    m_MecanicalStructure_Steel->SetVisAttributes(m_VisInox);

    m_MecanicalStructure->AddPlacedVolume(m_MecanicalStructure_Steel,Pos,Rot);


    // *** Aluminium part *** //
    path = basepath + "/NPSimulation/Detectors/Vendeta/Structure_meca_stl/Structure_Alu.stl";
    mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
    mesh->SetScale(mm);

    cad_solid = mesh->GetSolid();
    m_MecanicalStructure_Al = new G4LogicalVolume(cad_solid,m_Al,"Structure_Al",0,0,0);
    m_MecanicalStructure_Al->SetVisAttributes(m_VisAl);

    m_MecanicalStructure->AddPlacedVolume(m_MecanicalStructure_Al,Pos,Rot);

  }

  return m_MecanicalStructure;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Vendeta::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Vendeta");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Vendeta " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Vendeta " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
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
void Vendeta::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    //Det_pos+=Det_pos.unit()*Vendeta_NS::Thickness*0.5;
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

    BuildVendetaDetector()->MakeImprint(world,Det_pos,Rot,i+1);
  }

  if(m_Build_MecanicalStructure==1){
    G4RotationMatrix* RotMeca = new G4RotationMatrix();
    G4ThreeVector PosMeca = G4ThreeVector(0,0,0);

    BuildMecanicalStructure()->MakeImprint(world,PosMeca,RotMeca);
  }

  G4double platformWidth = 3.0 * m;
  G4double platformLength = 3.0 * m;
  G4double platformThickness = 1.2 * cm;
  G4ThreeVector *platformPosition = new G4ThreeVector(0,0,-1000);
  G4Box* platformSolid = new G4Box("PlatformSolid", platformWidth / 2.0, platformThickness / 2.0, platformLength / 2.0);
  G4LogicalVolume* platformLogical = new G4LogicalVolume(platformSolid, m_Al, "PlatformLogical");
  G4PVPlacement* platformPhysical = new G4PVPlacement(0, G4ThreeVector(0,-1000,0), platformLogical, "PlatformPhysical",world,0,false);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Vendeta::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Vendeta")){
    pTree->Branch("Vendeta", "TVendetaData", &m_Event) ;
  }
  pTree->SetBranchAddress("Vendeta", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Vendeta::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_VendetaScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    /* double Energy =Scorer->GetEnergy(i); */
    double EnergyLG = RandGauss::shoot(Scorer->GetEnergy(i),Vendeta_NS::ResoEnergyLG);
    double EnergyHG = RandGauss::shoot(Scorer->GetEnergy(i),Vendeta_NS::ResoEnergyHG);
    //Convert enrgy deposit to Light Output thanks the F. Pino et al. formula 
    double LightOutHG = 0.62*EnergyHG-1.3*(1-exp(-0.39*pow(EnergyHG,0.97))); // F. Pino et al
    double LightOutLG = 0.62*EnergyLG-1.3*(1-exp(-0.39*pow(EnergyLG,0.97))); // F. Pino et al
    // Apply Resolution on Charge measured on Vendeta data
    /* LightOut = RandGauss::shoot(LightOut, Vendeta_NS::ResoEnergy); */
    if( EnergyHG > Vendeta_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Vendeta_NS::ResoTime);
      int DetectorNbr = level[0]-1;

      // Filling HG
      m_Event->SetHGDetectorNbr(DetectorNbr);
      m_Event->SetHGQ1(LightOutHG*50000);
      m_Event->SetHGQ2(LightOutHG*50000);
      m_Event->SetHGTime(Time); 
      m_Event->SetHGQmax(0); 

      // Filling LG
      m_Event->SetLGDetectorNbr(DetectorNbr);
      m_Event->SetLGQ1(LightOutLG*20000);
      m_Event->SetLGQ2(LightOutLG*20000);
      m_Event->SetLGTime(Time); 
      m_Event->SetLGQmax(0); 

    }
  }

  ///////////
  // Process scorer
  ProcessScorers::PS_Process* Process_scorer = (ProcessScorers::PS_Process*) m_VendetaScorer->GetPrimitive(2);
  unsigned int ProcessMult = Process_scorer->GetMult();
  if(ProcessMult>0){
    string particle_name = Process_scorer->GetParticleName(0);
    if(particle_name=="gamma"){
      m_Event->SetHGIsSat(1); 
      m_Event->SetLGIsSat(1); 
    }

    if(particle_name=="neutron"){
      m_Event->SetHGIsSat(0); 
      m_Event->SetLGIsSat(0); 
    }

  }
  //for(unsigned int i=0; i<ProcessMult; i++){
  //  string particle_name = Process_scorer->GetParticleName(i);
  //}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Vendeta::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_VendetaScorer = CheckScorer("VendetaScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  G4VPrimitiveScorer* Process= new ProcessScorers::PS_Process("Process", 0) ;
  //and register it to the multifunctionnal detector
  m_VendetaScorer->RegisterPrimitive(Calorimeter);
  m_VendetaScorer->RegisterPrimitive(Interaction);
  m_VendetaScorer->RegisterPrimitive(Process);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_VendetaScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Vendeta::Construct(){
  return  (NPS::VDetector*) new Vendeta();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Vendeta{
    public:
      proxy_nps_Vendeta(){
        NPS::DetectorFactory::getInstance()->AddToken("Vendeta","Vendeta");
        NPS::DetectorFactory::getInstance()->AddDetector("Vendeta",Vendeta::Construct);
      }
  };

  proxy_nps_Vendeta p_nps_Vendeta;
}
