#include <TCanvas.h>
void BeamID(){
    gStyle->SetOptStat("nmeruoi");
  // Load the Main Tree
  TFile* file = new TFile("./Data/Analysis/e748_Physics_12Be.root");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");
  //tree->Print();

  TFile* fileR = new TFile("./Data/Analysis/e748_12Be.root");
  TTree* treeR = (TTree*) fileR->FindObjectAny("ResultTree");
  tree->AddFriend(treeR);
  //treeR->Print();

  // Load the IC chain
  TChain* IC = new TChain("Numexo2");
  IC->Add("./Data/Merged/run_0315.root");
  IC->Add("./Data/Merged/run_0316.root");
  IC->Add("./Data/Merged/run_0317.root");
  IC->Add("./Data/Merged/run_0318.root");
  IC->Add("./Data/Merged/run_0320.root");
  IC->Add("./Data/Merged/run_0321.root");
  IC->Add("./Data/Merged/run_0323.root");
  IC->Add("./Data/Merged/run_0325.root");
  IC->Add("./Data/Merged/run_0326.root");
  IC->Add("./Data/Merged/run_0327.root");
  IC->Add("./Data/Merged/run_0328.root");
  IC->Add("./Data/Merged/run_0329.root");
  IC->Add("./Data/Merged/run_0330.root");
  IC->Add("./Data/Merged/run_0331.root");
  IC->Add("./Data/Merged/run_0339.root");
  IC->Add("./Data/Merged/run_0341.root");
  IC->Add("./Data/Merged/run_0342.root");
  IC->Add("./Data/Merged/run_0346.root");
  IC->Add("./Data/Merged/run_0347.root");
  IC->Add("./Data/Merged/run_0348.root");
  // Friend the two trees
  tree->AddFriend(IC);


  TCanvas* c = new TCanvas("cc","xx",1000,500);
  c->Divide(2,1);
  c->cd(1);
  tree->Draw("IC_E:QPlast>>hIC(4000,0,16000,900,0,900)","TelescopeNumber@.size()==0","colz");
  c->cd(2);
  tree->Draw("QPlast:T_CATS1_CAV>>hIC2(400,0,400,4000,0,16000)","TelescopeNumber@.size()==0","colz");

  auto* c1 {new TCanvas("c1")};
  tree->Draw("T_CATS1_CAV");

  
  }
