// c++ includes
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

using namespace std;

int NumberOfDetectors= 72;
int NumberOfAnodes= 11;
int nentries=1e6;

/////////////////////////////////////
int main(int argc, char** argv){
  int run = atof(argv[1]);

  auto chain = new TChain("PhysicsTree");
  chain->Add(Form("/home/cyril/Workflow/nptool/Outputs/Analysis/run%i/run%i.root",run,run));

  nentries = chain->GetEntries();
  /* nentries = 50000000; */
  cout << "Number of entries: " << nentries << endl;

  TFile* ofile = new TFile(Form("histos_ToF/histo_tof_file_run%i.root",run),"recreate");
  TH2F* hLG_LSci[11];
  TH2F* hHG_LSci[11];
  TH1F* hLG[791];
  TH1F* hHG[791];
  

  TGraph *gFlyPath = new TGraph();
  gFlyPath->SetTitle(" ; Vendeta-Anode Index ; fly length (mm);");
  
  nentries = chain->GetEntries();
  cout << "Number of entries: " << nentries << endl;

  vector<double>* FC_Q1 = new vector<double>();
  vector<int>* FC_Anode_ID = new vector<int>();
  vector<bool>* FC_FakeFission = new vector<bool>();

  vector<double>* LG_Tof = new vector<double>();
  vector<int>* LG_ID = new vector<int>();
  vector<double>* LG_Q1 = new vector<double>();
  vector<double>* LG_Q2 = new vector<double>();
  vector<double>* LG_FlyPath = new vector<double>();

  vector<double>* HG_Tof = new vector<double>();
  vector<int>* HG_ID = new vector<int>();
  vector<double>* HG_Q1 = new vector<double>();
  vector<double>* HG_Q2 = new vector<double>();
  vector<double>* HG_FlyPath = new vector<double>();

  chain->SetBranchAddress("FC_Q1",&FC_Q1);
  chain->SetBranchAddress("FC_Anode_ID",&FC_Anode_ID);
  chain->SetBranchAddress("FC_FakeFission",&FC_FakeFission);
  chain->SetBranchAddress("LG_Tof",&LG_Tof);
  chain->SetBranchAddress("LG_ID",&LG_ID);
  chain->SetBranchAddress("LG_Q1",&LG_Q1);
  chain->SetBranchAddress("LG_Q2",&LG_Q2);
  chain->SetBranchAddress("LG_FlyPath",&LG_FlyPath);
  chain->SetBranchAddress("HG_Tof",&HG_Tof);
  chain->SetBranchAddress("HG_ID",&HG_ID);
  chain->SetBranchAddress("HG_Q1",&HG_Q1);
  chain->SetBranchAddress("HG_Q2",&HG_Q2);
  chain->SetBranchAddress("HG_FlyPath",&HG_FlyPath);

  for(int j=0; j<NumberOfAnodes; j++){

    TString histo_name = Form("hLG_LSci_Anode%i",j+1);
    hLG_LSci[j] = new TH2F(histo_name,histo_name,72,1,73,2000,-100,300);

    histo_name = Form("hHG_LSci_Anode%i",j+1);
    hHG_LSci[j] = new TH2F(histo_name,histo_name,72,1,73,2000,-100,300);

    for(int i=0; i<NumberOfDetectors; i++){
      int index = j + i*NumberOfAnodes;
      histo_name = Form("hLG_Det%i_Anode%i",i+1,j+1);
      hLG[index] = new TH1F(histo_name,histo_name,2000,-100,300);

      histo_name = Form("hHG_Det%i_Anode%i",i+1,j+1);
      hHG[index] = new TH1F(histo_name,histo_name,2000,-100,300);
    }
  }
  for(int i=0; i<nentries; i++){
    chain->GetEntry(i);

    if(i%100000==0){
      cout << setprecision(2) << "\033[34m\r Processing run "<< run <<" ..." << (double)i/nentries*100 << "\% done" << flush;
      cout << setprecision(2) << "\033[34m\r Processing tree..." << (double)i/nentries*100 << "\% done" << flush;
    }

    if(FC_Anode_ID->size()>0){
      bool Fake = FC_FakeFission->at(0);
      int anode = FC_Anode_ID->at(0);


      int mysize = LG_Tof->size();
      for(int j=0; j<mysize; j++){
        // LG //
        int index_LG = (anode-1) + (LG_ID->at(j)-1)*NumberOfAnodes;
        double PSD = LG_Q2->at(j)/LG_Q1->at(j);
        if(LG_ID->at(j)>0 && anode>0 && Fake==0 && LG_Q1->at(j)>2000){
          hLG[index_LG]->Fill(LG_Tof->at(j));
          hLG_LSci[anode-1]->Fill(LG_ID->at(j),LG_Tof->at(j));
          gFlyPath->SetPoint(index_LG,index_LG,LG_FlyPath->at(j));
          hLG[index_LG]->Fill(LG_Tof->at(j));		
        }
      }

      mysize = HG_Tof->size();
      for(int j=0; j<mysize; j++){
        // HG //
        int index_HG = (anode-1) + (HG_ID->at(j)-1)*NumberOfAnodes;
        double PSD = HG_Q2->at(j)/HG_Q1->at(j);
        if(HG_ID->at(j)>0 && anode>0 && Fake==0 && HG_ID->size()==1){
          hHG[index_HG]->Fill(HG_Tof->at(j));
          hHG_LSci[anode-1]->Fill(HG_ID->at(j),HG_Tof->at(j));
          hHG[index_HG]->Fill(HG_Tof->at(j));			
        }
      }
    }
  }

  ofile->WriteObject(gFlyPath,"gFlyPath");
  ofile->Write();
  ofile->Close();

}
