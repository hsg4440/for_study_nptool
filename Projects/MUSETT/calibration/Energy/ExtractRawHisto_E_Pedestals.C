////////////////////////////////////////////////////////////////////////////
// This macro takes a converted ROOT file from the GANIL e569 exp. and    //
// create an histogram for each strip (X and Y) of the MUST2 array filled //
// with the energy. The histograms are dumped in an output ROOT file.     //
//                                                                        //
// This is especially usefull for calibration purposes when there is no   //
// need to work directly on the TChain                                     //
////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"
#include "TBranch.h"
#include "TH1F.h"

#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Detectors/MUST2/TMust2Data.h"
#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/Detectors/Mugast/TMugastData.h"
#include "/home/mugast/Analysis_MUGAST_LISE_test/nptool/NPLib/include/RootInput.h"

#define NBTELESCOPE 4	
#define	NBSTRIPS	128
#define NBSILI     0

class pedestal{
  public:
    pedestal(){};
    ~pedestal(){};
  public: 
    map<int, vector<int> > pedestalX;
    map<int, vector<int> > pedestalY;
  
    int get_pedestalX(int det,int strip){
      return pedestalX[det][strip-1];
      }
    int get_pedestalY(int det,int strip){
      return pedestalY[det][strip-1];
      }

  public:
    int SubstractX(int det, int strip,pedestal other){
        return (get_pedestalX(det,strip)-other.get_pedestalX(det,strip));
      }
    
    int SubstractY(int det, int strip,pedestal other){
        return (get_pedestalY(det,strip)-other.get_pedestalY(det,strip));
      }


    void print(){
      cout << "X" << endl;
      for(auto& it: pedestalX){
        cout << " Telescope: " << it.first << endl;
        unsigned int size = it.second.size();
        for(unsigned int i = 0 ; i < size ; i++){
          std::cout << " strip: " << i+1 << " pedestal: " << it.second[i] << "\n"; 
        } 
      }

      cout << "Y" << endl;
      for(auto& it: pedestalY){
        cout << " Telescope: " << it.first << endl;
        unsigned int size = it.second.size();
        for(unsigned int i = 0 ; i < size ; i++){
          std::cout << " strip: " << i+1 << " pedestal: " << it.second[i] << "\n"; 
        } 
      }
    }

    void LoadPedestals(const std::string& path, int num){
      ifstream file(path.c_str());
      std::string buffer;
      // ignore the first line of comments
      getline(file,buffer);
      int pedestalE,pedestalT;
      int mate=0;
      int count=0;

      while (file >> pedestalE >> pedestalT){
        if(count++%16==0){
          mate++; 
        }
        if(mate<9){ // X
          pedestalX[num].push_back(pedestalE);
        } 
        else if(mate>8&&mate<17){ // Y
          pedestalY[num].push_back(pedestalE);
        }
        else{ // SiLi/CsI
          break; 
        }
      }
    }
};
////////////////////////////////////////////////////////////////////////////////
void Test(){
  pedestal q;
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI3_Sector1",1);

  pedestal p;
  p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI3_Sector1",1);
 
  for(unsigned int i = 0 ; i < 128 ; i++){
    cout << p.SubstractX(1,i+1,q) << " " << p.SubstractY(1,i+1,q)  << endl;
  }
}
///////////////////////////////////////////////////////////////////////////////
void ExtractRawHistos(const char* fname = "run_XXXX_X", std::string fDet="MUST2")
{

  pedestal p,q;
  if(fDet=="MUST2"){
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector1.pedestal",1);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector2.pedestal",2);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector3.pedestal",3);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector4.pedestal",4); 

   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector1.pedestal",1);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector2.pedestal",2);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector3.pedestal",3);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/test2023_sourceON_sector4.pedestal",4); 

   // q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI3_Sector1",1);
   // q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI3_Sector2",2);
   // q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI3_Sector3",3);
   // q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI3_Sector4",4); 

    }

 else if(fDet=="Mugast"){
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI1_Sector1",1);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI1_Sector2",3);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI1_Sector3",4);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI1_Sector4",5); 
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI2_Sector1",7);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI2_Sector2",9);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI2_Sector3",10);
   p.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/run165/pedestal_das_run165_MUVI2_Sector4",11); 

   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI1_Sector1",1);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI1_Sector2",3);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI1_Sector3",4);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI1_Sector4",5); 
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI2_Sector1",7);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI2_Sector2",9);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI2_Sector3",10);
   q.LoadPedestals("/home/mugast/ganacq_manip/mugast/das-save/pedestal-save/pedestal_read_3april_15h18/pedestal_MUVI2_Sector4",11); 
 
    }

  else{
    cout << "Error in detector name " << fDet << endl;
    exit(1);
    }


  RootInput* Input = RootInput::getInstance("RunToTreat.txt");
  TChain* Chain = Input->GetChain();
  // connect the TfDet.c_str()Data branch  
  Chain->SetBranchStatus(fDet.c_str(),true);
  Chain->SetBranchStatus("fMM*",true);

  TMust2Data *rawMust2;
  TMugastData *rawMugast;
  if (fDet=="MUST2"){
    rawMust2 = new TMust2Data();
    Chain->SetBranchAddress(fDet.c_str(),&rawMust2);
  }
  else if (fDet=="Mugast"){
    rawMugast = new TMugastData();
    Chain->SetBranchAddress(fDet.c_str(),&rawMugast);
  }

  // open the output ROOT file
  TString outFileName = "./Histograms/";
  outFileName += fname;
  outFileName += "_Raw" + fDet + "Histos_CorrectionPedestals.root";
  cout<< outFileName<< endl;
  TFile *outFile = new TFile(outFileName, "recreate");

  // prepare output histograms for Must2
  TH1F* hStripXEnergy[NBTELESCOPE][NBSTRIPS];
  TH1F* hStripYEnergy[NBTELESCOPE][NBSTRIPS];
  //from i=4 because I wanted only the T5,T8
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSTRIPS; j++) {
      // strips XE
      TString hnameXE     = Form("hMM%d_STRX_E%d", i+1, j+1);
      TString htitleXE    = Form("MM%d_STRX_E%d", i+1, j+1);
      hStripXEnergy[i][j] = new TH1F(hnameXE, htitleXE, 16384, 0, 16384);
      // strips YE
      TString hnameYE     = Form("hMM%d_STRY_E%d", i+1, j+1);
      TString htitleYE    = Form("MM%d_STRY_E%d", i+1, j+1);
      hStripYEnergy[i][j] = new TH1F(hnameYE, htitleYE, 16384, 0, 16384);
    }
  }
  TH1F* hSiLiEnergy[NBTELESCOPE][NBSILI];
  for (Int_t i = 0; i < NBTELESCOPE; i++) {
    for (Int_t j = 0; j < NBSILI; j++) {
      TString hnameSiLiE     = Form("hMM%d_SILI_E%d", i+1, j+1);
      TString htitleSiLiE    = Form("MM%d_SILI_E%d", i+1, j+1);
      hSiLiEnergy[i][j] = new TH1F(hnameSiLiE, htitleSiLiE, 16384, 0, 16384);
    }
  }

  // read the data
  Int_t nentries = Chain->GetEntries();
  cout << "TChain contains " << nentries << " events" << endl;
  for (Int_t i = 0; i < nentries; i++) {
    if (i%1000 == 0) cout << "\rEntry " << i << "\t" << i/(double)nentries * 100 << " %" << flush;
    ///////////////////////
    // read Must2 branch //
    ///////////////////////
    //brMust2->GetEntry(i);
    Chain->GetEntry(i);
    // get StripXE multiplicity
    Int_t multXE = 0; 
    if (fDet=="MUST2"){
      multXE = rawMust2->GetMMStripXEMult();
    }else if (fDet=="Mugast"){
      multXE = rawMugast->GetDSSDXEMult();
    }
    //cout<< "multXE" << multXE<< endl;
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    //
    Int_t det = 0;   
    Int_t stripX = 0;
    Int_t energy = 0;

    for (Int_t j = 0; j < multXE; j++) {
      if (fDet=="MUST2"){
        det    = rawMust2->GetMMStripXEDetectorNbr(j);
        stripX = rawMust2->GetMMStripXEStripNbr(j);
        energy = rawMust2->GetMMStripXEEnergy(j);
      } else if (fDet=="Mugast"){
        det    = rawMugast->GetDSSDXEDetectorNbr(j);
        stripX = rawMugast->GetDSSDXEStripNbr(j);
        energy = rawMugast->GetDSSDXEEnergy(j);
      }
      if ((det < NBTELESCOPE+1) && (stripX < NBSTRIPS+1)) {
        std::cout << ">>>>>>>>>> DEBUG : " << 0 << std::endl;
        std::cout << det << " " << stripX << " " << q  << std::endl;
        hStripXEnergy[det-1][stripX-1]->Fill(energy + p.SubstractX(det, stripX, q));
        std::cout << ">>>>>>>>>> DEBUG : " << 1 << std::endl;
      } else {
        cout << "Error filling histograms: X_E" << endl;
        if (fDet=="MUST2"){
          rawMust2->Dump();
        } else if (fDet=="Mugast"){
          rawMugast->Dump();
        }
      }
    }
    std::cout << ">>>>>>>>>> DEBUG : " << 2 << std::endl;
    // get StripYE multiplicity
    Int_t multYE=0;
    if (fDet=="MUST2"){
      multYE = rawMust2->GetMMStripYEMult();
    } else if (fDet=="Mugast"){
      multYE = rawMugast->GetDSSDYEMult();
    }
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms

    det = 0;   
    Int_t stripY = 0;
    energy = 0;
    for (Int_t j = 0; j < multYE; j++) {
      if (fDet=="MUST2"){
        det    = rawMust2->GetMMStripYEDetectorNbr(j);
        stripY = rawMust2->GetMMStripYEStripNbr(j);
        energy = rawMust2->GetMMStripYEEnergy(j);
      }else if (fDet=="Mugast"){
        det    = rawMugast->GetDSSDYEDetectorNbr(j);
        stripY = rawMugast->GetDSSDYEStripNbr(j);
        energy = rawMugast->GetDSSDYEEnergy(j);
      }

      if ((det < NBTELESCOPE+1) && (stripY < NBSTRIPS+1)) {
        hStripYEnergy[det-1][stripY-1]->Fill(energy + p.SubstractY(det, stripY, q));
      }
      else {
        cout << "Error filling histograms: Y_E" << endl;
        if (fDet=="MUST2"){
          rawMust2->Dump();
        }else if (fDet=="Mugast"){
          rawMugast->Dump();
        }
      }
    }
    // get SiLiE multiplicity

    //Int_t multSiLiE = rawMust2->GetMMSiLiEMult();
    // loop on StripXE multiplicity and get information concerning the event
    // fill histograms
    /* for (Int_t j = 0; j < multSiLiE; j++) {
       Int_t det    = rawMust2->GetMMSiLiEDetectorNbr(j);
       Int_t pad = rawMust2->GetMMSiLiEPadNbr(j);
       Int_t energy = rawMust2->GetMMSiLiEEnergy(j);
       if ((det < NBTELESCOPE+1) && (pad < NBSTRIPS+1)) {
       hSiLiEnergy[det-1][pad-1]->Fill(energy);
       }
       else {
       cout << "Error filling histograms: Must2_SiLi_E" << endl;
       rawMust2->Dump();
       }
       }*/

  }
  cout << '\n';

  // write on disk output file and close it
  outFile->Write();
  outFile->Close();
}




