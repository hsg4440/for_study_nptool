#include "TCanvas.h"
#include "TChain.h"
#include "TCutG.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TProfile.h"
#include <fstream>
#include <iostream>

// TChain* chain = NULL;


int CsI_Calib(){
//  TFile* fcuts[4][16];
//  TCutG* cuts[5][16];
//  for (int i = 0; i < 4; i++) {
//    for (int j = 0; j < 16; j++) {
//      // if (i != 3) {
//      fcuts[i][j] = new TFile(
//          Form("./CUTS/cuts_csi/CUT_%i_%i.root", i + 1, j + 1), "READ");
//      cuts[i][j]
//          = (TCutG*)fcuts[i][j]->FindObjectAny("CUTG");
//    }
//  }


	Int_t LocalMaxBinFit;
	Int_t LocalMinBinFit;
	Int_t GlobalMaxBinFit;
	Int_t GlobalMinBinFit;
	Int_t BinCounter = 0;
	Int_t NbBinMax;
	Int_t NbBin;

  
  TF1 *Gaus = new TF1("Gaus","gaus",-1000,1000);
	TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x^2",-1000,1000);
  TChain* chain_tmp = new TChain("PhysicsTree");
  chain_tmp->Add("../../Outputs/Analysis/r0153_cal_CsI.root");
  
  
  TFile* fcut[4][16];
  TCutG* cut[4][16];
  double params[4][16][3];
  for(int Telescope_Nb = 2; Telescope_Nb < 5; Telescope_Nb++){
    for(int CsI_Nb = 1; CsI_Nb < 17; CsI_Nb++){
      if(CsI_Nb < 17){
        std::cout << "Telescope : " << Telescope_Nb << "    CsI : " << CsI_Nb << std::endl;
        fcut[Telescope_Nb-1][CsI_Nb-1] = new TFile(Form("./CUT/CUT_CSI/CUT_%i_%i.root",Telescope_Nb,CsI_Nb));
        cut[Telescope_Nb-1][CsI_Nb-1] = (TCutG*)fcut[Telescope_Nb-1][CsI_Nb-1]->FindObjectAny("CUTG");
        cut[Telescope_Nb-1][CsI_Nb-1]->SetName(Form("CUT_%i_%i",CsI_Nb,Telescope_Nb));
        chain_tmp->Draw(Form("M2_ECsI_from_deltaE:CsI_E>>h_%i_%i(100,,,100,,)",CsI_Nb,Telescope_Nb),Form("TelescopeNumber == %i && CsI_N == %i &&CUT_%i_%i",Telescope_Nb,CsI_Nb,CsI_Nb,Telescope_Nb),"col"); 
	    
        TH2F *histogram = (TH2F*)fcut[Telescope_Nb-1][CsI_Nb-1]->FindObjectAny(Form("h_%i_%i", CsI_Nb,Telescope_Nb));
	      histogram->FitSlicesY(Gaus, -1000, 1000,0, "QG5",0);
	      TH1D *histogram_1 = (TH1D*)fcut[Telescope_Nb-1][CsI_Nb-1]->Get(Form("h_%i_%i_1", CsI_Nb,Telescope_Nb));
	      TH1D *histogram_0 = (TH1D*)fcut[Telescope_Nb-1][CsI_Nb-1]->Get(Form("h_%i_%i_0", CsI_Nb,Telescope_Nb));
	      TH1D *histogram_2 = (TH1D*)fcut[Telescope_Nb-1][CsI_Nb-1]->Get(Form("h_%i_%i_2", CsI_Nb,Telescope_Nb));

        LocalMaxBinFit = 1;
	      LocalMinBinFit = 1;
	      GlobalMaxBinFit = 1;
	      GlobalMinBinFit = 1;
        BinCounter = 0;
	      NbBinMax = 0;
	      NbBin = histogram_1->GetSize() - 2;
	   
        for(int k = 1; k <= NbBin; k++){
	      	std::cout <<histogram_2->GetBinContent(k)<< std::endl;
	      	if(histogram_2->GetBinContent(k) < 20){
          
	      		BinCounter += 1;
	      		if(LocalMaxBinFit - LocalMinBinFit < BinCounter){
	      			LocalMaxBinFit = k;
	      		}
	      	}
	      	else{
	      		if(BinCounter > NbBinMax){
	      			NbBinMax = BinCounter;
	      			GlobalMaxBinFit = LocalMaxBinFit;
	      			GlobalMinBinFit = LocalMinBinFit;
	      		}
	      		BinCounter = 0;
	      		LocalMaxBinFit = k;
	      		LocalMinBinFit = k+1;

	      	}
          
	      std::cout << "BinCounter : " << BinCounter << "\n";
	      std::cout << "NbBinMax : " << NbBinMax << "\n";
	      std::cout << "LocalMinBinFit " << LocalMinBinFit << "\n";
	      std::cout << "LocalMaxBinFit " << LocalMaxBinFit << "\n";
	      std::cout << "GlobalMinBinFit " << GlobalMinBinFit << "\n";
	      std::cout << "GlobalMaxBinFit " << GlobalMaxBinFit << "\n";
	    }
	      if(BinCounter > NbBinMax){
	      	NbBinMax = BinCounter;
	      	GlobalMaxBinFit = LocalMaxBinFit;
	      	GlobalMinBinFit = LocalMinBinFit;
	    }
	      std::cout << "BinCounter : " << BinCounter << "\n";
	      std::cout << "NbBinMax : " << NbBinMax << "\n";
	      std::cout << "LocalMinBinFit " << LocalMinBinFit << "\n";
	      std::cout << "LocalMaxBinFit " << LocalMaxBinFit << "\n";
	      std::cout << "GlobalMinBinFit " << GlobalMinBinFit << "\n";
	      std::cout << "GlobalMaxBinFit " << GlobalMaxBinFit << "\n";



	
	
	      histogram_1->Fit(f1,"","",histogram_1->GetBinCenter(GlobalMinBinFit),histogram_1->GetBinCenter(GlobalMaxBinFit));
       params[Telescope_Nb-1][CsI_Nb-1][0] = f1->GetParameter(0);
       params[Telescope_Nb-1][CsI_Nb-1][1] = f1->GetParameter(1);
       params[Telescope_Nb-1][CsI_Nb-1][2] = f1->GetParameter(2);


	
	      TFile *fout = new TFile(Form("./calibration/CsI/histograms/CsI_%i_MM%i.root",CsI_Nb,Telescope_Nb),"RECREATE");
	      histogram_0->Write("");
	      histogram_1->Write("");
	      histogram_2->Write("");
	      fout->Close();
      
      }
    }
    std::ofstream myfile;
    myfile.open(Form("./calibration/CsI/calib/CsI_MM%i.txt",Telescope_Nb),std::ios::in);
    for(int CsI = 1; CsI <=16; CsI++){
      myfile << Form("MUST2_T%i_CsI%i_E ",Telescope_Nb,CsI) << params[Telescope_Nb-1][CsI-1][0] << " " << params[Telescope_Nb-1][CsI-1][1] << " " << params[Telescope_Nb-1][CsI-1][2]  << "\n" ;
    }
    myfile.close();
  }
   
  return 0;
}

//int CsI_Calib_test(){
////  TFile* fcuts[4][16];
////  TCutG* cuts[5][16];
////  for (int i = 0; i < 4; i++) {
////    for (int j = 0; j < 16; j++) {
////      // if (i != 3) {
////      fcuts[i][j] = new TFile(
////          Form("./CUTS/cuts_csi/CUT_%i_%i.root", i + 1, j + 1), "READ");
////      cuts[i][j]
////          = (TCutG*)fcuts[i][j]->FindObjectAny("CUTG");
////    }
////  }
//
//
//	Int_t LocalMaxBinFit;
//	Int_t LocalMinBinFit;
//	Int_t GlobalMaxBinFit;
//	Int_t GlobalMinBinFit;
//	Int_t BinCounter = 0;
//	Int_t NbBinMax;
//	Int_t NbBin;
//
//  
//  TF1 *Gaus = new TF1("Gaus","gaus",-1000,1000);
//	TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x^2",-1000,1000);
//  
//  
//  TFile* fcut[2];
//  TCutG* cut[2];
//  fcut[0] = new TFile("./CUT/CUT_CSI/CUT_1_3.root");
//  fcut[1] = new TFile("./CUT/CUT_CSI/CUT_2_3.root");
//  cut[0] = (TCutG*)fcut[0]->FindObjectAny("CUTG");
//  cut[0]->SetName("CUT_1_3");
//  cut[1] = (TCutG*)fcut[1]->FindObjectAny("CUTG");
//  cut[1]->SetName("CUT_2_3");
//  //cut[1] = (TCutG*)fcut[1]->FindObjectAny("CUTG");
//  
//  TChain* chain_tmp = new TChain("PhysicsTree");
//  chain_tmp->Add("../../Outputs/Analysis/ana_r0161_calib.root");
//  
//  //auto c1 = new TCanvas();
//  //chain_tmp->Draw("M2_ECsI_from_deltaE:CsI_E>>h1(1000,,,1000,,)","TelescopeNumber == 3 && CsI_N == 1","col",1000000);
//  //auto c2 = new TCanvas();
//  //chain_tmp->Draw("M2_ECsI_from_deltaE:CsI_E>>h2(1000,,,1000,,)","TelescopeNumber == 3 && CsI_N == 1 &&CUT_1_3","col",1000000); 
//  //auto c3 = new TCanvas();
//  //chain_tmp->Draw("M2_ECsI_from_deltaE:CsI_E>>h3(1000,,,1000,,)","TelescopeNumber == 3 && CsI_N == 2","col",1000000);
//  auto c4 = new TCanvas();
//  chain_tmp->Draw("M2_ECsI_from_deltaE:CsI_E>>h4(100,,,100,,)","TelescopeNumber == 3 && CsI_N == 2 &&CUT_2_3","col"); 
//  
//	TH2F *histogram = (TH2F*)fcut[1]->FindObjectAny(Form("h%i", 4));
//	histogram->FitSlicesY(Gaus, -1000, 1000,0, "QG3",0);
//	TH1D *histogram_1 = (TH1D*)fcut[1]->Get(Form("h%i_1",4));
//	TH1D *histogram_0 = (TH1D*)fcut[1]->Get(Form("h%i_0",4));
//	TH1D *histogram_2 = (TH1D*)fcut[1]->Get(Form("h%i_2",4));
//	
//  
//  LocalMaxBinFit = 1;
//	LocalMinBinFit = 1;
//	GlobalMaxBinFit = 1;
//	GlobalMinBinFit = 1;
//	NbBinMax = 0;
//	NbBin = histogram_1->GetSize() - 2;
//
//	for(int k = 1; k <= NbBin; k++){
//		std::cout <<histogram_2->GetBinContent(k)<< std::endl;
//		if(TMath::Abs(histogram_2->GetBinContent(k)) < 10){
//		
//			BinCounter += 1;
//			if(LocalMaxBinFit < k){
//				LocalMaxBinFit = k;
//			}
//		}
//		else{
//			if(BinCounter > NbBinMax){
//				NbBinMax = BinCounter;
//				GlobalMaxBinFit = LocalMaxBinFit;
//				GlobalMinBinFit = LocalMinBinFit;
//			}
//			BinCounter = 0;
//			LocalMaxBinFit = k;
//			LocalMinBinFit = k + 1;
//
//		}
//	}
//
//	std::cout << GlobalMaxBinFit << "   " << GlobalMinBinFit << "   "<< histogram_1->GetSize() << "  " << histogram_1->GetBinError(1) << "  " << histogram_1->GetBinCenter(1)  << std::endl;
//
//	
//	
//	histogram_1->Fit(f1,"","",histogram_1->GetBinCenter(GlobalMinBinFit),histogram_1->GetBinCenter(GlobalMaxBinFit));
//  auto c5 = new TCanvas();
//  histogram_1->Draw("");
//	
//  std::ofstream myfile;
//  myfile.open ("./calibration/CsI/calib/CsI_MM3.txt");
//  myfile << "Writing this to a file.\n";
//  myfile << f1->GetParameter(0) << " " << f1->GetParameter(1) << " " << f1->GetParameter(2)  << "\n" ;
//  myfile.close();
//	TFile *fout = new TFile(Form("./calibration/CsI/histograms/CsI_%i_MM%i.root",2,3),"RECREATE");
//	histogram_1->Write("");
//	fout->Close();
//  return 0;
//}
//


////////////////////////////////////////////////////////////////////////////////
/*void LoadCutsCsI() {
  TFile* fcuts[4][16];
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 16; j++) {
      // if (i != 3) {
      fcuts[i][j] = new TFile(
          Form("./CUTS/cuts_csi/c_%i_%i.root", i + 1, j + 1), "READ");
      cuts[i][j]
          = (TCutG*)fcuts[i][j]->FindObjectAny(Form("c_%i_%i", i + 1, j + 1));
    }
    // else {
    //   cuts[i][j] = NULL;
    // }
    // }
  }
}

////////////////////////////////////////////////////////////////////////////////

void CreateChain() {
  TChain* chain_tmp = new TChain("PhysicsTree");
  // chain->Add("../nptool_outputs/Analysis/testNoCalibCsI.root");
  // chain->Add("../nptool_outputs/Analysis/CalibNoCsI14Nrun74.root");
  // chain->Add("../nptool_outputs/Analysis/CalibFiles/CalibNoCsI14Orun122_4_corrected.root");
  chain_tmp->Add("./NPOutput/PhyE744_0119_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0120_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0121_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0122_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0122_1_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0122_2_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0122_3_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0122_4_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0123_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0126_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0128_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0129_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0130_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0130_1_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0130_2_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0131_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_0_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_1_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_2_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_3_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_4_pixel.root");
  chain_tmp->Add("./NPOutput/PhyE744_0132_5_pixel.root");
  TMust2Physics*       M2 = new TMust2Physics();
  TBranch*             b_CsI_Ec;
  std::vector<double>* CsI_Ec = 0;
  TBranch*             b_CsI_th;
  double               CsI_th = 0;

  chain_tmp->SetBranchStatus("*", false);
  chain_tmp->SetBranchStatus("MUST2", true);
  chain_tmp->SetBranchStatus("CsI_th", true);
  chain_tmp->SetBranchAddress("MUST2", &M2);
  chain_tmp->SetBranchAddress("CsI_th", &CsI_th);
  // chain_tmp->SetBranchAddress("CsI_Ec", &CsI_Ec, &b_CsI_Ec);

  TFile* file = new TFile("CsIth_only.root", "RECREATE");
  TTree* T    = new TTree("PhysicsTree", "CsIthOnly");

  double csith, csiE, siE, csiN, six, siy, tele;
  T->Branch("CsI_th", &csith);
  T->Branch("CsI_E", &csiE);
  T->Branch("Si_E", &siE);
  T->Branch("CsI_N", &csiN);
  T->Branch("Si_X", &six);
  T->Branch("Si_Y", &siy);
  T->Branch("TelescopeNumber", &tele);

  unsigned int N = chain_tmp->GetEntries();

  for (int i = 0; i < N; i++) {
    chain_tmp->GetEntry(i);

    // if (i % N == 0)
    printf("%f\n", i / (double)N * 100);
    if (M2->Si_E.size() == 1 && M2->CsI_E[0] > 0) {
      csith = CsI_th;
      csiE  = M2->CsI_E[0];
      siE   = M2->Si_E[0];
      csiN  = M2->CsI_N[0];
      six   = M2->Si_X[0];
      siy   = M2->Si_Y[0];
      tele  = M2->TelescopeNumber[0];
      T->Fill();
    }
  }
  file->Write();
}

void LoadChain() {

  // CreateChain();
  chain = new TChain("PhysicsTree");
  chain->Add("./CsIth_only.root");
}

////////////////////////////////////////////////////////////////////////////////

void DrawCsIth_CsI_E() {

  LoadChain();
  LoadCutsCsI();

  int N_Telescope = 4;
  int N_CsI       = 16;

  TCanvas* canvas[N_Telescope];

  TFile* f = new TFile("./root_files/CsIthvsCsIE.root", "recreate");

  for (int i = 0; i < 4; i++) {
    canvas[i] = new TCanvas();
    canvas[i]->Divide(4, 4);
    for (int j = 0; j < 16; j++) {
      std::cout << i + 1 << " " << j + 1 << std::endl;
      canvas[i]->cd(j + 1);
      // cuts[i];
      // // if (i == 1 && j == 15) {
      //   chain->Draw(
      //       Form("CsI_th:CsI_E>>T%i_CsIth_CsI_%i(500,8000,10000, 200,0,40)",
      //            i + 1, j + 1),
      //       Form("TelescopeNumber==%i&&CsI_N==%i && "
      //            "c_%i_%i&&ThetaLab<11&&TelescopeNumber@.size()==1",
      //            i + 1, j + 1, i + 1, j + 1),
      //       "");
      // } else {
      chain->Draw(
          Form("CsI_th:CsI_E>>T%i_CsIth_CsI_%i(500,8000,10000, 200,0,40)",
               i + 1, j + 1),
          Form("TelescopeNumber==%i&&CsI_N==%i && "
               "c_%i_%i&&TelescopeNumber@.size()==1",
               i + 1, j + 1, i + 1, j + 1),
          "");
      // }
    }
  }
  f->Write();
}

void DrawCsIth_CsI_E_Simu() {

  TFile* file        = new TFile("NPOutput/Simu2CalibAnalysis.root", "OPEN");
  TTree* PhysicsTree = (TTree*)file->FindObjectAny("PhysicsTree");
  int    N_Telescope = 4;
  int    N_CsI       = 16;

  TCanvas* canvas[N_Telescope];

  TFile* f = new TFile("./root_files/CsIthvsCsIE_Simu.root", "recreate");

  for (int i = 0; i < 4; i++) {
    canvas[i] = new TCanvas();
    canvas[i]->Divide(4, 4);
    for (int j = 0; j < 16; j++) {
      std::cout << i + 1 << " " << j + 1 << std::endl;
      canvas[i]->cd(j + 1);
      // cuts[i];
      PhysicsTree->Draw(
          Form("CsI_th:CsI_E>>T%i_CsIth_CsI_%i(500,8000,10000, 200,0,40)",
               i + 1, j + 1),
          Form("TelescopeNumber==%i&&CsI_N==%i && "
               "TelescopeNumber@.size()==1",
               i + 1, j + 1),
          "");
    }
  }
  f->Write();
}


void InterpolCsI_E() {

  TFile*   FileCsIth   = new TFile("./root_files/CsIthvsCsIE.root", "READ");
  TFile*   FileInter   = new TFile("./root_files/CsI_Inter.root", "recreate");
  int      N_Telescope = 4;
  int      N_CsI       = 16;
  TCanvas* canvas[N_Telescope];
  TCanvas* c_residue[N_Telescope];

  for (int i = 0; i < N_Telescope; i++) {
    canvas[i] = new TCanvas();
    canvas[i]->Divide(4, 4);

    TH2F*     gh[4][16];
    TProfile* gh_ProfileX[4][16];

    for (int j = 0; j < N_CsI; j++) {
      canvas[i]->cd();
      canvas[i]->cd(j + 1);

      gh[i][j] = (TH2F*)FileCsIth->FindObjectAny(
          Form("T%i_CsIth_CsI_%i", i + 1, j + 1));
      // gh[i][j]->Rebin2D(2, 2);

      unsigned int nbinsX = gh[i][j]->GetXaxis()->GetNbins();
      unsigned int nbinsY = gh[i][j]->GetYaxis()->GetNbins();
      for (unsigned int r = 0; r < nbinsX; r++) {
        for (unsigned int l = 0; l < nbinsY; l++) {
          int counts = gh[i][j]->GetBinContent(r, l);
          if (counts < 4) {
            if (counts > 0)
              gh[i][j]->SetBinContent(r, l, 0, 0);
          }
        }
      }

      gh[i][j]->ProfileX(Form("hT%i_ProfileX%i", i + 1, j + 1));

      gh_ProfileX[i][j] = (TProfile*)gROOT->FindObjectAny(
          Form("hT%i_ProfileX%i", i + 1, j + 1));
      gh_ProfileX[i][j]->Draw();

      TGraph*      gCal = new TGraph();
      unsigned int Nx   = gh_ProfileX[i][j]->GetXaxis()->GetNbins();
      int          iter = 0;
      for (unsigned int k = 0; k < Nx; k++) {
        double binCenter = gh_ProfileX[i][j]->GetXaxis()->GetBinCenter(k);
        double val       = gh_ProfileX[i][j]->GetBinContent(k);
        if (val > 0) {
          gCal->SetPoint(iter, binCenter, val);
          iter++;
        }
      }
      gCal->SetName(Form("g%i_%i", i + 1, j + 1));
      gCal->Write();
    } // Loop CsI
  } // Loop Telescope
  FileInter->Write();
}

void FitCsI_E(string FIT_TYPE) {

  TFile*   FileCsIth   = new TFile("./root_files/CsIthvsCsIE.root", "READ");
  int      N_Telescope = 4;
  int      N_CsI       = 16;
  TCanvas* canvas[N_Telescope];
  TCanvas* c_residue[N_Telescope];

  for (int i = 0; i < N_Telescope; i++) {
    canvas[i] = new TCanvas();
    canvas[i]->Divide(4, 4);

    c_residue[i] = new TCanvas;
    c_residue[i]->Divide(4, 4);
    TH2F*     gh[5][16];
    TProfile* gh_ProfileX[5][16];
    TH2F*     gh_PX[5][16];

    double range_inf = 8000;
    double range_sup = 10000;

    TF1* f1;
    TF1* fpol1;
    TF1* fpol2;
    TF1* fpol3;
    TF1* fpol4;

    TF1* f2[5][16];

    fpol1 = new TF1("fpol1", "[0] + [1]*x", range_inf, range_sup);
    fpol2 = new TF1("fpol2", "[0] + [1]*x + [2]*x*x", range_inf, range_sup);
    fpol3 = new TF1("fpol3", "[0] + [1]*x + [2]*x*x + [3]*x*x*x", range_inf,
                    range_sup);
    fpol4 = new TF1("fpol4", "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x",
                    range_inf, range_sup);

    // fpol1->SetParLimits(0, 0, 10);

    double   p_pol1[2] = {0, 0};
    double   p_pol2[3] = {0, 0, 0};
    double   p_pol3[4] = {0, 0, 0, 0};
    double   p_pol4[5] = {0, 0, 0, 0, 0};
    ofstream Calib[5];

    Calib[i].open(Form("Cal_CsI_E_MM%i.cal", i + 1));
    for (int j = 0; j < N_CsI; j++) {
      canvas[i]->cd();
      canvas[i]->cd(j + 1);

      if (FIT_TYPE == "pol1") {
        f2[i][j] = new TF1(Form("f2%i", j + 1), "[0] + [1]*x", 0, 10000);
      } else if (FIT_TYPE == "pol2") {
        f2[i][j]
            = new TF1(Form("f2%i", j + 1), "[0] + [1]*x + [2]*x*x", 0, 10000);
      } else if (FIT_TYPE == "pol3") {
        f2[i][j] = new TF1(Form("f2%i", j + 1),
                           "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, 10000);
      } else if (FIT_TYPE == "pol4") {
        f2[i][j] = new TF1(Form("f2%i", j + 1),
                           "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x", 0,
                           10000);
      }

      gh[i][j] = (TH2F*)gDirectory->FindObjectAny(
          Form("T%i_CsIth_CsI_%i", i + 1, j + 1));
      gh[i][j]->Rebin2D(2, 2);
      // gh[i][j]->SetMinimum(100);
      // gh[i][j]->Draw();
      // gPad->WaitPrimitive();
      unsigned int nbinsX = gh[i][j]->GetXaxis()->GetNbins();
      unsigned int nbinsY = gh[i][j]->GetYaxis()->GetNbins();
      for (unsigned int r = 0; r < nbinsX; r++) {
        for (unsigned int l = 0; l < nbinsY; l++) {
          int counts = gh[i][j]->GetBinContent(r, l);
          // if (counts < 1 && counts > 0) {
          //   gh[i][j]->SetBinContent(r, l, 0, 0);
          // }
        }
      }

      gh[i][j]->ProfileX(Form("hT%i_ProfileX%i", i + 1, j + 1));
      gh[i][j]->Clone(Form("hT%i_PX%i", i + 1, j + 1));

      gh_ProfileX[i][j] = (TProfile*)gROOT->FindObjectAny(
          Form("hT%i_ProfileX%i", i + 1, j + 1));

      gh_PX[i][j]
          = (TH2F*)gROOT->FindObjectAny(Form("hT%i_PX%i", i + 1, j + 1));
      gh_PX[i][j]->GetXaxis()->SetRangeUser(0, 10000);

      if (gh_PX[i][j]->GetEntries() > 0) {

        if (FIT_TYPE == "pol1") {
          gh_PX[i][j]->Fit(fpol1, "Q", "", range_inf, range_sup);
          fpol1->GetParameters(&p_pol1[0]);
        } else if (FIT_TYPE == "pol2") {
          gh_PX[i][j]->Fit(fpol1, "Q", "", range_inf, range_sup);
          fpol1->GetParameters(&p_pol1[0]);
          fpol2->SetParameters(p_pol1[0], p_pol1[1], 0);
          fpol2->SetParLimits(1, 0, p_pol1[1] + 5);
          gh_PX[i][j]->Fit(fpol2, "Q", "", range_inf, range_sup);
          fpol2->GetParameters(&p_pol2[0]);
        } else if (FIT_TYPE == "pol3") {
          gh_PX[i][j]->Fit(fpol1, "Q", "", range_inf, range_sup);
          fpol1->GetParameters(&p_pol1[0]);
          fpol2->SetParameters(p_pol1[0], p_pol1[1], 0);
          fpol2->SetParLimits(1, 0, p_pol1[1] + 5);
          gh_PX[i][j]->Fit(fpol2, "Q", "", range_inf, range_sup);
          fpol2->GetParameters(&p_pol2[0]);
          fpol3->SetParameters(p_pol2[0], p_pol2[1], p_pol2[2], 0);
          gh_PX[i][j]->Fit(fpol3, "Q", "", range_inf, range_sup);
          fpol3->GetParameters(&p_pol3[0]);
        } else if (FIT_TYPE == "pol4") {
          gh_PX[i][j]->Fit(fpol1, "Q", "", range_inf, range_sup);
          fpol1->GetParameters(&p_pol1[0]);
          fpol2->SetParameters(p_pol1[0], p_pol1[1], 0);
          fpol2->SetParLimits(1, 0, p_pol1[1] + 5);
          gh_PX[i][j]->Fit(fpol2, "Q", "", range_inf, range_sup);
          fpol2->GetParameters(&p_pol2[0]);
          fpol3->SetParameters(p_pol2[0], p_pol2[1], p_pol2[2], 0);
          gh_PX[i][j]->Fit(fpol3, "Q", "", range_inf, range_sup);
          fpol3->GetParameters(&p_pol3[0]);
          fpol4->SetParameters(p_pol3[0], p_pol3[1], p_pol3[2], p_pol3[3], 0);
          gh_PX[i][j]->Fit(fpol4, "Q", "", range_inf, range_sup);
          fpol4->GetParameters(&p_pol4[0]);
        }

        if (FIT_TYPE == "pol1") {
          f2[i][j]->SetParameters(&p_pol1[0]);
        } else if (FIT_TYPE == "pol2") {
          f2[i][j]->SetParameters(&p_pol2[0]);
        } else if (FIT_TYPE == "pol3") {
          f2[i][j]->SetParameters(&p_pol3[0]);
        } else if (FIT_TYPE == "pol4") {
          f2[i][j]->SetParameters(&p_pol4[0]);
        }

        gh_PX[i][j]->Draw();
        // h[j]->Draw("col");
        f2[i][j]->SetLineColor(kGreen);
        f2[i][j]->Draw("same");

        if (FIT_TYPE == "pol1") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol1[0] << " " << p_pol1[1] << " " << endl;
        } else if (FIT_TYPE == "pol2") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol2[0] << " " << p_pol2[1] << " " << p_pol2[2]
                   << endl;
        } else if (FIT_TYPE == "pol3") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol3[0] << " " << p_pol3[1] << " " << p_pol3[2]
                   << " " << p_pol3[3] << endl;
        } else if (FIT_TYPE == "pol4") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol4[0] << " " << p_pol4[1] << " " << p_pol4[2]
                   << " " << p_pol4[3] << " " << p_pol4[4] << endl;
        }
      } else {
        if (FIT_TYPE == "pol1") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol1[0] << " " << p_pol1[1] << " " << endl;
        } else if (FIT_TYPE == "pol2") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol2[0] << " " << p_pol2[1] << " " << p_pol2[2]
                   << endl;
        } else if (FIT_TYPE == "pol3") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol3[0] << " " << p_pol3[1] << " " << p_pol3[2]
                   << " " << p_pol3[3] << endl;
        } else if (FIT_TYPE == "pol4") {
          Calib[i] << "MUST2_T" << i + 1 << "_CsI" << j + 1 << "_E"
                   << " " << p_pol4[0] << " " << p_pol4[1] << " " << p_pol4[2]
                   << " " << p_pol4[3] << " " << p_pol4[4] << endl;
        }
      }

      std::vector<double> v_iterator;
      std::vector<double> v_residue;
      // for (int k = 8200; k < 9000; k++) {
      for (int k = range_inf; k < range_sup + 500; k++) {
        double fit_val   = f2[i][j]->Eval(k);
        double histo_bin = gh_ProfileX[i][j]->GetXaxis()->FindBin(k);
        double histo_val = gh_ProfileX[i][j]->GetBinContent(histo_bin);
        if (histo_val != 0) {
          // v_iterator.push_back(k);
          double calib_val = 0;
          if (FIT_TYPE == "pol1") {
            calib_val = p_pol1[0] + p_pol1[1] * k;
          } else if (FIT_TYPE == "pol2") {
            calib_val = p_pol2[0] + p_pol2[1] * k + p_pol2[2] * k * k;
          } else if (FIT_TYPE == "pol3") {
            calib_val = p_pol3[0] + p_pol3[1] * k + p_pol3[2] * k * k
                        + p_pol3[3] * k * k * k;
          } else if (FIT_TYPE == "pol4") {
            calib_val = p_pol4[0] + p_pol4[1] * k + p_pol4[2] * k * k
                        + p_pol4[3] * k * k * k + p_pol4[4] * k * k * k * k;
          }
          if (calib_val < 26 && calib_val > 6) {
            v_iterator.push_back(calib_val);
            v_residue.push_back((histo_val - fit_val) / calib_val * 100);
            if (abs((fit_val - histo_val) / calib_val * 100) > 10) {
              std::cout << "Telescope_" << i + 1 << " Csi_" << j + 1 << " FIT "
                        << fit_val << " HISTO " << histo_val << " Residue :  "
                        << ((histo_val - fit_val) / calib_val * 100) << " %"
                        << std::endl;
            }
          }
        }
      }
      c_residue[i]->cd(j + 1);
      c_residue[i]->SetName(FIT_TYPE.c_str());
      c_residue[i]->SetTitle(FIT_TYPE.c_str());
      TGraph* g1 = new TGraph(v_iterator.size(), &v_iterator[0], &v_residue[0]);
      double  min = g1->GetXaxis()->GetXmin();
      double  max = g1->GetXaxis()->GetXmax();
      TLine*  l   = new TLine(min, 0, max, 0);
      l->SetLineWidth(1);
      g1->SetLineWidth(2);
      g1->SetLineColor(kRed);
      g1->SetTitle(Form("T_%i CsI_%i", i + 1, j + 1));
      // g1->GetXaxis()->SetLabelSize(0.07);
      g1->GetYaxis()->SetLabelSize(0.07);
      // g1->GetXaxis()->SetLabelOffset(0.03);
      g1->GetYaxis()->SetLabelOffset(0.02);
      // g1->GetYaxis()->SetRangeUser(-4, 4);
      g1->Draw("");
      l->Draw("same");
    } // Loop CsI
    Calib[i].close();
  } // Loop Telescope
}
*/
