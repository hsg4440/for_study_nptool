TH1F* fitPSp = 0;
TH1F* fitExA = 0;
TH1F* fitExB = 0;
TH1F* fitExC = 0;
TH1F* fitExD = 0;


// REMEMBER TO CHANGE FIT PARAMETER NUMNER!!!!!!!!!!!!!!!!!!!!!!!!!1
double parryFitFunc(Double_t *x, Double_t *par){
  //int simulate_rebin = 8;
  double result =0;
  //result+= par[0]*fitPSp->GetBinContent((int)(x[0]));///simulate_rebin)+1);
  //result+= par[1]*fitExA->GetBinContent((int)(x[0]));///simulate_rebin)+1);
  //result+= par[2]*fitExB->GetBinContent((int)(x[0]));///simulate_rebin)+1);
  //result+= par[3]*fitExC->GetBinContent((int)(x[0]));///simulate_rebin)+1);
  //result+= par[4]*fitExD->GetBinContent((int)(x[0]));///simulate_rebin)+1);

  result+= fitPSp->GetBinContent((int)(x[0])) / par[0];///simulate_rebin)+1);
  result+= fitExA->GetBinContent((int)(x[0])) / par[1];///simulate_rebin)+1);
  result+= fitExB->GetBinContent((int)(x[0])) / par[2];///simulate_rebin)+1);
  result+= fitExC->GetBinContent((int)(x[0])) / par[3];///simulate_rebin)+1);
  result+= fitExD->GetBinContent((int)(x[0])) / par[4];///simulate_rebin)+1);




  return result;
}


void FillParryHistograms(double binning){
  double bw;
  fitPSp = 0;
  fitExA = 0; fitExB = 0;
  fitExC = 0; fitExD = 0;

  TCanvas* tempcanv = new TCanvas("tempcanv","tempcanv",100,100);

  TFile *fPS = new TFile("PhaseSpaceForParryFit_14Jul23.root","READ");
    fitPSp = (TH1F*)fPS->Get("hfitPSp");
    bw = fitPSp->GetXaxis()->GetBinWidth(2);
    if(bw!=binning){
      fitPSp->Rebin(binning/bw);
    }

  TFile *fA  = new TFile("../../../Outputs/Analysis/Sim_18Oct22_47Kdp_4850.root","READ");
    TTree* treeA = (TTree*) fA->FindObjectAny("PhysicsTree");
    treeA->Draw("Ex>>hfitA(10000,-2,8)","Mugast.TelescopeNumber>0","");
    fitExA = (TH1F*)gDirectory->Get("hfitA");
    bw = fitExA->GetXaxis()->GetBinWidth(2);
    if(bw!=binning){
      fitExA->Rebin(binning/bw);
    }

  TFile *fB  = new TFile("../../../Outputs/Analysis/Sim_18Oct22_47Kdp_5150.root","READ");
    TTree* treeB = (TTree*) fB->FindObjectAny("PhysicsTree");
    treeB->Draw("Ex>>hfitB(10000,-2,8)","Mugast.TelescopeNumber>0","");
    fitExB = (TH1F*)gDirectory->Get("hfitB");
    bw = fitExB->GetXaxis()->GetBinWidth(2);
    if(bw!=binning){
      fitExB->Rebin(binning/bw);
    }

  TFile *fC  = new TFile("../../../Outputs/Analysis/Sim_18Oct22_47Kdp_5750.root","READ");
    TTree* treeC = (TTree*) fC->FindObjectAny("PhysicsTree");
    treeC->Draw("Ex>>hfitC(10000,-2,8)","Mugast.TelescopeNumber>0","");
    fitExC = (TH1F*)gDirectory->Get("hfitC");
    bw = fitExC->GetXaxis()->GetBinWidth(2);
    if(bw!=binning){
      fitExC->Rebin(binning/bw);
    }

  TFile *fD  = new TFile("../../../Outputs/Analysis/Sim_18Oct22_47Kdp_6050.root","READ");
    TTree* treeD = (TTree*) fD->FindObjectAny("PhysicsTree");
    treeD->Draw("Ex>>hfitD(10000,-2,8)","Mugast.TelescopeNumber>0","");
    fitExD = (TH1F*)gDirectory->Get("hfitD");
    bw = fitExD->GetXaxis()->GetBinWidth(2);
    if(bw!=binning){
      fitExD->Rebin(binning/bw);
    }

  delete tempcanv;

  cout << "LOADED SIMULATION FILES FOR UNBOUND PARRY FIT" << endl;

}


void FitUnboundParryMethod(TH1F* hist, double binning){
  

  FillParryHistograms(binning);


cout << "hereStart" << endl;
  string settings = "RBSN";

cout << "hereBeforeDefineFunc" << endl;
  TF1* func = new TF1("fitUnbnd",parryFitFunc, 4.0, 7.5, 5);
  func->SetParameter(0,5.);
  func->SetParLimits(0,0.,100.);
  func->SetParameter(1,5.);
  func->SetParLimits(1,0.,100.);
  func->SetParameter(2,5.);
  func->SetParLimits(2,0.,100.);
  func->SetParameter(3,5.);
  func->SetParLimits(3,0.,100.);
  func->SetParameter(4,5.);
  func->SetParLimits(4,0.,100.);



cout << "hereBeforeFit" << endl;
  TFitResultPtr fit = hist->Fit(func, settings.c_str(), "", 4.0, 7.5);


  hist->SetLineColor(kBlack);
  hist->SetLineWidth(2);
  hist->Draw("hist");
  fitPSp->Scale(1./fit->Parameter(0));
  fitPSp->SetLineColor(kRed);
  fitPSp->Draw("same hist");
  fitExA->Scale(1./fit->Parameter(1));
  fitExA->SetLineColor(kOrange);
  fitExA->Draw("same hist");
  fitExB->Scale(1./fit->Parameter(2));
  fitExB->SetLineColor(kGreen);
  fitExB->Draw("same hist");
  fitExC->Scale(1./fit->Parameter(3));
  fitExC->SetLineColor(kBlue);
  fitExC->Draw("same hist");
  fitExD->Scale(1./fit->Parameter(4));
  fitExD->SetLineColor(kViolet);
  fitExD->Draw("same hist");

}

