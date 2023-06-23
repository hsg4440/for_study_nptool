string reactionName; /* defined by choice of dp or dt */
#include "DefineColours.h"
#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

//#include "CS2_dt.h"
#include "CS2_master.h"
//#include "ThreeBodyBreakup.h"
//#include "ThreeBodyBreakup_FitPhaseSpace.h"


void AddGammaLines(TH1F* hist, double particle, double ymax){
//  string base = "sub ";
//
//  for(int i=1; i<means.size();i++){
//    string name = base + to_string(means.at(i));
//    TLine *line = new TLine(particle-means.at(i), 0.0, particle-means.at(i), ymax);
//    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
//    line->Draw();
//    TText *text = new TText((1.-(means.at(i)/particle))*particle,0.8*ymax,name.c_str());
//    text->SetTextAngle(90);
//    //text->SetTextSize(40);
//    text->Draw();
//  }
}

void AddPlacedGammas(TH1F* hist, double ymax){
//  hist->Draw();
//  for(int i=0; i<knowngammas.size();i++){
//    TLine *line = new TLine(knowngammas.at(i), 0.0, knowngammas.at(i), ymax);
//    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
//    line->Draw();
//  }
}

void Figure_ELabThetaLabAll(){

  TCanvas* cELabTLab = new TCanvas("cELabTLab","cELabTLab",1000,1000);

  chain->Draw("ELab:ThetaLab>>kEl(360,0,180,500,0,10)","abs(T_MUGAST_VAMOS-2750)<200 && MUST2.TelescopeNumber==5","colz");
  TH2F* kEl = (TH2F*) gDirectory->Get("kEl");
  kEl->SetTitle("");
  kEl->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  kEl->GetYaxis()->SetTitle("E_{lab} [MeV]");


  chain->Draw("ELab:ThetaLab>>kdp(360,0,180,500,0,10)","abs(T_MUGAST_VAMOS-2700)<400 && Mugast.TelescopeNumber>0","colz same");
  TH2F* kdp = (TH2F*) gDirectory->Get("kdp");

  chain->Draw("ELab:ThetaLab>>kdt(360,0,180,500,0,10)","abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber<5 && cutTime && cutTritons","colz same");
  TH2F* kdt = (TH2F*) gDirectory->Get("kdt");

  kEl->Draw("");
  kdp->Draw("same");
  kdt->Draw("same");

}

void CS_Diagnosis(){
  auto majorCanv = new TCanvas("CompareCanv","CompareCanv",1500,1500);
  majorCanv->Divide(3,3);
  canclone(majorCanv, 1, "c_peakFits_110_112"); 
  canclone(majorCanv, 2, "c_peakFits_115_118"); 
  canclone(majorCanv, 3, "c_peakFits_120_122"); 
  canclone(majorCanv, 4, "c_peakFits_125_128"); 
  canclone(majorCanv, 5, "c_peakFits_130_132"); 
  canclone(majorCanv, 6, "c_peakFits_135_138"); 
  canclone(majorCanv, 7, "c_peakFits_140_142"); 
  canclone(majorCanv, 8, "c_peakFits_145_148"); 
  canclone(majorCanv, 9, "c_peakFits_150_152"); 
}

void CS(){
/* Overload function */
  cout << "- CS(stateE, stateSp, orb_l, orb_j, nodes) "<< endl;
  cout << "---- 0.000, f7/2 = CS(0.000, 4, 3, 3.5, 0.10, \"\") "<< endl;
  cout << "---- 0.579, f7/2 = CS(0.579, 4, 3, 3.5, 0.10, \"\") "<< endl;
  cout << "---- 0.691, f7/2 = CS(0.691, 4, 3, 3.5, 0.10, \"\") "<< endl;
  cout << "---- 1.945, s1/2 = CS(1.945, 1, 0, 0.5, 0.10, \"\") "<< endl;
  cout << "---- 2.233, s1/2 = CS(2.233, 1, 0, 0.5, 0.10, \"\") "<< endl;
  cout << "---- 2.732, dX/2 = CS(2.732, 2, 2, 1.5, 0.10, \"\") "<< endl;
  cout << "---- 3.344 (x-?) "                                << endl;
  cout << "    FIT TOGETHER = CS(3.377, 1, 0, 0.5, 0.10, \"\")" << endl;
  cout << "---- 3.410 (?-?) "                                << endl;
  cout << "---- 4.297, d5/2 = CS(4.297, 3, 2, 2.5, 0.10, \"\") "<< endl;
  cout << "NOT SURE OF THE L TRANSFERS OR SPINS OF THESE STATES!!! CORRECT AS I GO ALONG!!!!" << endl;
}

/* MAIN FUNCTION */

void MM_Timing_Comparison(){
  TCanvas* c = new TCanvas("cMMT","cMMT",1000,1000);


  chain->Draw("MUST2.Si_T>>h1(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==1 && MUST2.CsI_E>0",
	      "same");
               TH1F* h1 = (TH1F*) gDirectory->Get("h1");
  chain->Draw("MUST2.Si_T>>h2(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==2 && MUST2.CsI_E>0",
	      "same");
               TH1F* h2 = (TH1F*) gDirectory->Get("h2");
  chain->Draw("MUST2.Si_T>>h3(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==3 && MUST2.CsI_E>0",
	      "same");
               TH1F* h3 = (TH1F*) gDirectory->Get("h3");
  chain->Draw("MUST2.Si_T>>h4(200,600,700)",
	      "abs(T_MUGAST_VAMOS-2750)<350 && MUST2.TelescopeNumber==4 && MUST2.CsI_E>0",
	      "same");
               TH1F* h4 = (TH1F*) gDirectory->Get("h4");

  h1->SetLineColor(kRed);
  h2->SetLineColor(kGreen);
  h3->SetLineColor(kBlue);
  h4->SetLineColor(kViolet);

//  h1->Draw();
//  h2->Draw("same");
//  h3->Draw("same");
//  h4->Draw("same");
}

void Plots_47Kdt(){

  /* Load graphical cut */
  //TFile gcIn("GraphicalCut_22Jun22.root");
  //TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritons");
  
  //TFile gcIn("cutTritonsWide.root");
  TFile gcIn("cutTritons_26Aug22Long.root");
  //TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritonsWide");
  TCutG* cutTritons = (TCutG*) gcIn.FindObjectAny("cutTritons");
  cutTritons->SetName("cutTritons");

  //TFile gcIn2("cutTime.root");
  TFile gcIn2("cutTimeNew.root");
  TCutG* cutTime = (TCutG*) gcIn2.FindObjectAny("cutTime");
  cutTime->SetName("cutTime");

  TFile gcIn3("cutDoublePeakGarbage.root");
  TCutG* cutGarbage = (TCutG*) gcIn3.FindObjectAny("cutDoublePeakGarbage");
  cutGarbage->SetName("cutGarbage");


  /**************/
  TFile gcInA("cutTritons_HighTLowE.root");
  TCutG* cutHighTLowE = (TCutG*) gcInA.FindObjectAny("cutTritons");
  cutHighTLowE->SetName("cutHighTLowE");

  TFile gcInB("cutTritons_SlimGate.root");
  TCutG* cutSlim = (TCutG*) gcInB.FindObjectAny("cutTritons");
  cutSlim->SetName("cutSlim");

  TFile gcInC("cutTritons_HighELowT.root");
  TCutG* cutHighELowT = (TCutG*) gcInC.FindObjectAny("cutTritons");
  cutHighELowT->SetName("cutHighELowT");


  /**************/

  LoadChain47Kdt();
  gStyle->SetOptStat("nemMrRi");

  tCentre = 280;  tRange = 30; //post-calibration
  //tCentre = 2750;  tRange = 350; //Wide is fine because I use the 2D time gate
  //tCentre = 2550;  tRange = 150;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "MUST2.TelescopeNumber>0 && MUST2.TelescopeNumber<5";
  reactionName = "47K(d,t)";
  
  cout << "==============================================" << endl;
  cout << "=============== (d,t) reaction ===============" << endl;
  cout << "==============================================" << endl;
  cout << " (d,t) selection: 'cutTritons' && 'cutTime'   " << endl;
  cout << "==============================================" << endl;
  cout << ""<< endl;
  CS();
  cout << ""<< endl;
  cout << "==============================================" << endl;

}
