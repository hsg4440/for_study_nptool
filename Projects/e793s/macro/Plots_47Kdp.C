string reactionName; /* defined by choice of dp or dt */
#include "DefineColours.h"
#include "GausFit.h"
#include "KnownPeakFitter.h"
#include "DrawPlots.h"

#include "CS2_master.h"
//#include "CS2.h"
//#include "CS2_MGX.h"
#include "ThreeBodyBreakup.h"
#include "ThreeBodyBreakup_FitPhaseSpace.h"
//#include "20Oct22_CompareYield.h"


void AddGammaLines(TH1F* hist, double particle, double ymax){
	
  string base = "sub ";

  for(int i=0; i<means.size();i++){
    string name = base + to_string(means.at(i));
    TLine *line = new TLine(particle-means.at(i), 0.0, particle-means.at(i), ymax);
    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
    line->Draw();
    TText *text = new TText((1.-(means.at(i)/particle))*particle,0.8*ymax,name.c_str());
    text->SetTextAngle(90);
    //text->SetTextSize(40);
    text->Draw();
  }
}

void AddPlacedGammas(TH1F* hist, double ymax){
  hist->Draw();
  for(int i=0; i<knowngammas.size();i++){
    TLine *line = new TLine(knowngammas.at(i), 0.0, knowngammas.at(i), ymax);
    line->SetLineColor(kBlack); line->SetLineStyle(kDotted);
    line->Draw();
  }
}

void CompareCountsInThetaLab(){
  auto canv = new TCanvas("cCompareExpSim","cCompareExpSim",1500,1500);
  auto file = new TFile("../../../Outputs/Analysis/Sim_47Kdp_10Aug22_TrueStripRemoval.root");  
  auto tree = (TTree*) file->FindObjectAny("PhysicsTree");

  canv->Divide(2,3);
  canv->cd(1);
    chain->Draw("ThetaLab>>exp1(120,100,160)","Mugast.TelescopeNumber==1","");
    tree->Draw("ThetaLab>>sim1(120,100,160)","Mugast.TelescopeNumber==1","same hist");
    auto exp1 = (TH1F*) gDirectory->Get("exp1");
    auto sim1 = (TH1F*) gDirectory->Get("sim1");
    exp1->SetLineColor(kRed);
    exp1->GetYaxis()->SetRangeUser(0.,700.);
    sim1->SetLineColor(kBlue);
    sim1->Scale(0.05);
  canv->cd(2);
    chain->Draw("ThetaLab>>exp2(120,100,160)","Mugast.TelescopeNumber==2","");
    tree->Draw("ThetaLab>>sim2(120,100,160)","Mugast.TelescopeNumber==2","same hist");
    auto exp2 = (TH1F*) gDirectory->Get("exp2");
    auto sim2 = (TH1F*) gDirectory->Get("sim2");
    exp2->SetLineColor(kRed);
    exp2->GetYaxis()->SetRangeUser(0.,700.);
    sim2->SetLineColor(kBlue);
    sim2->Scale(0.05);
  canv->cd(3);
    chain->Draw("ThetaLab>>exp3(120,100,160)","Mugast.TelescopeNumber==3","");
    tree->Draw("ThetaLab>>sim3(120,100,160)","Mugast.TelescopeNumber==3","same hist");
    auto exp3 = (TH1F*) gDirectory->Get("exp3");
    auto sim3 = (TH1F*) gDirectory->Get("sim3");
    exp3->SetLineColor(kRed);
    exp3->GetYaxis()->SetRangeUser(0.,700.);
    sim3->SetLineColor(kBlue);
    sim3->Scale(0.05);
  canv->cd(4);
    chain->Draw("ThetaLab>>exp4(120,100,160)","Mugast.TelescopeNumber==4","");
    tree->Draw("ThetaLab>>sim4(120,100,160)","Mugast.TelescopeNumber==4","same hist");
    auto exp4 = (TH1F*) gDirectory->Get("exp4");
    auto sim4 = (TH1F*) gDirectory->Get("sim4");
    exp4->SetLineColor(kRed);
    exp4->GetYaxis()->SetRangeUser(0.,700.);
    sim4->SetLineColor(kBlue);
    sim4->Scale(0.05);
  canv->cd(5);
    chain->Draw("ThetaLab>>exp5(120,100,160)","Mugast.TelescopeNumber==5","");
    tree->Draw("ThetaLab>>sim5(120,100,160)","Mugast.TelescopeNumber==5","same hist");
    auto exp5 = (TH1F*) gDirectory->Get("exp5");
    auto sim5 = (TH1F*) gDirectory->Get("sim5");
    exp5->SetLineColor(kRed);
    exp5->GetYaxis()->SetRangeUser(0.,700.);
    sim5->SetLineColor(kBlue);
    sim5->Scale(0.05);
  canv->cd(6);
    chain->Draw("ThetaLab>>exp7(120,100,160)","Mugast.TelescopeNumber==7","");
    tree->Draw("ThetaLab>>sim7(120,100,160)","Mugast.TelescopeNumber==7","same hist");
    auto exp7 = (TH1F*) gDirectory->Get("exp7");
    auto sim7 = (TH1F*) gDirectory->Get("sim7");
    exp7->SetLineColor(kRed);
    exp7->GetYaxis()->SetRangeUser(0.,700.);
    sim7->SetLineColor(kBlue);
    sim7->Scale(0.05);

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
// Overload function //
  cout << "- CS(stateE, stateSp, orb_l, orb_j, binEx, options) "      << endl;
  cout << "|-----------------| GOOD |-----------------|"     << endl;
  cout << "---- 0.143, p3/2 = CS(0.143, 2, 1, 1.5, 0.05, \"\") "   << endl;    
  cout << "---- 0.279, NOT POPULATED "                      << endl;
  cout << "---- 0.728, WEAKLY POPULATED "                   << endl;
  cout << "---- 0.967, p1/2 = CS(0.967, 0, 1, 0.5, 0.05, \"\") "   << endl;
  cout << "---- 1.409, p3/2 = CS(1.409, 1, 1, 1.5, 0.05, \"\") "   << endl;
  cout << "---- 1.978, p3/2 = CS(1.978, 1, 1, 0.5, 0.05, \"\") "   << endl;
  cout << "---- 2.407, p3/2 = CS(2.407, 0, 1, 0.5, 0.05, \"\") "   << endl;
  cout << "---- 2.908, p3/2 = CS(2.908, 2, 1, 1.5, 0.05, \"mixed\") "   << endl;
  cout << "---- 2.908, f5/2 = CS(2.908, 2, 3, 2.5, 0.05, \"mixed\") "   << endl;
  cout << "---- 3.254, f5/2 = CS(3.254, 3, 3, 2.5, 0.05, \"\") "   << endl;
  cout << "---- 3.601, f5/2 = CS(3.601, 3, 3, 2.5, 0.05, \"\") "   << endl;
  cout << endl;
  cout << "---- 3.792 (x-?) "                                << endl;
  cout << "       FIT TOGETHER = CS(3.830, 2, 3, 2.5, 0.05, \"\")" << endl;
  cout << "---- 3.868 (2-?) "                                << endl;
  cout << endl;
  cout << "|----------------| UNSURE |----------------|"     << endl;
  cout << "---- 4.061, f5/2 = CS(4.061, 3, 3, 2.5, 0.05, \"\") "   << endl;
  cout << "---- 4.387, f5/2 = CS(4.387, 3, 3, 2.5, 0.05, \"\") "   << endl;
}

void CompareGammas_ParticleRegions_48K(double binwidth, double minEx, double maxEx, double stepEx){
  
  vector<TH1F*>  Eg;
  vector<string> Names;
  int numHists = (int) ((maxEx-minEx)/stepEx);
  cout << numHists << endl;

  string draw = "AddBack_EDC>>EgTemp(" + to_string((int) (5./binwidth)) + ",0,5)";
  cout << draw << endl;

  auto cTEMP = new TCanvas("cTEMP","cTEMP",500,500);
  for(int i=0; i<numHists; i++){
    cout << "====================================" << endl;
    cout << i << endl;
    double centre = minEx + (stepEx*(double)i) + (0.5*stepEx);

    string gate = "Mugast.TelescopeNumber>0 && abs(T_MUGAST_VAMOS-2700)<400 && Ex@.size()==1 && abs(Ex-" 
  	      + to_string(centre) 
  	      + ")<" 
  	      + to_string(stepEx/2.0);
 
    string nameTemp = "Ex = " 
	            + to_string(centre-(0.5*stepEx)) 
		    + " to "
	            + to_string(centre+(0.5*stepEx));

    string histTemp = "Ex" 
	            + to_string(i);
//	            + to_string((int)centre-(0.5*stepEx)) 
//		    + "to"
//	            + to_string(centre+(0.5*stepEx));

    cout << gate << endl;

    chain->Draw(draw.c_str(),gate.c_str(),"");
    auto EgTemp = (TH1F*) gDirectory->Get("EgTemp");
    EgTemp->SetLineColor(i+1);
    EgTemp->SetNameTitle(histTemp.c_str(),histTemp.c_str());
    Eg.push_back(EgTemp);
    Names.push_back(nameTemp);
  
  }
  delete cTEMP;

  cout << "out" << endl;
 

  auto canv = new TCanvas("cCompareGammas","cCompareGammas",1000,1000);
  Eg[0]->Draw();
  for (int i=1; i<numHists; i++){
    Eg[i]->Draw("same");
  }

  cout << "out2" << endl;

  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  for (int i=0; i<numHists; i++){
    cout << "----------" << endl;
    cout << i << endl;
    legend->AddEntry(Eg[i],Names[i].c_str(),"f");
  }
  legend->Draw("same");

}


/* MAIN FUNCTION */

void Plots_47Kdp(){

cout << "test" << endl;

  LoadChain47Kdp();
  gStyle->SetOptStat("nemMrRi");

cout << "test" << endl;
  tCentre = 2700;  tRange = 400;
  timegate = "abs(T_MUGAST_VAMOS-" + to_string(tCentre) + ")<" + to_string(tRange);
  det_gate = "Mugast.TelescopeNumber>0 && Mugast.TelescopeNumber<8";
  reactionName = "47K(d,p)";

cout << "test" << endl;
  cout << "==============================================" << endl;
  cout << "=============== (d,p) reaction ===============" << endl;
  cout << "==============================================" << endl;
  cout << ""<< endl;
  CS();
  cout << ""<< endl;
  cout << "==============================================" << endl;

}
