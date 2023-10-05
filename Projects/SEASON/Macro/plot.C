#include "TGraph.h"
#include "TCanvas.h"
#include "TGraphErrors.h"

void PlotEfficiencyAlpha(string peaksfile){
  
  ifstream fin(peaksfile.c_str(), ios::in | ios::binary);
  
  string line;
  stringstream linestream;
  int TempDist;
  double TempEff_DSSD1, TempEff_Tunnel, temp;
  
  vector<double> Dist;
  vector<double> Eff_DSSD1, Eff_Tunnel, Eff_Tot;
  
  while(fin.good()){
    getline(fin,line);
    linestream.clear();
    linestream.str(line);
    if(linestream.fail() or line.empty()) continue;
    if(line.find("#") != string::npos) continue;
    
    linestream >> TempDist >> TempEff_DSSD1 >> temp >> TempEff_Tunnel >> temp;
    
    Dist.push_back(TempDist);
    Eff_DSSD1.push_back(TempEff_DSSD1*100.);
    Eff_Tunnel.push_back(TempEff_Tunnel*100.);
    Eff_Tot.push_back(TempEff_Tunnel*100.+TempEff_DSSD1*100.);
  }
  
  fin.close();
  
  
  TGraph *g1 = new TGraph(16,Dist.data(),Eff_DSSD1.data());

  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(50);
  g1->SetMarkerSize(4);
  g1->SetTitle(";DSSD1-Foil (mm);#epsilon_{#alpha} (%)");
  
  TCanvas *c = new TCanvas("EffDSSD1","EffDSSD1");
  g1->GetXaxis()->CenterTitle();
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->SetTitleOffset(0.75);
  g1->Draw("ap");
  c->Update();
  
  TGraph *g2 = new TGraph(16,Dist.data(),Eff_Tunnel.data());

  g2->SetMarkerColor(1);
  g2->SetMarkerStyle(50);
  g2->SetMarkerSize(4);
  g2->SetTitle(";Tunnel-Foil (mm);#epsilon_{#alpha} (%)");
  
  TCanvas *c2 = new TCanvas("EffTunnel","EffTunnel");
  g2->GetXaxis()->CenterTitle();
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->CenterTitle();
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->SetTitleOffset(0.75);
  g2->Draw("ap");
  c2->Update();
  
  TGraph *g3 = new TGraph(16,Dist.data(),Eff_Tot.data());

  g3->SetMarkerColor(1);
  g3->SetMarkerStyle(50);
  g3->SetMarkerSize(4);
  g3->SetTitle(";Detector-Foil (mm);#epsilon_{#alpha} (%)");
  
  TCanvas *c3 = new TCanvas("EffTot","EffTot");
  g3->GetXaxis()->CenterTitle();
  g3->GetXaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->CenterTitle();
  g3->GetYaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->SetTitleOffset(0.75);
  g3->Draw("ap");
  c3->Update();
}

void PlotEfficiencyElectrons(string peaksfile){
  
  ifstream fin(peaksfile.c_str(), ios::in | ios::binary);
  
  string line;
  stringstream linestream;
  int TempDist;
  double TempEff_DSSD1, TempEff_Tunnel;
  
  vector<double> Dist;
  vector<double> Eff_DSSD1, Eff_Tunnel, Eff_Tot;
  
  while(fin.good()){
    getline(fin,line);
    linestream.clear();
    linestream.str(line);
    if(linestream.fail() or line.empty()) continue;
    if(line.find("#") != string::npos) continue;
    
    linestream >> TempDist >> TempEff_DSSD1 >> TempEff_Tunnel;
    
    Dist.push_back(TempDist);
    Eff_DSSD1.push_back(TempEff_DSSD1*100.);
    Eff_Tunnel.push_back(TempEff_Tunnel*100.);
    Eff_Tot.push_back(TempEff_Tunnel*100.+TempEff_DSSD1*100.);
  }
  
  fin.close();
  
  
  TGraph *g1 = new TGraph(Dist.size(),Dist.data(),Eff_DSSD1.data());

  g1->SetMarkerColor(1);
  g1->SetMarkerStyle(50);
  g1->SetMarkerSize(4);
  g1->SetTitle(";DSSD1-Foil (mm);#epsilon_{e-} (%)");
  
  TCanvas *c = new TCanvas("EffDSSD1","EffDSSD1");
  g1->GetXaxis()->CenterTitle();
  g1->GetXaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitleSize(0.05);
  g1->GetYaxis()->SetTitleOffset(0.75);
  g1->Draw("ap");
  c->Update();
  
  TGraph *g2 = new TGraph(Dist.size(),Dist.data(),Eff_Tunnel.data());

  g2->SetMarkerColor(1);
  g2->SetMarkerStyle(50);
  g2->SetMarkerSize(4);
  g2->SetTitle(";Tunnel-Foil (mm);#epsilon_{e-} (%)");
  
  TCanvas *c2 = new TCanvas("EffTunnel","EffTunnel");
  g2->GetXaxis()->CenterTitle();
  g2->GetXaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->CenterTitle();
  g2->GetYaxis()->SetTitleSize(0.05);
  g2->GetYaxis()->SetTitleOffset(0.75);
  g2->Draw("ap");
  c2->Update();
  
  TGraph *g3 = new TGraph(Dist.size(),Dist.data(),Eff_Tot.data());

  g3->SetMarkerColor(1);
  g3->SetMarkerStyle(50);
  g3->SetMarkerSize(4);
  g3->SetTitle(";Detector-Foil (mm);#epsilon_{e-} (%)");
  
  TCanvas *c3 = new TCanvas("EffTot","EffTot");
  g3->GetXaxis()->CenterTitle();
  g3->GetXaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->CenterTitle();
  g3->GetYaxis()->SetTitleSize(0.05);
  g3->GetYaxis()->SetTitleOffset(0.75);
  g3->Draw("ap");
  c3->Update();
}

void PlotMultipleElectronsEfficiency(string peaksfile){
  
  ifstream fin(peaksfile.c_str(), ios::in | ios::binary);
  
  string line;
  stringstream linestream;
  int TempDist, TempNPart;
  double TempEffDSSD1, TempEffDSSD1Err, TempEffTunnel, TempEffTunnelErr, TempEffTot, TempEffTotErr;
  
  vector<int> Dist, NPart;
  vector<double> EffDSSD1, EffDSSD1Err, EffTunnel, EffTunnelErr, EffTot, EffTotErr;
  
  while(fin.good()){
    getline(fin,line);
    linestream.clear();
    linestream.str(line);
    if(linestream.fail() or line.empty()) continue;
    if(line.find("#") != string::npos) continue;
    
    linestream >> TempNPart >> TempDist >> TempEffDSSD1 >> TempEffDSSD1Err >> TempEffTunnel >> TempEffTunnelErr >> TempEffTot >> TempEffTotErr;
    
    NPart.push_back(TempNPart);
    Dist.push_back(TempDist);
    EffDSSD1.push_back(TempEffDSSD1*100.);
    EffDSSD1Err.push_back(TempEffDSSD1Err*100.);
    EffTunnel.push_back(TempEffTunnel*100.);
    EffTunnelErr.push_back(TempEffTunnelErr*100.);
    EffTot.push_back(TempEffTot*100.);
    EffTotErr.push_back(TempEffTotErr*100.);
  }
  
  fin.close();
  
  double Dist_1[7];
  
  double EffTot_1[7];
  double EffTot_2[7];
  double EffTot_3[7];  
  double EffTot_4[7];
  double EffTot_5[7];
  double EffTot_6[7];  
  
  double EffTotErr_1[7];
  double EffTotErr_2[7];
  double EffTotErr_3[7];  
  double EffTotErr_4[7];
  double EffTotErr_5[7];
  double EffTotErr_6[7];  
  
  double EffDSSD1_1[7];
  double EffDSSD1_2[7];
  double EffDSSD1_3[7];  
  double EffDSSD1_4[7];
  double EffDSSD1_5[7];
  double EffDSSD1_6[7];  
  
  double EffDSSD1Err_1[7];
  double EffDSSD1Err_2[7];
  double EffDSSD1Err_3[7];  
  double EffDSSD1Err_4[7];
  double EffDSSD1Err_5[7];
  double EffDSSD1Err_6[7];  
  
  double EffTunnel_1[7];
  double EffTunnel_2[7];
  double EffTunnel_3[7];  
  double EffTunnel_4[7];
  double EffTunnel_5[7];
  double EffTunnel_6[7];  
  
  double EffTunnelErr_1[7];
  double EffTunnelErr_2[7];
  double EffTunnelErr_3[7];  
  double EffTunnelErr_4[7];
  double EffTunnelErr_5[7];
  double EffTunnelErr_6[7];  
  
  for(int i=0;i<7;i++){
    Dist_1[i] = Dist[i];
    
    EffTot_1[i] = EffTot[i];
    EffTot_2[i] = EffTot[i+7];
    EffTot_3[i] = EffTot[i+14];
    EffTot_4[i] = EffTot[i+21];
    EffTot_5[i] = EffTot[i+28];
    EffTot_6[i] = EffTot[i+35];
    
    EffTotErr_1[i] = EffTotErr[i];
    EffTotErr_2[i] = EffTotErr[i+7];
    EffTotErr_3[i] = EffTotErr[i+14];
    EffTotErr_4[i] = EffTotErr[i+21];
    EffTotErr_5[i] = EffTotErr[i+28];
    EffTotErr_6[i] = EffTotErr[i+35];
    
    EffDSSD1_1[i] = EffDSSD1[i];
    EffDSSD1_2[i] = EffDSSD1[i+7];
    EffDSSD1_3[i] = EffDSSD1[i+14];
    EffDSSD1_4[i] = EffDSSD1[i+21];
    EffDSSD1_5[i] = EffDSSD1[i+28];
    EffDSSD1_6[i] = EffDSSD1[i+35];
    
    EffDSSD1Err_1[i] = EffDSSD1Err[i];
    EffDSSD1Err_2[i] = EffDSSD1Err[i+7];
    EffDSSD1Err_3[i] = EffDSSD1Err[i+14];
    EffDSSD1Err_4[i] = EffDSSD1Err[i+21];
    EffDSSD1Err_5[i] = EffDSSD1Err[i+28];
    EffDSSD1Err_6[i] = EffDSSD1Err[i+35];
    
    EffTunnel_1[i] = EffTunnel[i];
    EffTunnel_2[i] = EffTunnel[i+7];
    EffTunnel_3[i] = EffTunnel[i+14];
    EffTunnel_4[i] = EffTunnel[i+21];
    EffTunnel_5[i] = EffTunnel[i+28];
    EffTunnel_6[i] = EffTunnel[i+35];
    
    EffTunnelErr_1[i] = EffTunnelErr[i];
    EffTunnelErr_2[i] = EffTunnelErr[i+7];
    EffTunnelErr_3[i] = EffTunnelErr[i+14];
    EffTunnelErr_4[i] = EffTunnelErr[i+21];
    EffTunnelErr_5[i] = EffTunnelErr[i+28];
    EffTunnelErr_6[i] = EffTunnelErr[i+35];
  }
  
  
  TGraphErrors *gTot1 = new TGraphErrors(7,Dist_1,EffTot_1,NULL,EffTotErr_1);
  TGraphErrors *gTot2 = new TGraphErrors(7,Dist_1,EffTot_2,NULL,EffTotErr_2);
  TGraphErrors *gTot3 = new TGraphErrors(7,Dist_1,EffTot_3,NULL,EffTotErr_3);
  TGraphErrors *gTot4 = new TGraphErrors(7,Dist_1,EffTot_4,NULL,EffTotErr_4);
  TGraphErrors *gTot5 = new TGraphErrors(7,Dist_1,EffTot_5,NULL,EffTotErr_5);
  TGraphErrors *gTot6 = new TGraphErrors(7,Dist_1,EffTot_6,NULL,EffTotErr_6);
  
  gTot1->SetMarkerColor(1);
  gTot2->SetMarkerColor(2);
  gTot3->SetMarkerColor(3);
  gTot4->SetMarkerColor(4);
  gTot5->SetMarkerColor(6);
  gTot6->SetMarkerColor(7);
  
  gTot1->SetMarkerStyle(50);
  gTot2->SetMarkerStyle(52);
  gTot3->SetMarkerStyle(50);
  gTot4->SetMarkerStyle(52);
  gTot5->SetMarkerStyle(50);
  gTot6->SetMarkerStyle(52);
  
  gTot1->SetMarkerSize(4);
  gTot2->SetMarkerSize(4);
  gTot3->SetMarkerSize(4);
  gTot4->SetMarkerSize(4);
  gTot5->SetMarkerSize(4);
  gTot6->SetMarkerSize(4);
  
  
  TGraphErrors *gDSSD1_1 = new TGraphErrors(7,Dist_1,EffDSSD1_1,NULL,EffDSSD1Err_1);
  TGraphErrors *gDSSD1_2 = new TGraphErrors(7,Dist_1,EffDSSD1_2,NULL,EffDSSD1Err_2);
  TGraphErrors *gDSSD1_3 = new TGraphErrors(7,Dist_1,EffDSSD1_3,NULL,EffDSSD1Err_3);
  TGraphErrors *gDSSD1_4 = new TGraphErrors(7,Dist_1,EffDSSD1_4,NULL,EffDSSD1Err_4);
  TGraphErrors *gDSSD1_5 = new TGraphErrors(7,Dist_1,EffDSSD1_5,NULL,EffDSSD1Err_5);
  TGraphErrors *gDSSD1_6 = new TGraphErrors(7,Dist_1,EffDSSD1_6,NULL,EffDSSD1Err_6);
  
  gDSSD1_1->SetMarkerColor(1);
  gDSSD1_2->SetMarkerColor(2);
  gDSSD1_3->SetMarkerColor(3);
  gDSSD1_4->SetMarkerColor(4);
  gDSSD1_5->SetMarkerColor(6);
  gDSSD1_6->SetMarkerColor(7);
  
  gDSSD1_1->SetMarkerStyle(50);
  gDSSD1_2->SetMarkerStyle(52);
  gDSSD1_3->SetMarkerStyle(50);
  gDSSD1_4->SetMarkerStyle(52);
  gDSSD1_5->SetMarkerStyle(50);
  gDSSD1_6->SetMarkerStyle(52);
  
  gDSSD1_1->SetMarkerSize(4);
  gDSSD1_2->SetMarkerSize(4);
  gDSSD1_3->SetMarkerSize(4);
  gDSSD1_4->SetMarkerSize(4);
  gDSSD1_5->SetMarkerSize(4);
  gDSSD1_6->SetMarkerSize(4);
  
  TGraphErrors *gTunnel_1 = new TGraphErrors(7,Dist_1,EffTunnel_1,NULL,EffTunnelErr_1);
  TGraphErrors *gTunnel_2 = new TGraphErrors(7,Dist_1,EffTunnel_2,NULL,EffTunnelErr_2);
  TGraphErrors *gTunnel_3 = new TGraphErrors(7,Dist_1,EffTunnel_3,NULL,EffTunnelErr_3);
  TGraphErrors *gTunnel_4 = new TGraphErrors(7,Dist_1,EffTunnel_4,NULL,EffTunnelErr_4);
  TGraphErrors *gTunnel_5 = new TGraphErrors(7,Dist_1,EffTunnel_5,NULL,EffTunnelErr_5);
  TGraphErrors *gTunnel_6 = new TGraphErrors(7,Dist_1,EffTunnel_6,NULL,EffTunnelErr_6);
  
  gTunnel_1->SetMarkerColor(1);
  gTunnel_2->SetMarkerColor(2);
  gTunnel_3->SetMarkerColor(3);
  gTunnel_4->SetMarkerColor(4);
  gTunnel_5->SetMarkerColor(6);
  gTunnel_6->SetMarkerColor(7);
  
  gTunnel_1->SetMarkerStyle(50);
  gTunnel_2->SetMarkerStyle(52);
  gTunnel_3->SetMarkerStyle(50);
  gTunnel_4->SetMarkerStyle(52);
  gTunnel_5->SetMarkerStyle(50);
  gTunnel_6->SetMarkerStyle(52);
  
  gTunnel_1->SetMarkerSize(4);
  gTunnel_2->SetMarkerSize(4);
  gTunnel_3->SetMarkerSize(4);
  gTunnel_4->SetMarkerSize(4);
  gTunnel_5->SetMarkerSize(4);
  gTunnel_6->SetMarkerSize(4);
  
  
  TMultiGraph *mgTot = new TMultiGraph();
  mgTot->SetTitle(";Detector-Foil distance (mm);#epsilon_{e-} (%)");
  mgTot->Add(gTot1,"p");
  mgTot->Add(gTot2,"p");
  mgTot->Add(gTot3,"p");
  mgTot->Add(gTot4,"p");
  mgTot->Add(gTot5,"p");
  mgTot->Add(gTot6,"p");
  
  TMultiGraph *mgDSSD1 = new TMultiGraph();
  mgDSSD1->SetTitle(";Detector-Foil distance (mm);#epsilon_{e-} (%)");
  mgDSSD1->Add(gDSSD1_1,"p");
  mgDSSD1->Add(gDSSD1_2,"p");
  mgDSSD1->Add(gDSSD1_3,"p");
  mgDSSD1->Add(gDSSD1_4,"p");
  mgDSSD1->Add(gDSSD1_5,"p");
  mgDSSD1->Add(gDSSD1_6,"p");
  
  TMultiGraph *mgTunnel = new TMultiGraph();
  mgTunnel->SetTitle(";Detector-Foil distance (mm);#epsilon_{e-} (%)");
  mgTunnel->Add(gTunnel_1,"p");
  mgTunnel->Add(gTunnel_2,"p");
  mgTunnel->Add(gTunnel_3,"p");
  mgTunnel->Add(gTunnel_4,"p");
  mgTunnel->Add(gTunnel_5,"p");
  mgTunnel->Add(gTunnel_6,"p");
  
  
  
  auto legend1 = new TLegend(0.65,0.55,0.9,0.9);
  legend1->AddEntry(gTot1,"Emission of 1 electron","p");
  legend1->AddEntry(gTot2,"Emission of 1 #alpha and 1 electron","p");
  legend1->AddEntry(gTot3,"Emission of 1 #alpha and 2 electrons","p");
  legend1->AddEntry(gTot4,"Emission of 1 #alpha and 3 electrons","p");
  legend1->AddEntry(gTot5,"Emission of 1 #alpha and 4 electrons","p");
  legend1->AddEntry(gTot6,"Emission of 1 #alpha and 5 electrons","p");
  
  auto legend2 = new TLegend(0.65,0.55,0.9,0.9);
  legend2->AddEntry(gDSSD1_1,"Emission of 1 electron","p");
  legend2->AddEntry(gDSSD1_2,"Emission of 1 #alpha and 1 electron","p");
  legend2->AddEntry(gDSSD1_3,"Emission of 1 #alpha and 2 electrons","p");
  legend2->AddEntry(gDSSD1_4,"Emission of 1 #alpha and 3 electrons","p");
  legend2->AddEntry(gDSSD1_5,"Emission of 1 #alpha and 4 electrons","p");
  legend2->AddEntry(gDSSD1_6,"Emission of 1 #alpha and 5 electrons","p");
  
  auto legend3 = new TLegend(0.65,0.55,0.9,0.9);
  legend3->AddEntry(gTunnel_1,"Emission of 1 electron","p");
  legend3->AddEntry(gTunnel_2,"Emission of 1 #alpha and 1 electron","p");
  legend3->AddEntry(gTunnel_3,"Emission of 1 #alpha and 2 electrons","p");
  legend3->AddEntry(gTunnel_4,"Emission of 1 #alpha and 3 electrons","p");
  legend3->AddEntry(gTunnel_5,"Emission of 1 #alpha and 4 electrons","p");
  legend3->AddEntry(gTunnel_6,"Emission of 1 #alpha and 5 electrons","p");
  

  TCanvas *c1 = new TCanvas("EffTot","EffTot");
  mgTot->GetXaxis()->CenterTitle();
  mgTot->GetXaxis()->SetTitleSize(0.05);
  mgTot->GetYaxis()->CenterTitle();
  mgTot->GetYaxis()->SetTitleSize(0.05);
  mgTot->GetYaxis()->SetTitleOffset(0.75);
  mgTot->Draw("ap");
  legend1->Draw();
  c1->Update();
  
  TCanvas *c2 = new TCanvas("EffDSSD1","EffDSSD1");
  mgDSSD1->GetXaxis()->CenterTitle();
  mgDSSD1->GetXaxis()->SetTitleSize(0.05);
  mgDSSD1->GetYaxis()->CenterTitle();
  mgDSSD1->GetYaxis()->SetTitleSize(0.05);
  mgDSSD1->GetYaxis()->SetTitleOffset(0.75);
  mgDSSD1->Draw("ap");
  legend2->Draw();
  c2->Update();
  
  TCanvas *c3 = new TCanvas("EffTunnel","EffTunnel");
  mgTunnel->GetXaxis()->CenterTitle();
  mgTunnel->GetXaxis()->SetTitleSize(0.05);
  mgTunnel->GetYaxis()->CenterTitle();
  mgTunnel->GetYaxis()->SetTitleSize(0.05);
  mgTunnel->GetYaxis()->SetTitleOffset(0.75);
  mgTunnel->Draw("ap");
  legend3->Draw();
  c3->Update();

}
