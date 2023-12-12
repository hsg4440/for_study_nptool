void DrawExELabThetaLab() {

  TFile* f1 = new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bed6Li.root");
  TTree* t1 = (TTree*)f1->FindObjectAny("PhysicsTree");
  TFile* f2 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise12Bepalpha.root");
  TTree* t2 = (TTree*)f2->FindObjectAny("PhysicsTree");
  TFile* f3 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bepalpha.root");
  TTree* t3 = (TTree*)f3->FindObjectAny("PhysicsTree");

  TCanvas* c1 = new TCanvas;
  c1->Divide(2, 1);
  c1->cd(1);
  t1->Draw("Ex+0.177>>hEx10Bed6Li(200,-20,20)", "CsI_E<150");
  c1->cd(2);
  t1->Draw("ELab:ThetaLab>>hET10Bed6Li(1000,0,50,1000,0,500)", "", "col");
  NPL::Reaction r1("10Be(d,6Li)6He@300MeV");
  r1.GetKinematicLine3()->Draw("same");
  r1.GetKinematicLine4()->Draw("same");

  TCanvas* c2 = new TCanvas;
  c2->Divide(2, 1);
  c2->cd(1);
  t2->Draw("Ex+0.177>>hEx12Bepalpha(200,-10,10)", "CsI_E>0&&CsI_E<150");
  c2->cd(2);
  t2->Draw("ELab:ThetaLab>>hET12Bepalpha(1000,0,50,1000,0,500)", "CsI_E>0", "col");
  NPL::Reaction r2("12Be(p,4He)9Li@360MeV");
  r2.GetKinematicLine3()->Draw("same");
  r2.GetKinematicLine4()->Draw("same");

  TCanvas* c3 = new TCanvas;
  c3->Divide(2, 1);
  c3->cd(1);
  t3->Draw("Ex+0.177>>hEx10Bepalpha(200,-10,10)", "CsI_E>0&&CsI_E<150");
  c3->cd(2);
  t3->Draw("ELab:ThetaLab>>hET10Bepalpha(1000,0,50,1000,0,500)", "CsI_E>0", "col");
  NPL::Reaction r3("10Be(p,4He)7Li@360MeV");
  r3.GetKinematicLine3()->Draw("same");
  r3.GetKinematicLine4()->Draw("same");
}

void Draw2() {
  TFile* f1 = new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bed6Li.root");
  TTree* t1 = (TTree*)f1->FindObjectAny("PhysicsTree");
  TFile* f2 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise12Bepalpha.root");
  TTree* t2 = (TTree*)f2->FindObjectAny("PhysicsTree");
  TFile* f3 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bepalpha.root");
  TTree* t3 = (TTree*)f3->FindObjectAny("PhysicsTree");
  TGraph* g1 = new TGraph("cross_sections/12CpAlphaNptoolCut20_50.txt");
  TGraph* g2 = new TGraph("cross_sections/12Cd6LiNptool03_70.txt");

  TCanvas* c4 = new TCanvas;
  c4->Divide(4, 1);

  c4->cd(1);
  g1->Draw("");
  g2->Draw("same");

  c4->cd(2);
  t1->Draw("ELab:ThetaLab>>hET10Bed6Li(1000,0,50,1000,0,500)", "", "col");
  t1->Draw("ELab:ThetaLab>>hET10Bed6Li(1000,0,50,1000,0,500)", "", "col");
  NPL::Reaction r1("10Be(d,6Li)6He@300MeV");
  r1.GetKinematicLine3()->Draw("same");
  r1.GetKinematicLine4()->Draw("same");
  c4->cd(3);
  t2->Draw("ELab:ThetaLab>>hET12Bepalpha(1000,0,50,1000,0,500)", "CsI_E>0", "col");
  NPL::Reaction r2("12Be(p,4He)9Li@456MeV");
  r2.GetKinematicLine3()->Draw("same");
  r2.GetKinematicLine4()->Draw("same");
  c4->cd(4);
  t3->Draw("ELab:ThetaLab>>hET10Bepalpha(1000,0,50,1000,0,500)", "CsI_E>0", "col");
  NPL::Reaction r3("10Be(p,4He)7Li@380MeV");
  r3.GetKinematicLine3()->Draw("same");
  r3.GetKinematicLine4()->Draw("same");
}

void DrawEff() {

  TFile* f1 = new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bed6Li.root");
  TTree* t1 = (TTree*)f1->FindObjectAny("PhysicsTree");
  TFile* f2 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise12Bepalpha.root");
  TTree* t2 = (TTree*)f2->FindObjectAny("PhysicsTree");
  TFile* f3 =
      new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bepalpha.root");
  TTree* t3 = (TTree*)f3->FindObjectAny("PhysicsTree");
  TGraph* g1 = new TGraph("cross_sections/12CpAlphaNptoolCut20_50.txt");
  TGraph* g2 = new TGraph("cross_sections/12Cd6LiNptool03_70.txt");

  // TFile* fEff = new TFile("hEffCM.root");
  // TH1F* heff = (TH1F*) fEff->FindObjectAny("hDetecThetaCM");
  TFile* fEff = new TFile("hEffMUGASTLISE.root");
  TH1F* heff = (TH1F*)fEff->FindObjectAny("hEfficiency");

  t3->Draw("ThetaCM>>h3(180,0,180)", "CsI_E<150", "");
  TH1F* h3 = (TH1F*)gDirectory->FindObjectAny("h3");

  h3->Divide(heff);
  h3->Scale(1. / 10.);
  h3->Scale(2);
  h3->Scale(2);
  g1->Draw("same");
  // TCanvas* c2 = new TCanvas();
}
