void DrawELab10Be() {
  TCanvas* c1 = new TCanvas;
  c1->cd();
  TFile* f3 = new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/test.root");
  // new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bepalpha.root");
  TTree* t3 = (TTree*)f3->FindObjectAny("PhysicsTree");
  t3->Draw("ELab:ThetaLab>>hET10Bepalpha(1000,0,50,1000,0,500)", "CsI_E>0", "");
  TH2F* hBe1 = (TH2F*)gDirectory->FindObject("hET10Bepalpha");
  hBe1->SetMarkerStyle(1);
  hBe1->SetMarkerColor(kRed);
  hBe1->Draw();
  NPL::Reaction r3("10Be(p,4He)7Li@360MeV");
  r3.GetKinematicLine3()->Draw("same");
  r3.GetKinematicLine4()->Draw("same");

  // TFile* fa = new TFile("ha.root", "open");
  // TH2F* ha = (TH2F*)fa->FindObjectAny("h");
  // ha->SetMarkerColor(kBlack);
  // TFile* fli = new TFile("hli.root", "open");
  // TH2F* hli = (TH2F*)fli->FindObjectAny("h2");
  // hli->SetMarkerColor(kBlue);

  // ha->Draw("same");
  // hli->Draw("same");

  /*
  TCanvas *c2 = new TCanvas;
  c2->cd();
  TFile* f1 = new TFile("/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Analysis/AnaMugastAtLise10Bed6Li.root");
  TTree* t1 = (TTree*)f1->FindObjectAny("PhysicsTree");
  t1->Draw("ELab:ThetaLab>>hET10Bed6Li(1000,0,50,1000,0,500)", "", "");
  TH2F* hBe = (TH2F*) gDirectory->FindObject("hET10Bed6Li");
  hBe->SetMarkerStyle(2);
  hBe->SetMarkerColor(kBlue);
  hBe->Draw();
  NPL::Reaction r1("10Be(d,6Li)6He@300MeV");
  r1.GetKinematicLine3()->Draw("same");
  r1.GetKinematicLine4()->Draw("same");

  TFile * fa2 = new TFile("h6he.root","open");
  TH2F* ha2 = (TH2F*) fa2->FindObjectAny("h");
  TFile * fli2 = new TFile("h6li.root","open");
  TH2F* hli2 = (TH2F*) fli2->FindObjectAny("h2");

  ha2->Draw("same");
  hli2->Draw("same");
  */
}
