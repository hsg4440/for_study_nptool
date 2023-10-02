
//////////////////////////////////////////
void Draw_cut()
{

  TChain* chain = new TChain("PhysicsTree");
  chain->Add("../../../root/analysis/run_55.root");
  //TFile* histofile = new TFile("histo/histo_EQ1_lg.root");
  TFile* cutfile = new TFile("cut10Be.root","update");

  TCanvas* c1 = new TCanvas("c1","c1",1200,1200);


  for(int det_number=1; det_number<9; det_number++){
    TString to_draw = Form("DeltaEcorr:Elab>>h%i(1000,0,200,1000,0,60)",det_number);
    TString condition = Form("fTS_TMW>0 && Telescope==%i",det_number);

    chain->Draw(to_draw,condition,"colz",5e7);

    TH2F *h = (TH2F*)gDirectory->FindObjectAny(Form("h%i",det_number));
    h->Draw("colz");

    TCutG* cutg;
    cutg = (TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    cutg->SetName(Form("cut10Be_det%i",det_number));

    cutfile->cd();
    cutg->Write();
  }
  //histofile->Close();
  cutfile->Close();

}

/////////////////////////////////////////
