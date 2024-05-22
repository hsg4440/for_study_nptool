void compare2(){

 auto fz = new TFile("root/zaihong/run0582_RIG20210424_6He.root");
 auto tz = (TTree*) fz->FindObjectAny("rig");
 auto fl = new TFile("root/analysis/result582.root");
 auto tl = (TTree*) fl->FindObjectAny("PhysicsTree");
 tl->AddFriend(tz);
 auto cfdc0= new TCanvas();
 cfdc0->Divide(2,1);
 cfdc0->cd(1);
 tl->Draw("FDC0_X>>h1(1000,-100,100)","");
 tl->Draw("SamuraiFDC0.PosX>>hh1(1000,-100,100)","","same");
 cfdc0->cd(2);
 tl->Draw("FDC0_Y>>h2(1000,-100,100)","");
 tl->Draw("SamuraiFDC0.PosY>>hh2(1000,-100,100)","Minos.Tracks_P0@.size()>1","same");
  auto h1=(TH1*)gDirectory->FindObjectAny("h1");
  h1->SetLineColor(kOrange+7);
  auto hh1=(TH1*)gDirectory->FindObjectAny("hh1");
  hh1->SetLineColor(kAzure+7);
  auto h2=(TH1*)gDirectory->FindObjectAny("h2");
  h2->SetLineColor(kOrange+7);
  auto hh2=(TH1*)gDirectory->FindObjectAny("hh2");
  hh2->SetLineColor(kAzure+7);
  

 auto cfdc2= new TCanvas();
 cfdc2->Divide(2,1);
 cfdc2->cd(1);
 tl->Draw("FDC2_X>>h3(1000,-2000,2000)","");
 tl->Draw("SamuraiFDC2.PosX+252>>hh3(1000,-2000,2000)","Minos.Tracks_P0@.size()>1","same");
 cfdc2->cd(2);
 tl->Draw("FDC2_Y>>h4(1000,-2000,2000)","");
 tl->Draw("SamuraiFDC2.PosY>>hh4(1000,-2000,2000)","Minos.Tracks_P0@.size()>1","same");

  auto h3=(TH1*)gDirectory->FindObjectAny("h3");
  h3->SetLineColor(kOrange+7);
  auto hh3=(TH1*)gDirectory->FindObjectAny("hh3");
  hh3->SetLineColor(kAzure+7);
  hh3->Scale(h3->Integral()/hh3->Integral());
  auto h4=(TH1*)gDirectory->FindObjectAny("h4");
  h4->SetLineColor(kOrange+7);
  auto hh4=(TH1*)gDirectory->FindObjectAny("hh4");
  hh4->SetLineColor(kAzure+7);
  hh4->Scale(h4->Integral()/hh4->Integral());



}
