void windows(){

  auto fileM= new TFile("root/analysis/test824.root");
  auto treeM= (TTree*) fileM->FindObjectAny("PhysicsTree");
  auto fileR= new TFile("root/mrdc/run824/run824_MINOS.root");
  auto treeR= (TTree*) fileR->FindObjectAny("RawTree");

  treeM->AddFriend(treeR);
  auto c = new TCanvas();
  c->Divide(2,1);
  c->cd(1);
  treeM->Draw("Y_Vertex:X_Vertex>>h1(100,-100,100,100,-100,100)","","colx");
  c->cd(2);
  treeM->Draw("Y_Vertex:X_Vertex>>h2(100,-100,100,100,-100,100)","(Trigger&0x00ff)<64","col");


}
