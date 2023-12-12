void show(){

  // TFile* file = new TFile("../../Outputs/Analysis/PhysicsTree.root","READ");
  TFile* file = new TFile("../../Outputs/Analysis/AnaMugastAtLise12Bepalpha.root","READ");
  TTree* tree = (TTree*) file->FindObjectAny("PhysicsTree");
  
  tree->Draw("ELab:ThetaLab>>h(1000,0,50,1000,0,500)","ELab>0","colz");
  NPL::Reaction r;
  r.ReadConfigurationFile("12Bepalpha.reaction");
  r.GetKinematicLine3()->Draw("c");
  r.GetKinematicLine4()->Draw("c");
  TCanvas*c =new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  tree->Draw("ELab:OriginalELab>>hE(1000,0,300,1000,0,300)","ELab>0","colz");
  TLine* lineE = new TLine(0,0,300,300);
  lineE->Draw();
  c->cd(2);
  tree->Draw("ThetaLab:OriginalThetaLab>>hT(1000,0,50,1000,0,50)","ELab>0","colz");
  TLine* lineT = new TLine(0,0,50,50);
  lineT->Draw();
  c->cd(3);
  tree->Draw("Y:X>>hY(300,-150,150,300,-150,150)","ELab>0","colz");
  c->cd(4);
   tree->Draw("Ex>>hEx(300,-10,10)","ELab>0","colz");        

}
