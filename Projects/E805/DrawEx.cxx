#include "DrawEx.h"


bool TestEntry(){
  GATCONFMASTER  = **GATCONFMASTER_;
  if(!(GATCONFMASTER.size() > 0 && GATCONFMASTER[0] == 1))
    return false;
  TACPhysics    = **TACPhysics_ ;
  ZDDPhysics = **ZDDPhysics_;
  
  bool TACCondition1 = false;
  for(unsigned int i = 0; i < TACPhysics.TAC_Time.size(); i++){
    // if(TACPhysics.TAC_Name[i].compare("TAC_CATS_PL") == 0 && abs((int)TACPhysics.TAC_Time[i] - 28000) < 1500 && ZDDPhysics.PL_TS.size() > 0 && abs((int)(TACPhysics.TAC_TS[i] - ZDDPhysics.PL_TS[0]) -61) <4)
    if(TACPhysics.TAC_Name[i].compare("TAC_CATS_PL") == 0 && cut_Cr->IsInside(TACPhysics.TAC_Time[i],ZDDPhysics.ICSum) && ZDDPhysics.PL_TS.size() > 0)
      TACCondition1 = true;
  }
  if(!TACCondition1)
    return false;
  
//
  //if(!( abs(ZDDPhysics.ICSum - 6750) < 200)) // deuton et triton
  //  return false;
  
  return true;
}

void UnallocateVariables(){
  //for(int i = 0; i < 10; i++){
  //PlasticRaw[i] = (*PlasticRaw_   )[i];
  //PlasticRawTS[i] = (*PlasticRaw_TS_)[i];
  //}
  //TAC_CATS_HF    = **TAC_CATS_HF_;
  //TAC_CATS_HFTS = **TAC_CATS_HF_TS_;
  //TAC_CATS_EXOGAM = **TAC_CATS_EXOGAM_;
  //TAC_CATS_EXOGAM = **TAC_CATS_EXOGAM_TS_;
  //TAC_MMG_CATS2  = **TAC_MMG_CATS2_;
  //TAC_MMG_CATS2TS = **TAC_MMG_CATS2_TS_;
  //TAC_MMG_CATS1  = **TAC_MMG_CATS1_;
  //TAC_MMG_CATS1TS = **TAC_MMG_CATS1_TS_;
  //TAC_MMG_EXOGAM = **TAC_MMG_EXOGAM_;
  //TAC_MMG_EXOGAMTS = **TAC_MMG_EXOGAM_TS_;
  //TAC_CATS1_CATS2 = **TAC_CATS1_CATS2_;
  //TAC_CATS1_CATS2TS = **TAC_CATS1_CATS2_TS_;
  //TAC_D4_CATS1   = **TAC_D4_CATS1_;
  //TAC_D4_CATS1TS = **TAC_D4_CATS1_TS_;
  //TAC_PL_1       = **TAC_PL_1_;
  //TAC_PL_1TS    = **TAC_PL_1_TS_;
  //TAC_PL_2       = **TAC_PL_2_;
  //TAC_PL_2TS    = **TAC_PL_2_TS_;
  //TAC_PL_3       = **TAC_PL_3_;
  //TAC_PL_3TS    = **TAC_PL_3_TS_;
  //TAC_PL_4       = **TAC_PL_4_;
  //TAC_PL_4TS    = **TAC_PL_4_TS_;
  //TAC_PL_5       = **TAC_PL_5_;
  //TAC_PL_5TS    = **TAC_PL_5_TS_;
  M2_CsI_E_p     = **M2_CsI_E_p_ ;
  M2_CsI_E_d     = **M2_CsI_E_d_ ;
  M2_CsI_E_t     = **M2_CsI_E_t_ ;
  M2_CsI_E_a     = **M2_CsI_E_a_ ;
  M2_ELab       = **M2_ELab_ ;
  M2_ThetaLab   = **M2_ThetaLab_ ;
  M2_ThetaCM     = **M2_ThetaCM_ ;
  M2_X           = **M2_X_ ;
  M2_Y           = **M2_Y_ ;
  M2_Z           = **M2_Z_ ;
  // M2_dE          = **M2_dE_ ;
  Must2Physics   = **Must2Physics_ ;
  CATSPhysics    = **CATSPhysics_ ;
  ExogamPhysics   = **ExogamPhysics_ ;

}

void DrawEx(){

  // c->Add("./data/NPRoot/Analysis/npa362.root");
  //  c->Add("./ssd/NPA_33S.root");
  //  c->Add("./ssd/NPA_34S.root");
  //  c->Add("./ssd/NPA_359.root");
   //c->Add("./ssd/NPA_360.root");
   //c->Add("./ssd/NPA_361.root");
   //c->Add("./ssd/NPA_362.root");
   c->Add("./ssd/NPA_364.root");
  //  c->Add("./data/NPRoot/Analysis/NPA_34S.root");
  //  c->Add("./data/NPRoot/Analysis/NPA_359.root");
   //c->Add("./data/NPRoot/Analysis/NPA_360.root");
   //c->Add("./data/NPRoot/Analysis/NPA_361.root");
  //  c->Add("./data/NPRoot/Analysis/NPA_364.root");
    // c->Add("./data_local/NPA_364.root");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_32*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_33*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_34*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_359*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_360*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_361*");
  //c->Add("/gold/nestar/data/GANIL/e805/NPRoot/Analysis/MUST_CSI_CATS_362*");
  TreeReader = new TTreeReader(c);


  
  GATCONFMASTER_ = new TTreeReaderValue<std::vector<unsigned int>>(*TreeReader,"GATCONF");
  // M2_TelescopeM_ = new TTreeReaderValue<unsigned short>(*TreeReader,"M2_TelescopeM");
  M2_CsI_E_p_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_CsI_E_p");
  M2_CsI_E_d_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_CsI_E_d");
  M2_CsI_E_t_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_CsI_E_t");
  M2_CsI_E_a_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_CsI_E_a");
  M2_ELab_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_ELab");
  M2_ThetaLab_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_ThetaLab");
  M2_ThetaCM_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_ThetaCM");
  M2_X_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_X");
  M2_Y_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_Y");
  M2_Z_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_Z");
  // M2_dE_ = new TTreeReaderValue<std::vector<double>>(*TreeReader,"M2_dE");
  Must2Physics_ = new TTreeReaderValue<TMust2Physics>(*TreeReader,"MUST2");
  CATSPhysics_ = new TTreeReaderValue<TCATSPhysics>(*TreeReader,"CATS");
  ZDDPhysics_ = new TTreeReaderValue<TZDDPhysics>(*TreeReader,"ZDD");
  ExogamPhysics_ = new TTreeReaderValue<TExogamPhysics>(*TreeReader,"Exogam");
  TACPhysics_ = new TTreeReaderValue<TTACPhysics>(*TreeReader,"TAC");


  // Cr48_pd.SetExcitationHeavy(0.472);

  TH2F* hExd_EGam = new TH2F("hExd_EGam","hExd_EGam",1000,-20,20,10000,0,10000);
  TH1F* hExd_WithCATS = new TH1F("hExd_w","hExd_w",450,-5,25);
  TH1F* hExd_WithCATSb = new TH1F("hExd_wb","hExd_wb",450,-5,25);
  TH1F* hExd_WithCATSs = new TH1F("hExd_ws","hExd_ws",450,-5,25);
  TH1F* hExd_NoCATS = new TH1F("hExd_n","hExd_n",450,-5,25);
  TH2F* hELTLd_s = new TH2F("hELTL_s","hELTL_s",800,0,40,1000,0,200);
  TH2F* hExTLd = new TH2F("hExTLd","hExTLd",800,0,40,200,-5,15);
  TH1F* hExt = new TH1F("hExt","hExt",450,-5,25);
  TH1F* hExd_gatedgama = new TH1F("hExd_gatedgama","hExd_ws",450,-5,25);
  std::map<unsigned int,TH1F*> hExd_WithCATS_m;
  std::map<unsigned int,TH1F*> hExd_WithCATSb_m;
  std::map<unsigned int,TH1F*> hExd_WithCATSs_m;
  std::map<unsigned int,std::map<unsigned int,TH1F*>*> hExd_WithCATSs_CSI;
  std::map<unsigned int,TH1F*> hExd_NoCATS_m;
  std::map<unsigned int,TH2F*> hELTLd_s_m;
  std::map<unsigned int,TH2F*> hExTLd_m;
  std::map<unsigned int,TH1F*> hExt_m;
  for(unsigned int i = 1; i <= 4; i++){
    hExd_WithCATS_m[i] = new TH1F(Form("hExd_w_%i",i),Form("hExd_w_%i",i),450,-5,25);
    hExd_WithCATSb_m[i] = new TH1F(Form("hExd_wb_%i",i),Form("hExd_wb_%i",i),450,-5,25);
    hExd_WithCATSs_m[i] = new TH1F(Form("hExd_ws_%i",i),Form("hExd_ws_%i",i),450,-5,25);
    hExd_NoCATS_m[i] = new TH1F(Form("hExd_n_%i",i),Form("hExd_n_%i",i),450,-5,25);
    hELTLd_s_m[i] = new TH2F(Form("hELTLd_s_%i",i),Form("hELTLd_s_%i",i),800,0,40,1000,0,200);
    hExTLd_m[i] = new TH2F(Form("hExTLd_%i",i),Form("hExTLd_%i",i),800,0,40,200,-5,15);
    hExt_m[i] = new TH1F(Form("hExt_%i",i),Form("hExt_%i",i),450,-5,25);
    hExd_WithCATSs_CSI[i] = new std::map<unsigned int, TH1F*>;
    for(unsigned int j = 1; j <= 16; j++){
      (*hExd_WithCATSs_CSI[i])[j] = new TH1F(Form("hExd_ws_T%i_CSI%i",i,j),Form("hExd_ws_T%i_CSI%i",i,j),300,-5,25);
    }
  }
  // TH1F* hExdNoEloss = new TH1F("hExdNoEloss","hExdNoEloss",450,-5,25);
  TH2F * hKine = new TH2F("hKine","hKine",1000,0,100,1000,0,100);
  
  c->LoadTree(0);
  ULong64_t GetEntries = c->GetEntries();
  if (!GetEntries){
      printf("ERROR in the Chain !! \n");
      return;
  }
  clock_t start = clock(), current, end;
  while(TreeReader->Next()){
        Long64_t cEntry = TreeReader->GetCurrentEntry();;
        if (cEntry%100000 == 0){
            current = clock();
            Double_t Frac = 1.0*cEntry/GetEntries;
            Double_t Time = ((double) (current-start)/CLOCKS_PER_SEC);
            Double_t TimeLeft = Time*(1/Frac - 1.);

            std::cout << "\rEntry : " << cEntry
                 << "/" << GetEntries
                 << " --- "
                 << Form("%.2f",100.*cEntry/GetEntries) <<" %"
                 << " --- "
                 <<Form("%.00f RunEvt/sec",cEntry/Time)
                <<" --- "
               << " Time Left : "<< Form("%d min ",(int)TimeLeft/60)
               << Form("%01.00d sec",(int)TimeLeft%60)
               << std::flush;
        }

    if(TestEntry()){
    UnallocateVariables();
    // std::cout << "TESTEntry" << std::endl;
    int mult = Must2Physics.Si_E.size();
    // std::cout << mult << std::endl;
    if (M2_CsI_E_d.size() == 1 &&  mult==1){
      double ELab = Must2Physics.Si_E[0];
      double ELab_t = Must2Physics.Si_E[0];
      double ELab_w, ELab_n, ELab_wb, ELab_ws;
      double TLab_w, TLab_n, TLab_wb, TLab_ws;
      double Ex_w, Ex_n, Ex_wb, Ex_ws;
      if(M2_CsI_E_d.at(0)>0)
        ELab+=M2_CsI_E_d.at(0);
      if(M2_CsI_E_t.at(0)>0)
        ELab_t+=M2_CsI_E_t.at(0);


        double t1;
        double PositionOnTargetX1;
        double PositionOnTargetY1;
        double t2;
        double PositionOnTargetX2;
        double PositionOnTargetY2;


      // Add energy Loss...
      if(CATSPhysics.PositionX.size() == 2 && CATSPhysics.PositionY.size() == 2){
        TVector3 InteractionPosition(M2_X.at(0), M2_Y.at(0), M2_Z.at(0)); 
        TVector3 BeamDirection(0,0,1);

        if(CATSPhysics.PositionX.size() == 2){
        t1 = (-23-1090.1)/(-23-1587.1);
        t2 = (0-1090.1)/(0-1587.1);
          if(CATSPhysics.DetNumber[0] == 1 && CATSPhysics.DetNumber[1] == 2)
        {
          double t = (0-(-1587.1))/(-1090.1-(-1587.1));
          // double t = (0-(-1587.1))/(1090.1);
          PositionOnTargetX1= CATSPhysics.PositionX[0]-2 + (CATSPhysics.PositionX[1] - CATSPhysics.PositionX[0]-2)*t;
          PositionOnTargetY1= CATSPhysics.PositionY[0] + (CATSPhysics.PositionY[1] - CATSPhysics.PositionY[0])*t;

        // std::cout << PositionOnTargetX2 << " " << CATSPhysics.PositionOnTargetX << std::endl;
        TVector3 BeamImpactCATS(CATSPhysics.PositionOnTargetX,CATSPhysics.PositionOnTargetY,0);
        // TVector3 BeamImpactCATS_s(CATSPhysics.PositionOnTargetX-2,CATSPhysics.PositionOnTargetY,0);
        // TVector3 BeamImpactCATS(PositionOnTargetX2,PositionOnTargetY2,0);
        TVector3 BeamImpactCATS_s(PositionOnTargetX1,PositionOnTargetY1,0);
        TVector3 BeamImpactCATS_NoCATS(0,0,0);
        TVector3 BeamDirection_CATS(CATSPhysics.PositionX[1] - CATSPhysics.PositionX[0],CATSPhysics.PositionY[1] - CATSPhysics.PositionY[0],CATSPhysics.PositionZ[1] - CATSPhysics.PositionZ[0]);
        //std::cout << BeamDirection_CATS.X() << " " << BeamDirection_CATS.Y() << " " << BeamDirection_CATS.Z() << std::endl;
        //TVector3 HitDirection_WithCATS = Must2Physics.GetPositionOfInteraction(0) - BeamImpactCATS;
        //TVector3 HitDirection_NoCATS = Must2Physics.GetPositionOfInteraction(0);
        TVector3 HitDirection_WithCATS = InteractionPosition - BeamImpactCATS;
        TVector3 HitDirection_WithCATS_s = InteractionPosition - BeamImpactCATS_s;
        TVector3 HitDirection_NoCATS = InteractionPosition - BeamImpactCATS_NoCATS;
        double ThetaNormalTarget_WithCATS = HitDirection_WithCATS.Angle(TVector3(0,0,1));
        double ThetaNormalTarget_WithCATS_s = HitDirection_WithCATS_s.Angle(TVector3(0,0,1));
        double ThetaNormalTarget_NoCATS = HitDirection_NoCATS.Angle(TVector3(0,0,1));
        // std::cout << ThetaNormalTarget_WithCATS << " " << ThetaNormalTarget_NoCATS << " " << CATSPhysics.PositionOnTargetX << " " << CATSPhysics.PositionOnTargetY<<  std::endl;
        if(cut_deuton->IsInside(M2_CsI_E_d.at(0), Must2Physics.Si_E[0]) && CATSPhysics.PositionX.size() == 2 && CATSPhysics.PositionY.size() == 2){
          // std::cout << "TEST" << std::endl;
          //hExdNoEloss->Fill(a.ReconstructRelativistic(ELab, M2_ThetaLab->at(0)*deg));
          // ELab = deuteron_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_WithCATS); 
          ELab_n = deuteron_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_NoCATS);
          ELab_w = deuteron_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_WithCATS); 
          ELab_ws =deuteron_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_WithCATS_s);
          TLab_n = HitDirection_NoCATS.Angle(BeamDirection);
          TLab_w = HitDirection_WithCATS.Angle(BeamDirection);
          TLab_ws =HitDirection_WithCATS_s.Angle(BeamDirection);
          Ex_n = Cr48_pd.ReconstructRelativistic( ELab_n,TLab_n);
          Ex_w = Cr48_pd.ReconstructRelativistic( ELab_w,TLab_w);
          Ex_ws = Cr48_pd.ReconstructRelativistic(ELab_ws,TLab_ws);
          Ex_wb = Cr48_pd.ReconstructRelativistic(deuteron_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_WithCATS), HitDirection_WithCATS.Angle(BeamDirection_CATS));
    
          TLorentzVector PHeavy_pd = Cr48_pd.LorentzAfterReaction(ELab_ws , TLab_ws);
          // TLorentzVector PHeavy_pt = Cr48_pt.LorentzAfterReaction(Energy["triton"] , M2_ThetaLab[countMust2]);
          // TLorentzVector PHeavy_p3He = Reaction_pd->LorentzAfterReaction(Energy["alpha"] , M2_ThetaLab[countMust2]);
          double Beta_pd = PHeavy_pd.Beta();
          // std::cout << Beta_pd << std::endl;
          // Beta_pt.push_back(PHeavy_pt.Beta());
          // Beta_p3He.push_back(PHeavy_p3He.Beta());
  
          int EXO_AB_size = ExogamPhysics.E_AB.size();
          for(unsigned int countExo = 0 ; countExo < EXO_AB_size; countExo++){
          // Doing Doppler correction only if one reaction occurs
              hExd_EGam->Fill(Ex_ws,Doppler_Correction(ExogamPhysics.Theta_D[countExo], ExogamPhysics.Phi_D[countExo], 0,0,Beta_pd,ExogamPhysics.E_AB[countExo]));
          }

          if(ExogamPhysics.E_AB.size()> 0)
            hExd_gatedgama->Fill(Ex_w);
          
          hExd_WithCATS->Fill(Ex_w);
          hExd_WithCATSb->Fill(Ex_wb);
          hExd_WithCATSs->Fill(Ex_ws);
          hExd_NoCATS->Fill(Ex_n);
          hELTLd_s->Fill(TLab_ws*180/3.14,ELab_ws);
          hExd_WithCATS_m[Must2Physics.TelescopeNumber[0]]->Fill(Ex_w);
          hExd_WithCATSb_m[Must2Physics.TelescopeNumber[0]]->Fill(Ex_wb);
          hExd_WithCATSs_m[Must2Physics.TelescopeNumber[0]]->Fill(Ex_ws);
          hExd_NoCATS_m[Must2Physics.TelescopeNumber[0]]->Fill(Ex_n);
          hELTLd_s_m[Must2Physics.TelescopeNumber[0]]->Fill(TLab_ws*180/3.14,ELab_ws);
          hExTLd->Fill(TLab_n*180/3.14,Ex_n);
          hExTLd_m[Must2Physics.TelescopeNumber[0]]->Fill(TLab_n*180/3.14,Ex_n);
          if(Must2Physics.CsI_N[0] > 0)
          (*hExd_WithCATSs_CSI[Must2Physics.TelescopeNumber[0]])[Must2Physics.CsI_N[0]]->Fill(Ex_ws);
          // hKine->Fill(M2_ThetaLab.,ELab);
        }
        if(cut_triton->IsInside(M2_CsI_E_t.at(0), Must2Physics.Si_E[0])){
          ELab_t = triton_CH2.EvaluateInitialEnergy(ELab_t, TargetThickness*0.5, ThetaNormalTarget_WithCATS); 
          ELab_t = Cr48_pt.ReconstructRelativistic(triton_CH2.EvaluateInitialEnergy(ELab, TargetThickness*0.5, ThetaNormalTarget_WithCATS), HitDirection_WithCATS_s.Angle(BeamDirection));
          hExt->Fill(ELab_t);
          hExt_m[Must2Physics.TelescopeNumber[0]]->Fill(ELab_t);
        }
        } 
        }
      }
    }
}
  }
  TFile * fout = new TFile("./data/NPRoot/Analysis/DrawEx10_3_shifted.root","RECREATE");
  hExd_gatedgama->Write("");
  hExd_EGam->Write("");
  hExd_WithCATS->Write("");
  hExd_WithCATSb->Write("");
  hExd_WithCATSs->Write("");
  hExd_NoCATS->Write("");
  hELTLd_s->Write("");
  hExTLd->Write("");
  hExt->Write("");
  for(unsigned int i = 1; i <= 4; i++){
    hExd_WithCATS_m[i]->Write("");
    hExd_WithCATSb_m[i]->Write("");
    hExd_WithCATSs_m[i]->Write("");
    hExd_NoCATS_m[i]->Write("");
    hELTLd_s_m[i]->Write("");
    hExTLd_m[i]->Write("");
    hExt_m[i]->Write("");
    for(unsigned int j = 1; j <= 16; j++){
      (*hExd_WithCATSs_CSI[i])[j]->Write("");
  }
  }
  // hKine->Write("");
  Cr48_pd.GetKinematicLine3()->Write("");
  fout->Close();
  TCanvas* c[4];
  TCanvas* csum[4];
  for(unsigned int i = 0; i <= 3; i++){
  c[i] = new TCanvas();
  csum[i] = new TCanvas();
  c[i]->Divide(4,4);
  csum[i]->Divide(2,2);
  for(unsigned int j = 1; j <= 4; j++){
    csum[i]->cd(j);
    if(j == 1){
    (*hExd_WithCATSs_CSI[i+1])[2]->Add((*hExd_WithCATSs_CSI[i+1])[5]);
    (*hExd_WithCATSs_CSI[i+1])[2]->Draw("");
    }
    else if(j == 2){
    (*hExd_WithCATSs_CSI[i+1])[3]->Add((*hExd_WithCATSs_CSI[i+1])[6]);
    (*hExd_WithCATSs_CSI[i+1])[3]->Add((*hExd_WithCATSs_CSI[i+1])[9]);
    (*hExd_WithCATSs_CSI[i+1])[3]->Draw("");
    }
    else if(j == 3){
    (*hExd_WithCATSs_CSI[i+1])[4]->Add((*hExd_WithCATSs_CSI[i+1])[7]);
    (*hExd_WithCATSs_CSI[i+1])[4]->Add((*hExd_WithCATSs_CSI[i+1])[10]);
    (*hExd_WithCATSs_CSI[i+1])[4]->Add((*hExd_WithCATSs_CSI[i+1])[13]);
    (*hExd_WithCATSs_CSI[i+1])[4]->Draw("");
    }
    else if(j == 4){
    (*hExd_WithCATSs_CSI[i+1])[8]->Add((*hExd_WithCATSs_CSI[i+1])[11]);
    (*hExd_WithCATSs_CSI[i+1])[8]->Add((*hExd_WithCATSs_CSI[i+1])[14]);
    (*hExd_WithCATSs_CSI[i+1])[8]->Draw("");
    }
  }
  for(unsigned int j = 1; j <= 16; j++){
    c[i]->cd(j);
    (*hExd_WithCATSs_CSI[i+1])[j]->Draw("");
  }
  c[i]->Draw("");
  csum[i]->Draw("");
  }

}