void DrawExInvariant() {
  ////////////////////////////////////////////////////////////////////////////////
  TLorentzVector LV_alpha, LV_7Li;
  // TLorentzVector LV_alpha, LV_proton, LV_fragment;
  double mu = 931.5;                   // MeV/c2
  double m_alpha = 4 * mu + 2.425;     // Mev/c2
  double m_7Li = 7 * mu + 14.91;       // MeV/c2
  // double m_proton = mu + 7.289;        // Mev/c2
  // double m_fragment = 44 * mu - 3.755; // Mev/c2
  // double m_5Li = 5 * mu + 11.68;       // MeV/c2

  TChain* t = new TChain("PhysicsTree");
  t->Add("../../Outputs/Analysis/testnew.root");     // Name of simulated file

  // YOUR ALPHA AND Li CUT
  TFile* fcuta = new TFile("cuts/CUT_alpha.root", "read");
  TCutG* cuta = (TCutG*)gDirectory->FindObjectAny("CUT_alpha");
  TFile* fcutli = new TFile("cuts/CUT_7Li.root", "read");
  TCutG* cutli = (TCutG*)gDirectory->FindObjectAny("CUT_7Li");

  ////////////////////////////////////////////////////////////////////////////////
  TMust2Physics* M2 = new TMust2Physics();
  std::vector<double>* ELab = 0;
  std::vector<double>* ThetaLab = 0;
  std::vector<double>* X = 0;
  std::vector<double>* Y = 0;
  std::vector<double>* Z = 0;
  std::vector<double>* Ex = 0;
  TBranch* b_ELab;
  TBranch* b_ThetaLab;
  TBranch* b_X;
  TBranch* b_Y;
  TBranch* b_Z;
  TBranch* b_Ex;

  // float CATS1_X, CATS1_Y, CATS2_X, CATS2_Y, Xf, Yf;

  t->SetBranchStatus("*", false);
  t->SetBranchStatus("MUST2", true);
  t->SetBranchAddress("MUST2", &M2);
  t->SetBranchStatus("ELab", true);
  t->SetBranchAddress("ELab", &ELab, &b_ELab);
  t->SetBranchStatus("ThetaLab", true);
  t->SetBranchAddress("ThetaLab", &ThetaLab, &b_ThetaLab);
  // t->SetBranchStatus("CATS1_X", "true");
  // t->SetBranchAddress("CATS1_X", &CATS1_X);
  // t->SetBranchStatus("CATS1_Y", "true");
  // t->SetBranchAddress("CATS1_Y", &CATS1_Y);
  // t->SetBranchStatus("CATS2_X", "true");
  // t->SetBranchAddress("CATS2_X", &CATS2_X);
  // t->SetBranchStatus("CATS2_Y", "true");
  // t->SetBranchAddress("CATS2_Y", &CATS2_Y);
  t->SetBranchStatus("X", "true");
  t->SetBranchAddress("X", &X, &b_X);
  t->SetBranchStatus("Y", "true");
  t->SetBranchAddress("Y", &Y, &b_Y);
  t->SetBranchStatus("Z", "true");
  t->SetBranchAddress("Z", &Z, &b_Z);
  t->SetBranchStatus("Ex", "true");
  t->SetBranchAddress("Ex", &Ex, &b_Ex);

  TH1F* hExTot = new TH1F("hExTot", "hExTot", 1000, -100, 100);
  TH1F* hEx = new TH1F("hEx", "hEx", 1000, -100, 100);
  TH2F* hE1E2 = new TH2F("hE1E2", "hE1E2", 1000, 0, 500, 1000, 0, 500);

  ////////////////////////////////////////////////////////////////////////////////
  // int n_entries = t->GetEntries();
  int n_entries = 1e7;
  for (unsigned int i = 0; i < n_entries; i++) {
    t->GetEntry(i);

    double e_alpha;
    double theta_alpha;
    double phi_alpha;

    double e_li;
    double theta_li;
    double phi_li;

    double xtarget = 0;
    double ytarget = 0;
    // double CATSX_diff = 0;
    // double CATSY_diff = 0;

    // if (CATS1_X > -1500 && CATS2_X > -1500 && Xf > -1500) {
    //   CATSX_diff = -(CATS2_X - CATS1_X);
    //   xtarget = -Xf;
    // }
    // if (CATS1_Y > -1500 && CATS2_X > -1500 && Yf > -1500) {
    //   ytarget = Yf;
    //   CATSY_diff = CATS2_X - CATS1_X;
    // }

    TVector3 BeamImpact(xtarget, ytarget, 0);
    TVector3 BeamDirection(0,0,1);  

    if (M2->Si_E.size() == 2) {
      if ((cuta->IsInside(M2->CsI_E[0], M2->Si_E[0]) && cutli->IsInside(M2->CsI_E[1], M2->Si_E[1]))) {
        // Particle1 = alpha
        TVector3 PositionOfInteraction_alpha(X->at(0), Y->at(0), Z->at(0));
        TVector3 HitDirection_alpha = PositionOfInteraction_alpha - BeamImpact;
        e_alpha = ELab->at(0);
        theta_alpha = HitDirection_alpha.Angle(BeamDirection);
        phi_alpha = HitDirection_alpha.Phi();

        // Particle2 = 7li
        TVector3 PositionOfInteraction_li(X->at(1), Y->at(1), Z->at(1));
        TVector3 HitDirection_li = PositionOfInteraction_li - BeamImpact;
        e_li = ELab->at(1);
        theta_li = HitDirection_li.Angle(BeamDirection);
        phi_li = HitDirection_li.Phi();
      }
      else if (cuta->IsInside(M2->CsI_E[1], M2->Si_E[1]) && cutli->IsInside(M2->CsI_E[0], M2->Si_E[0])) {
        // Particle2 = alpha
        TVector3 PositionOfInteraction_alpha(X->at(1), Y->at(1), Z->at(1));
        TVector3 HitDirection_alpha = PositionOfInteraction_alpha - BeamImpact;
        e_alpha = ELab->at(1);
        theta_alpha = HitDirection_alpha.Angle(BeamDirection);
        phi_alpha = HitDirection_alpha.Phi();

        // Particle1 = 7li
        TVector3 PositionOfInteraction_li(X->at(0), Y->at(0), Z->at(0));
        TVector3 HitDirection_li = PositionOfInteraction_li - BeamImpact;
        e_li = ELab->at(0);
        theta_li = HitDirection_li.Angle(BeamDirection);
        phi_li = HitDirection_li.Phi();
      }

      ////////////////////////////////////////////////////////////////////////////////
      double gamma_alpha = e_alpha / m_alpha + 1;
      double beta_alpha = sqrt(1 - 1 / (pow(gamma_alpha, 2.)));
      double p_alpha_mag = gamma_alpha * m_alpha * beta_alpha;
      double etot_alpha = sqrt(pow(p_alpha_mag, 2) + pow(m_alpha, 2));

      TVector3 p_alpha(sin(theta_alpha) * cos(phi_alpha), // px
                       sin(theta_alpha) * sin(phi_alpha), // py
                       cos(theta_alpha));                 // pz
      p_alpha.SetMag(p_alpha_mag);
      LV_alpha.SetPxPyPzE(p_alpha.x(), p_alpha.y(), p_alpha.z(), etot_alpha);

      ////////////////////////////////////////////////////////////////////////////////
      double gamma_li = e_li / m_7Li + 1;
      double beta_li = sqrt(1 - 1 / (pow(gamma_li, 2.)));
      double p_li_mag = gamma_li * m_7Li * beta_li;
      double etot_li = sqrt(pow(p_li_mag, 2) + pow(m_7Li, 2));

      TVector3 p_li(sin(theta_li) * cos(phi_li), // px
                        sin(theta_li) * sin(phi_li), // py
                        cos(theta_li));                  // pz
      p_li.SetMag(p_li_mag);
      LV_7Li.SetPxPyPzE(p_li.x(), p_li.y(), p_li.z(), etot_li);
      TLorentzVector LV_total = LV_alpha + LV_7Li;

      double ExTot = LV_total.Mag() - m_alpha - m_7Li;
      hExTot->Fill(ExTot);
      if (Ex->size() >= 2) {  
        hEx->Fill(Ex->at(0));
        hEx->Fill(Ex->at(1));
      }
      hE1E2->Fill(e_alpha, e_li);
    }
  }
  TCanvas* c1 = new TCanvas();
  gStyle->SetOptStat(0);
  hExTot->Draw();
  hEx->SetLineColor(kRed);
  hEx->Draw("same");

  TF1* fitExTot = new TF1("fitExTot", "gaus", hExTot->GetXaxis()->GetXmin(), hExTot->GetXaxis()->GetXmax());
  hExTot->Fit(fitExTot, "Q");
  double meanExTot = fitExTot->GetParameter(1);
  std::cout << "Center of hExTot: " << meanExTot << " MeV" << std::endl; 

  double fitRangeMin = 0.0; 
  double fitRangeMax = 15.0; 
  TF1* fitEx = new TF1("fitEx", "gaus", fitRangeMin, fitRangeMax);
  hEx->Fit(fitEx);
  double meanEx = fitEx->GetParameter(1);
  std::cout << "Center of hEx: " << meanEx << " MeV" << std::endl;

  TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
  legend->AddEntry(hExTot, "hExTot", "l"); // "l" for line
  legend->AddEntry(hEx, "hex", "l");
  legend->Draw();
}

