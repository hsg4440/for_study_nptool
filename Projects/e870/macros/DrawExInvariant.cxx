void DrawDalitz() {
  ////////////////////////////////////////////////////////////////////////////////
  TLorentzVector LV_alpha, LV_proton, LV_fragment;
  double mu = 931.5;                   // MeV/c2
  double m_proton = mu + 7.289;        // Mev/c2
  double m_alpha = 4 * mu + 2.425;     // Mev/c2
  double m_fragment = 44 * mu - 3.755; // Mev/c2
  double m_5Li = 5 * mu + 11.68;       // MeV/c2

  TChain* t = new TChain("PhysicsTree"); // SimulatedTree
  t->Add("NPRootA/NPr0371_*a.root");     // Name of simulated file

  // YOUR ALPHA AND Li CUT
  TFile* fcuta = new TFile("CUT_e805/CUT_alpha_must.root", "read");
  TCutG* cuta = (TCutG*)gDirectory->FindObjectAny("CUT_alpha_must");
  TFile* fcutp = new TFile("CUT_e805/CUT_proton_must.root", "read");
  TCutG* cutp = (TCutG*)gDirectory->FindObjectAny("CUT_proton_must");

  ////////////////////////////////////////////////////////////////////////////////
  TMust2Physics* M2 = new TMust2Physics();
  std::vector<double>* ELab = 0;
  std::vector<double>* ThetaLab = 0;
  std::vector<double>* M2_X = 0;
  std::vector<double>* M2_Y = 0;
  std::vector<double>* M2_Z = 0;
  TBranch* b_ELab;
  TBranch* b_ThetaLab;
  TBranch* b_M2_X;
  TBranch* b_M2_Y;
  TBranch* b_M2_Z;

  float CATS1_X, CATS1_Y, CATS2_X, CATS2_Y, Xf, Yf;

  t->SetBranchStatus("*", false);
  t->SetBranchStatus("MUST2", true);
  t->SetBranchAddress("MUST2", &M2);
  t->SetBranchStatus("M2_ELab", true);
  t->SetBranchAddress("M2_ELab", &ELab, &b_ELab);
  t->SetBranchStatus("M2_ThetaLab", true);
  t->SetBranchAddress("M2_ThetaLab", &ThetaLab, &b_ThetaLab);
  t->SetBranchStatus("CATS1_X", "true");
  t->SetBranchAddress("CATS1_X", &CATS1_X);
  t->SetBranchStatus("CATS1_Y", "true");
  t->SetBranchAddress("CATS1_Y", &CATS1_Y);
  t->SetBranchStatus("CATS2_X", "true");
  t->SetBranchAddress("CATS2_X", &CATS2_X);
  t->SetBranchStatus("CATS2_Y", "true");
  t->SetBranchAddress("CATS2_Y", &CATS2_Y);
  t->SetBranchStatus("M2_X", "true");
  t->SetBranchAddress("M2_X", &M2_X, &b_M2_X);
  t->SetBranchStatus("M2_Y", "true");
  t->SetBranchAddress("M2_Y", &M2_Y, &b_M2_Y);
  t->SetBranchStatus("M2_Z", "true");
  t->SetBranchAddress("M2_Z", &M2_Z, &b_M2_Z);

  TH1F* hExTot = new TH1F("hExTot", "hExTot", 1000, -100, 100);
  TH1F* hEx5Li = new TH1F("hEx5Li", "hEx5Li", 1000, -100, 100);

  TH2F* hE1E2 = new TH2F("hE1E2", "hE1E2", 1000, 0, 500, 1000, 0, 500);

  ////////////////////////////////////////////////////////////////////////////////
  // int n_entries = t->GetEntries();
  int n_entries = 1e7;
  for (unsigned int i = 0; i < n_entries; i++) {
    t->GetEntry(i);

    double e_alpha;
    double theta_alpha;
    double phi_alpha;

    double e_proton;
    double theta_proton;
    double phi_proton;

    double e_fragment;
    double theta_fragment;
    double phi_fragment;

    double xtarget = 0;
    double ytarget = 0;
    double CATSX_diff = 0;
    double CATSY_diff = 0;

    if (CATS1_X > -1500 && CATS2_X > -1500 && Xf > -1500) {
      CATSX_diff = -(CATS2_X - CATS1_X);
      xtarget = -Xf;
    }
    if (CATS1_Y > -1500 && CATS2_X > -1500 && Yf > -1500) {
      ytarget = Yf;
      CATSY_diff = CATS2_X - CATS1_X;
    }

    TVector3 BeamImpact(xtarget, ytarget, 0);
    TVector3 BeamDirection(CATSX_diff, CATSY_diff, 497.1);

    if (M2->Si_E.size() == 2) {
      if ((cuta->IsInside(M2->CsI_E[0], M2->Si_E[0]) && cutp->IsInside(M2->CsI_E[1], M2->Si_E[1]))) {
        // Particle1 = alpha
        TVector3 PositionOfInteraction_alpha(M2_X->at(0), M2_Y->at(0), M2_Z->at(0));
        TVector3 HitDirection_alpha = PositionOfInteraction_alpha - BeamImpact;
        e_alpha = ELab->at(0);
        theta_alpha = HitDirection_alpha.Angle(BeamDirection);
        phi_alpha = HitDirection_alpha.Phi();

        // Particle2 = proton
        TVector3 PositionOfInteraction_proton(M2_X->at(1), M2_Y->at(1), M2_Z->at(1));
        TVector3 HitDirection_proton = PositionOfInteraction_proton - BeamImpact;
        e_proton = ELab->at(1);
        theta_proton = HitDirection_proton.Angle(BeamDirection);
        phi_proton = HitDirection_proton.Phi();
      }
      else if (cuta->IsInside(M2->CsI_E[1], M2->Si_E[1]) && cutp->IsInside(M2->CsI_E[0], M2->Si_E[0])) {
        // Particle2 = alpha
        TVector3 PositionOfInteraction_alpha(M2_X->at(1), M2_Y->at(1), M2_Z->at(1));
        TVector3 HitDirection_alpha = PositionOfInteraction_alpha - BeamImpact;
        e_alpha = ELab->at(1);
        theta_alpha = HitDirection_alpha.Angle(BeamDirection);
        phi_alpha = HitDirection_alpha.Phi();

        // Particle1 = proton
        TVector3 PositionOfInteraction_proton(M2_X->at(0), M2_Y->at(0), M2_Z->at(0));
        TVector3 HitDirection_proton = PositionOfInteraction_proton - BeamImpact;
        e_proton = ELab->at(0);
        theta_proton = HitDirection_proton.Angle(BeamDirection);
        phi_proton = HitDirection_proton.Phi();
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
      double gamma_proton = e_proton / m_proton + 1;
      double beta_proton = sqrt(1 - 1 / (pow(gamma_proton, 2.)));
      double p_proton_mag = gamma_proton * m_proton * beta_proton;
      double etot_proton = sqrt(pow(p_proton_mag, 2) + pow(m_proton, 2));

      TVector3 p_proton(sin(theta_proton) * cos(phi_proton), // px
                        sin(theta_proton) * sin(phi_proton), // py
                        cos(theta_proton));                  // pz
      p_proton.SetMag(p_proton_mag);
      LV_proton.SetPxPyPzE(p_proton.x(), p_proton.y(), p_proton.z(), etot_proton);
      TLorentzVector LV_total = LV_alpha + LV_proton;

      double ExTot = LV_total.Mag() - m_alpha - m_proton;
      hExTot->Fill(ExTot);
      if (Ex.size() >= 2) {
        hEx->Fill(Ex.at(0));
        hEx->Fill(Ex.at(1));
      }
      hE1E2->Fill(e_alpha, e_proton);
    }
  }
  TCanvas* c1 = new TCanvas();
  hExTot->Draw();
  hEx->Draw("same");
}
~
