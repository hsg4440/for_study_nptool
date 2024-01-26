void DrawExInvariant() {
  ////////////////////////////////////////////////////////////////////////////////
  TLorentzVector LV_alpha, LV_7Li, LV_p;
  double mu = 931.5;               // MeV/c2
  double m_alpha = 4 * mu + 2.425; // Mev/c2
  double m_7Li = 7 * mu + 14.91;   // MeV/c2
  double m_10Be = 10 * mu + 12.61; // MeV/c2
  double mp = mu + 7.289;          // MeV/c2
  // double m_5Li = 5 * mu + 11.68;       // MeV/c2

  TChain* t = new TChain("PhysicsTree");
  t->Add("../../Outputs/Analysis/Ana10Bepalpha_38AMeV_CH2.root"); // Name of simulated file

  // YOUR ALPHA AND Li CUT
  TFile* fcuta = new TFile("cuts/CUT_alpha.root", "read");
  TCutG* cuta = (TCutG*)gDirectory->FindObjectAny("CUT_alpha");
  TFile* fcutli = new TFile("cuts/CUT_7Li.root", "read");
  TCutG* cutli = (TCutG*)gDirectory->FindObjectAny("CUT_7Li");

  ////////////////////////////////////////////////////////////////////////////////
  TMust2Physics* M2 = new TMust2Physics();
  std::vector<double>* ELab = 0;
  std::vector<double>* ThetaLab = 0;
  std::vector<double>* PhiLab = 0;
  std::vector<double>* Ex = 0;
  std::vector<double>* ExNoCATS = 0;
  TBranch* b_ELab;
  TBranch* b_ThetaLab;
  TBranch* b_PhiLab;
  TBranch* b_X;
  TBranch* b_Y;
  TBranch* b_Z;
  TBranch* b_Ex;
  TBranch* b_ExNoCATS;

  t->SetBranchStatus("*", false);
  t->SetBranchStatus("MUST2", true);
  t->SetBranchAddress("MUST2", &M2);
  t->SetBranchStatus("ELab", true);
  t->SetBranchAddress("ELab", &ELab, &b_ELab);
  t->SetBranchStatus("ThetaLab", true);
  t->SetBranchAddress("ThetaLab", &ThetaLab, &b_ThetaLab);
  t->SetBranchStatus("PhiLab", true);
  t->SetBranchAddress("PhiLab", &PhiLab, &b_PhiLab);
  t->SetBranchStatus("Ex", "true");
  t->SetBranchAddress("Ex", &Ex, &b_Ex);
  t->SetBranchStatus("ExNoCATS", "true");
  t->SetBranchAddress("ExNoCATS", &ExNoCATS, &b_ExNoCATS);

  TH1F* hExInv = new TH1F("hExInv", "hExInv", 200, -5, 5);
  TH1F* hExMM = new TH1F("hExMM", "hExMM", 200, -5, 5);
  TH1F* hExMMNoCATS = new TH1F("hExMMNoCATS", "hExMMNoCATS", 200, -5, 5);
  TH2F* hE1E2 = new TH2F("hE1E2", "hE1E2", 200, 0, 500, 1000, 0, 500);
  TH2F* hELabThetaLab = new TH2F("hELabThetLab", "hELabThetaLab", 1000, 0, 50, 1000, 0, 500);

  ////////////////////////////////////////////////////////////////////////////////
  int n_entries = t->GetEntries();
  // int n_entries = 1e7;
  for (unsigned int i = 0; i < n_entries; i++) {
    t->GetEntry(i);

    double e_alpha = -1000;
    double theta_alpha = -1000;
    double phi_alpha = -1000;

    double e_li = -1000;
    double theta_li = -1000;
    double phi_li = -1000;

    if (M2->Si_E.size() == 2) {
      if ((cuta->IsInside(M2->CsI_E[0], M2->Si_E[0])) && (cutli->IsInside(M2->CsI_E[1], M2->Si_E[1]))) {
        hExMM->Fill(Ex->at(0));
        hExMMNoCATS->Fill(ExNoCATS->at(0));
        // Particle1 = alpha
        e_alpha = ELab->at(0);
        theta_alpha = ThetaLab->at(0) * M_PI / 180.;
        phi_alpha = PhiLab->at(0);
        if (phi_alpha + 180. > 360)
          phi_alpha = phi_alpha - 180;
        else
          phi_alpha = phi_alpha + 180;
        phi_alpha = phi_alpha * M_PI / 180.;

        // Particle2 = 7li
        e_li = ELab->at(1);
        theta_li = ThetaLab->at(1) * M_PI / 180.;
        phi_li = PhiLab->at(1);
        if (phi_li + 180. > 360)
          phi_li = phi_li - 180;
        else
          phi_li = phi_li + 180;
        phi_li = phi_li * M_PI / 180.;
      }
      else if ((cuta->IsInside(M2->CsI_E[1], M2->Si_E[1])) && (cutli->IsInside(M2->CsI_E[0], M2->Si_E[0]))) {
        hExMM->Fill(Ex->at(1));
        hExMMNoCATS->Fill(ExNoCATS->at(1));
        // Particle2 = alpha
        e_alpha = ELab->at(1);
        theta_alpha = ThetaLab->at(1) * M_PI / 180.;
        phi_alpha = PhiLab->at(1);
        if (phi_alpha + 180. > 360)
          phi_alpha = phi_alpha - 180;
        else
          phi_alpha = phi_alpha + 180;
        phi_alpha = phi_alpha * M_PI / 180.;

        // Particle1 = 7li
        e_li = ELab->at(0);
        theta_li = ThetaLab->at(0) * M_PI / 180.;
        phi_li = PhiLab->at(0);
        if (phi_li + 180. > 360)
          phi_li = phi_li - 180;
        else
          phi_li = phi_li + 180;
        phi_li = phi_li * M_PI / 180.;
      }

      if (e_alpha > 0 && e_li > 0) {
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

        // double etot_alpha = m_alpha + e_alpha;
        // double p_alpha_mag = sqrt(etot_alpha * etot_alpha - m_alpha * m_alpha);
        // LV_alpha.SetPxPyPzE(p_alpha_mag * sin(theta_alpha) * cos(phi_alpha), p_alpha_mag * sin(phi_alpha),
        //                     p_alpha_mag * cos(theta_alpha), etot_alpha);

        ////////////////////////////////////////////////////////////////////////////////
        double gamma_li = e_li / m_7Li + 1;
        double beta_li = sqrt(1 - 1 / (pow(gamma_li, 2.)));
        double p_li_mag = gamma_li * m_7Li * beta_li;
        double etot_li = sqrt(pow(p_li_mag, 2) + pow(m_7Li, 2));

        TVector3 p_li(sin(theta_li) * cos(phi_li), // px
                      sin(theta_li) * sin(phi_li), // py
                      cos(theta_li));              // pz
        p_li.SetMag(p_li_mag);
        LV_7Li.SetPxPyPzE(p_li.x(), p_li.y(), p_li.z(), etot_li);

        //         double etot_li = m_7Li + e_li;
        //         double p_li_mag = sqrt(etot_li * etot_li - m_7Li * m_7Li);
        //         LV_7Li.SetPxPyPzE(p_li_mag * sin(theta_li) * cos(phi_li), p_li_mag * sin(phi_li), p_li_mag *
        //         cos(theta_li),
        //                           etot_li);

        LV_p.SetPxPyPzE(0, 0, 0, mp);
        TLorentzVector LV_total = LV_alpha + LV_7Li - LV_p;

        double ExTot = LV_total.Mag() - m_10Be;
        hExInv->Fill(ExTot);
        hE1E2->Fill(e_alpha, e_li);
        hELabThetaLab->Fill(theta_alpha * 180 / M_PI, e_alpha);
        hELabThetaLab->Fill(theta_li * 180 / M_PI, e_li);
      }
    }
  }
  TCanvas* c1 = new TCanvas();
  gStyle->SetOptStat(0);
  hExInv->Draw();
  hExMM->SetLineColor(kRed);
  hExMM->Draw("same");
  hExMMNoCATS->SetLineColor(kGreen);
  hExMMNoCATS->Draw("same");

  TF1* fitExTot = new TF1("fitExTot", "gaus", hExInv->GetXaxis()->GetXmin(), hExInv->GetXaxis()->GetXmax());
  fitExTot->SetNpx(10000);
  hExInv->Fit(fitExTot, "Q");
  double sigmaExTot = fitExTot->GetParameter(2);
  std::cout << "FWHM of hExInv: " << sigmaExTot*2.35 << " MeV" << std::endl;

  double fitRangeMin = 0.0;
  double fitRangeMax = 15.0;
  TF1* fitEx = new TF1("fitEx", "gaus", fitRangeMin, fitRangeMax);
  fitEx->SetNpx(10000);
  hExMM->Fit(fitEx, "Q", "", -10, 10);
  double sigmaEx = fitEx->GetParameter(2);
  std::cout << "FWHM of hEx: " << sigmaEx*2.35 << " MeV" << std::endl;

  TF1* fitExNoCATS = new TF1("fitExNoCATS", "gaus", fitRangeMin, fitRangeMax);
  fitExNoCATS->SetNpx(10000);
  hExMMNoCATS->Fit(fitExNoCATS, "Q", "", -10, 10);
  double sigmaExNoCATS = fitExNoCATS->GetParameter(2);
  std::cout << "FWHM of hExNoCATS: " << sigmaExNoCATS*2.35 << " MeV" << std::endl;

  TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
  legend->AddEntry(hExInv, "hExInv", "l"); // "l" for line
  legend->AddEntry(hExMM, "hExMM", "l");
  legend->Draw();

  TCanvas* c2 = new TCanvas();
  hELabThetaLab->Draw("col");
}
