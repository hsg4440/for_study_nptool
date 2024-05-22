void DrawEx() {
  TChain* t = new TChain("PhysicsTree");

  TFile* fcutDEE3 = NULL;
  TCutG* cutDEE3 = NULL;
  TFile* fcutDEE4 = NULL;
  TCutG* cutDEE4 = NULL;

  TFile* fcutToF3 = NULL;
  TCutG* cutToF3 = NULL;
  bool ToF3 = false;
  TFile* fcutToF4 = NULL;
  TCutG* cutToF4 = NULL;
  bool ToF4 = false;

  ////////////////////////////////////////////////////////////////////////////////
  // t->Add("../../Outputs/Analysis/Ana10Bepalpha_38AMeV_CH2.root"); // Name of simulated file
  // NPL::Reaction r("10Be(p,4He)7Li@380");

  // // YOUR 3 AND 4 CUT
  // fcutDEE3 = new TFile("cuts/CUT_alpha.root", "read");
  // cutDEE3 = (TCutG*)gDirectory->FindObjectAny("CUT_alpha");
  // fcutDEE4 = new TFile("cuts/CUT_7Li.root", "read");
  // cutDEE4 = (TCutG*)gDirectory->FindObjectAny("CUT_7Li");

  ////////////////////////////////////////////////////////////////////////////////
  // t->Add("../../Outputs/Analysis/Ana10Bed6Li_20AMeV_CD2.root"); // Name of simulated file
  // NPL::Reaction r("10Be(d,6Li)6He@200");
  
  t->Add("../../Outputs/Analysis/Ana10Bed6Li_38AMeV_CD2.root"); // Name of simulated file
  NPL::Reaction r("10Be(d,6Li)6He@380");

  // YOUR 3 AND 4 CUT
  fcutDEE3 = new TFile("cuts/CUT_6Li.root", "read");
  cutDEE3 = (TCutG*)gDirectory->FindObjectAny("CUT_6Li");
  fcutDEE4 = new TFile("cuts/CUT_6He.root", "read");
  cutDEE4 = (TCutG*)gDirectory->FindObjectAny("CUT_6He");

  fcutToF3 = new TFile("cuts/CUTToF_6Li.root", "read");
  cutToF3 = (TCutG*)gDirectory->FindObjectAny("CUTToF_6Li");
  ToF3 = true;

  ////////////////////////////////////////////////////////////////////////////////
  TLorentzVector LV3, LV4, LVTarget;
  double mBeam = r.GetParticle1()->Mass();
  double mTarget = r.GetParticle2()->Mass();
  double mParticle3 = r.GetParticle3()->Mass();
  double mParticle4 = r.GetParticle4()->Mass();

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

  TH1F* hExInv = new TH1F("hExInv", "hExInv", 200, -10, 10);
  TH1F* hExMM = new TH1F("hExMM", "hExMM", 200, -10, 10);
  TH1F* hExMMNoCATS = new TH1F("hExMMNoCATS", "hExMMNoCATS", 200, -10, 10);
  TH2F* hE1E2 = new TH2F("hE1E2", "hE1E2", 200, 0, 500, 1000, 0, 500);
  TH2F* hELabThetaLab = new TH2F("hELabThetLab", "hELabThetaLab", 1000, 0, 50, 1000, 0, 500);
  TH2F* hExELab = new TH2F("hExELab", "hExELab", 1000, 0, 500, 1000, -10, 10);

  ////////////////////////////////////////////////////////////////////////////////
  int n_entries = t->GetEntries();
  // int n_entries = 1e7;
  for (unsigned int i = 0; i < n_entries; i++) {
    t->GetEntry(i);

    double e3 = -1000;
    double theta3 = -1000;
    double phi3 = -1000;

    double e4 = -1000;
    double theta4 = -1000;
    double phi4 = -1000;

    int Particle3 = -1000;
    int Particle4 = -1000;
    int number_of_match = 0;

    if (M2->Si_E.size() == 2) {
      if (M2->CsI_E[0] > 0) {
        if (cutDEE3->IsInside(M2->CsI_E[0], M2->Si_E[0])) {
          Particle3 = 0;
          number_of_match++;
        }
      }
      else if (M2->CsI_E[0] < 0 && ToF3 == true) {
        if ((cutDEE3->IsInside(M2->CsI_E[0], M2->Si_E[0])) || (cutToF3->IsInside(M2->Si_E[0], M2->Si_T[0]))) {
          Particle3 = 0;
          number_of_match++;
        }
      }
      if (M2->CsI_E[1] > 0) {
        if (cutDEE4->IsInside(M2->CsI_E[1], M2->Si_E[1])) {
          Particle4 = 1;
          number_of_match++;
        }
      }
      else if (M2->CsI_E[1] < 0 && ToF4 == true) {
        if ((cutDEE4->IsInside(M2->CsI_E[1], M2->Si_E[1])) || (cutToF4->IsInside(M2->Si_E[1], M2->Si_T[1]))) {
          Particle4 = 1;
          number_of_match++;
        }
      }
      if (cutDEE3->IsInside(M2->CsI_E[1], M2->Si_E[1])) {
        Particle3 = 1;
        number_of_match++;
      }
      else if (M2->CsI_E[1] < 0 && ToF3 == true) {
        if ((cutDEE3->IsInside(M2->CsI_E[1], M2->Si_E[1])) || (cutToF3->IsInside(M2->Si_E[1], M2->Si_T[1]))) {
          Particle3 = 1;
          number_of_match++;
        }
      }
      if (M2->CsI_E[0] > 0) {
        if (cutDEE4->IsInside(M2->CsI_E[0], M2->Si_E[0])) {
          Particle4 = 0;
          number_of_match++;
        }
      }
      else if (M2->CsI_E[0] < 0 && ToF4 == true) {
        if ((cutDEE4->IsInside(M2->CsI_E[0], M2->Si_E[0])) || (cutToF4->IsInside(M2->Si_E[0], M2->Si_T[0]))) {
          Particle4 = 0;
          number_of_match++;
        }
      }

      if (number_of_match == 2) {
        if (Particle3 == 0 && Particle4 == 1) {
          hExMM->Fill(Ex->at(0));
          hExMMNoCATS->Fill(ExNoCATS->at(0));
          e3 = ELab->at(0);
          theta3 = ThetaLab->at(0) * M_PI / 180.;
          phi3 = PhiLab->at(0);
          if (phi3 + 180. > 360)
            phi3 = phi3 - 180;
          else
            phi3 = phi3 + 180;
          phi3 = phi3 * M_PI / 180.;

          e4 = ELab->at(1);
          theta4 = ThetaLab->at(1) * M_PI / 180.;
          phi4 = PhiLab->at(1);
          if (phi4 + 180. > 360)
            phi4 = phi4 - 180;
          else
            phi4 = phi4 + 180;
          phi4 = phi4 * M_PI / 180.;
        }
        if (Particle3 == 1 && Particle4 == 0) {
          hExMM->Fill(Ex->at(1));
          hExMMNoCATS->Fill(ExNoCATS->at(1));
          e3 = ELab->at(1);
          theta3 = ThetaLab->at(1) * M_PI / 180.;
          phi3 = PhiLab->at(1);
          if (phi3 + 180. > 360)
            phi3 = phi3 - 180;
          else
            phi3 = phi3 + 180;
          phi3 = phi3 * M_PI / 180.;

          e4 = ELab->at(0);
          theta4 = ThetaLab->at(0) * M_PI / 180.;
          phi4 = PhiLab->at(0);
          if (phi4 + 180. > 360)
            phi4 = phi4 - 180;
          else
            phi4 = phi4 + 180;
          phi4 = phi4 * M_PI / 180.;
        }
      }

      if (e3 > 0 && e4 > 0) {
        ////////////////////////////////////////////////////////////////////////////////
        double gamma3 = e3 / mParticle3 + 1.;
        double beta3 = sqrt(1 - 1 / (pow(gamma3, 2.)));
        double p3_mag = gamma3 * mParticle3 * beta3;
        double etot3 = sqrt(pow(p3_mag, 2) + pow(mParticle3, 2));

        TVector3 p3(sin(theta3) * cos(phi3), // px
                    sin(theta3) * sin(phi3), // py
                    cos(theta3));            // pz
        p3.SetMag(p3_mag);
        LV3.SetPxPyPzE(p3.x(), p3.y(), p3.z(), etot3);

        ////////////////////////////////////////////////////////////////////////////////
        double gamma4 = e4 / mParticle4 + 1;
        double beta4 = sqrt(1 - 1 / (pow(gamma4, 2.)));
        double p4_mag = gamma4 * mParticle4 * beta4;
        double etot4 = sqrt(pow(p4_mag, 2) + pow(mParticle4, 2));

        TVector3 p4(sin(theta4) * cos(phi4), // px
                    sin(theta4) * sin(phi4), // py
                    cos(theta4));            // pz
        p4.SetMag(p4_mag);
        LV4.SetPxPyPzE(p4.x(), p4.y(), p4.z(), etot4);

        LVTarget.SetPxPyPzE(0, 0, 0, mTarget);
        TLorentzVector LV_total = LV3 + LV4 - LVTarget;

        double ExTot = LV_total.Mag() - mBeam;
        hExInv->Fill(ExTot);
        hE1E2->Fill(e3, e4);
        hELabThetaLab->Fill(theta3 * 180 / M_PI, e3);
        hELabThetaLab->Fill(theta4 * 180 / M_PI, e4);
        hExELab->Fill(e3, ExTot);
        hExELab->Fill(e4, ExTot);
      }
    }
  }
  TCanvas* c1 = new TCanvas();
  gStyle->SetOptStat(0);
  hExInv->Draw();
  hExMM->SetLineColor(kRed);
  hExMM->Draw("same");

  TF1* fitExTot = new TF1("fitExTot", "gaus", hExInv->GetXaxis()->GetXmin(), hExInv->GetXaxis()->GetXmax());
  fitExTot->SetNpx(10000);
  hExInv->Fit(fitExTot, "Q");
  double sigmaExTot = fitExTot->GetParameter(2);
  std::cout << "FWHM of hExInv: " << sigmaExTot * 2.35 << " MeV" << std::endl;

  double fitRangeMin = 0.0;
  double fitRangeMax = 15.0;
  TF1* fitEx = new TF1("fitEx", "gaus", fitRangeMin, fitRangeMax);
  fitEx->SetNpx(10000);
  hExMM->Fit(fitEx, "Q", "", -10, 10);
  double sigmaEx = fitEx->GetParameter(2);
  std::cout << "FWHM of hEx: " << sigmaEx * 2.35 << " MeV" << std::endl;

  TLegend* legend = new TLegend(0.75, 0.75, 0.85, 0.85);
  legend->AddEntry(hExInv, "hExInv", "l"); // "l" for 4ne
  legend->AddEntry(hExMM, "hExMM", "l");
  legend->Draw();

  TCanvas* c2 = new TCanvas();
  hELabThetaLab->Draw("col");
}
