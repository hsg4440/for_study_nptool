void DrawCarbonBackground() {

  TFile* f = new TFile(
      "/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10BeCH2LongLong.root", "open");
      // "/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bepalpha.root", "open");
      // "/Users/valerian/Software/nptool_gitlab/nptool/Outputs/Simulation/MugastAtLise10Bed6LiThick.root", "open");
  TTree* t = (TTree*)f->FindObjectAny("SimulatedTree");

  TFile* fcuta = new TFile("cuta.root", "open");
  TCutG* cuta = (TCutG*)fcuta->FindObjectAny("cuta");
  // TFile* fcuta = new TFile("cut6He.root", "open");
  // TCutG* cuta = (TCutG*)fcuta->FindObjectAny("cut6he");
  // TFile* fcuta = new TFile("cutapa.root", "open");
  // TCutG* cuta = (TCutG*)fcuta->FindObjectAny("cuta");
  // TFile* fcuta = new TFile("cut6hed6li.root", "open");
  // TCutG* cuta = (TCutG*)fcuta->FindObjectAny("cut6he");

  TFile* fcutli = new TFile("cutli.root", "open");
  TCutG* cutli = (TCutG*)fcutli->FindObjectAny("cutli");
  // TFile* fcutli = new TFile("cut6Li.root", "open");
  // TCutG* cutli = (TCutG*)fcutli->FindObjectAny("cut6li");
  // TFile* fcutli = new TFile("cutlipa.root", "open");
  // TCutG* cutli = (TCutG*)fcutli->FindObjectAny("cutli");
  // TFile* fcutli = new TFile("cut6lid6li.root", "open");
  // TCutG* cutli = (TCutG*)fcutli->FindObjectAny("cut6li");

  TInteractionCoordinates* TI = new TInteractionCoordinates();
  t->SetBranchStatus("*", false);
  t->SetBranchStatus("InteractionCoordinates", true);
  t->SetBranchAddress("InteractionCoordinates", &TI);
  TMust2Data* M2 = new TMust2Data();
  t->SetBranchStatus("MUST2", true);
  t->SetBranchAddress("MUST2", &M2);

  TH2F* h = new TH2F("h", "h", 1000, 0, 30, 1000, 0, 500);
  TH2F* h2 = new TH2F("h2", "h2", 1000, 0, 30, 1000, 0, 500);
  TH2F* hE = new TH2F("hE", "hE", 1000, 0, 500, 1000, 0, 500);
  TH2F* hT = new TH2F("hT", "hT", 1000, 0, 100, 1000, 0, 100);

  unsigned int entries = t->GetEntries();

  int n = 0;
  for (unsigned int i = 0; i < entries; i++) {
    t->GetEntry(i);

    // unsigned int DetSize = TI->fDetected_Particle_Name.size();
    // for (unsigned int j = 0; j < DetSize; j++) {
    //   if (TI->fDetected_Particle_Name[j] == "alpha" && TI->fDetected_Position_Z[j] >
    //   250&&TI->fDetected_Energy[j]>0.5) {
    //   // if (TI->fDetected_Particle_Name[j] == "6He" && TI->fDetected_Position_Z[j] > 250) {
    //     for (unsigned int k = 0; k < DetSize; k++) {
    //       if (TI->fDetected_Particle_Name[k] == "7Li" && TI->fDetected_Position_Z[k] >
    //       250&&TI->fDetected_Energy[k]>0.5) {
    //       // if (TI->fDetected_Particle_Name[k] == "6Li" && TI->fDetected_Position_Z[k] > 250) {
    //         h->Fill(TI->fDetected_Angle_Theta[j], TI->fDetected_Energy[j]);
    //         h2->Fill(TI->fDetected_Angle_Theta[k], TI->fDetected_Energy[k]);
    //         cout << j << " " << k << endl;
    //         n++;
    //         cout << n << " "
    //              << "alpha and 7Li" << endl;
    //       }
    //     }
    //   }
    // }
    // }

    unsigned int M2Size = M2->fMM_StripXE_Energy.size();
    // unsigned int M2Size = M2->fMM_CsIE_Energy.size();
    // cout << M2Size << endl;
    for (unsigned int j = 0; j < M2Size; j++) {
      if (M2->fMM_CsIE_Energy.size() > 0) {
        // cout << 1 << " " << M2->fMM_CsIE_Energy[j] << endl;
        if (cuta->IsInside(M2->fMM_CsIE_Energy[j] * 500. / 16384., M2->fMM_StripXE_Energy[j] * 63. / 8192 - 63.)) {
          for (unsigned int k = 0; k < M2Size; k++) {
            // cout << 2 << endl;
            if (cutli->IsInside(M2->fMM_CsIE_Energy[k] * 500. / 16384., M2->fMM_StripXE_Energy[k] * 63. / 8192 - 63.)) {
              h->Fill(TI->fDetected_Angle_Theta[j],
                      M2->fMM_CsIE_Energy[j] * 500. / 16384. + M2->fMM_StripXE_Energy[j] * 63. / 8192 - 63.);
              h2->Fill(TI->fDetected_Angle_Theta[k],
                       M2->fMM_CsIE_Energy[k] * 500. / 16384. + M2->fMM_StripXE_Energy[k] * 63. / 8192 - 63.);

              hT->Fill(TI->fDetected_Angle_Theta[j], TI->fDetected_Angle_Theta[k]);
              hE->Fill(M2->fMM_CsIE_Energy[j] * 500. / 16384. + M2->fMM_StripXE_Energy[j] * 63. / 8192 - 63.,
                       M2->fMM_CsIE_Energy[k] * 500. / 16384. + M2->fMM_StripXE_Energy[k] * 63. / 8192 - 63.);

              cout << j << " " << k << endl;
              n++;
              cout << n
                   << " "
                   // << "alpha and 7Li" << endl;
                   << "6He and 6Li" << endl;
            }
          }
        }
      }
    }
  }

  TCanvas* c1 = new TCanvas;

  // NPL::Reaction b("10Be(p,4He)7Li@380MeV");
  NPL::Reaction b("10Be(d,6Li)6He@300MeV");
  h->Draw("");
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kRed);
  h2->Draw("Same");
  h2->SetMarkerStyle(20);
  TGraph* g = b.GetKinematicLine3();
  g->SetLineColor(kRed);
  g->Draw("same");
  b.GetKinematicLine4()->Draw("same");

  TCanvas* c2 = new TCanvas;
  c2->Divide(2, 1);
  c2->cd(1);
  hE->Draw("col");
  c2->cd(2);
  hT->Draw("col");
}
