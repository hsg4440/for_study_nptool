void GenerateCS() {

  // TGraph* g = new TGraph("cross_sections/12Cd6LiNptool03_70.txt");
  // double scale = 1;
  // TGraph* g = new TGraph("cross_sections/12Cd6Li18MeV.txt");
  // TGraph* g = new TGraph("cross_sections/12Cd6Li18MeV_DWBA.txt");
  // double scale = 1.e3;
  TGraph* g = new TGraph("cross_sections/12CpAlphaNptool.txt");
  double scale = 1;
  g->Draw();

  std::vector<double> vx;
  std::vector<double> vy;

  ofstream file;
  // file.open("cross_sections/12Cd6LiNptool03_70.txt");
  file.open("cross_sections/12Cd6LiNptool18MeV.txt");
  // file.open("cross_sections/test.txt");

  TRandom3* r = new TRandom3();
  for (int i = 0; i < 1800; i++) {

    if (i / 10. < 10 || i / 10. > 90) {
      file << i / 10. << " " << 0 << endl;
      vx.push_back(i / 10.);
      vy.push_back(0);
    }
    else if (g->Eval(i / 10.) < 0) {
      file << i / 10. << " " << 0 << endl;
      vx.push_back(i / 10.);
      vy.push_back(0);
    }
    else {
      file << i / 10. << " " << g->Eval(i / 10.)*scale << endl;
      vx.push_back(i / 10.);
      vy.push_back(g->Eval(i / 10.)*scale);
    }
  }

  TGraph* g2 = new TGraph(vx.size(), &vx[0], &vy[0]);
  g2->Draw("same");
  g2->SetLineColor(kRed);

  double integral = 0;
  int npoints = 0;
  for (unsigned int i = 0; i < g2->GetN(); i++) {
    integral += i * g2->Eval(i);
    npoints++;
    cout << integral << " " << npoints << endl;
  }
  cout << integral / npoints << endl;
}
