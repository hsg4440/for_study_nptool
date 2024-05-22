void PlotYieldEvolution(int A_asked=147){
  int Energy[10] = {6,7,8,9,10,11,12,13,14,15};
  TGraphErrors* gevol = new TGraphErrors();

  for(int i = 0; i<10; i++){
    ifstream ifile;
    string input_filename = "Mass_Yield_" + to_string(Energy[i]) + "MeV.dat";
    ifile.open(input_filename.c_str());

    int A;
    double yield;
    double yield_err;
    while(ifile>>A>>yield>>yield_err){
      if(A==A_asked){
        gevol->SetPoint(i,Energy[i]-5.5,yield);
        gevol->SetPointError(i,0,yield_err);
      }
    }
  }

  gevol->SetMarkerStyle(8);
  gevol->SetMarkerSize(1);
  gevol->Draw("ap");
}
