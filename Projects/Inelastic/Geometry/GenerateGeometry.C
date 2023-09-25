double R = 2500; // mm

void GenerateGeometry()
{

  ofstream ofile;
  string ofilename = "Vendeta.detector";
  ofile.open(ofilename.c_str());
  ofile << "\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;

  ofile << "Target " << endl;
  ofile << "  THICKNESS= 10 micrometer" << endl;
  ofile << "  RADIUS= 20 mm" << endl;
  ofile << "  MATERIAL= 238U " << endl;
  ofile << "  ANGLE= 0 deg " << endl;
  ofile << "  X= 0 mm " << endl;
  ofile << "  Y= 0 mm " << endl;
  ofile << "  Z= 0 mm " << endl;

  double Theta = 20;
  double Phi[6] = {0, 20, 40, 140, 160, 180};
  for(int p=0; p<6; p++){
    if(p%2==0) Theta = 20;
    else Theta = 15;
    cout << p << " " << Theta << endl;
    for(int i = 0; i <15; i++){
      ofile << "\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      ofile << "Vendeta" << endl;
      ofile << "  R= " << R << " mm" << endl;
      ofile << "  THETA= " << Theta << " deg" << endl;
      ofile << "  PHI= " << Phi[p] << " deg" << endl;

      //ofile << "\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      //ofile << "Vendeta" << endl;
      //ofile << "  R= " << R << " mm" << endl;
      //ofile << "  THETA= " << -Theta << " deg" << endl;
      //ofile << "  PHI= " << Phi << " deg" << endl;


      Theta += 10;
    }
  }

  ofile.close();


}
