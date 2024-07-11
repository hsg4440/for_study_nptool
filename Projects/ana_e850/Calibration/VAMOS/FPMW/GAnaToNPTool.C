void GAnaToNPTool(){

  string rootfilename = "/home/morfouacep/Physics/ganil/pista/analysisenv-e850-2023/Calibs/";
  ifstream ifile;
  string list_filename[8];
  list_filename[0] = "TMW1_X.cal";
  list_filename[1] = "TMW1_Y.cal";
  list_filename[2] = "TMW2_X.cal";
  list_filename[3] = "TMW2_Y.cal";
  list_filename[4] = "FPMW0_X.cal";
  list_filename[5] = "FPMW0_Y.cal";
  list_filename[6] = "FPMW1_X.cal";
  list_filename[7] = "FPMW1_Y.cal";

  for(int i=0; i<4; i++){
    string input_filename = rootfilename + list_filename[i];
    ifile.open(input_filename.c_str());

    string XorY;
    if(i==0 || i==2 || i==4 || i==6) XorY = "X";
    else if(i==1 || i==3 || i==5 || i==7) XorY = "Y";
    int det;
    if(i==0 || i==1) det=1;
    else if(i==2 || i==3) det=2;
    else if(i==4 || i==5) det=3;
    else if(i==6 || i==7) det=4;



    ofstream ofile;
    string output_filename = "npformat/FPMW" + to_string(det) + "_" + XorY + ".cal";
    ofile.open(output_filename.c_str());

    string buffer;
    getline(ifile,buffer);
    getline(ifile,buffer);
    getline(ifile,buffer);
    getline(ifile,buffer);
    double p0, p1, p2, scale;
    int strip=0;
    while(ifile>>p0>>p1>>p2>>scale>>buffer){
      string token = "FPMW_DET" + to_string(det) + "_STRIP"+ XorY + to_string(strip);

      ofile << token << " " << p0 << " " << p1 << " " << p2 << " " << scale << endl;
      strip++;
    }

    ifile.close();
    ofile.close();
  }

}
