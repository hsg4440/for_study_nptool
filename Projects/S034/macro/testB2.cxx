void testB2(){

  SamuraiFieldMap field;
  double angle = 30*deg;
  double fdc2angle = (59.930-90.0)*deg;
  field.LoadMap(angle,"field_map/180702-2,40T-3000.table.bin",10);
  field.SetFDC2Angle(fdc2angle);
  field.SetFDC2R(3686.77 + 880.745/2.);
  //auto scan  = field.BrhoScan(2.5,10,0.1);
//  scan->Draw();
  //auto h=new TH1D("dB","dB",1000,-2e-1,2e-1); 
 //  new TCanvas();
  auto h=new TH1D("dB","dB",1000,-0.1,0.1); 
  auto hX=new TH2D("dBX","dBX",500,-100,100,500,-0.01,0.01); 
  auto hY=new TH2D("dBY","dBY",500,-100,100,500,-0.01,0.01); 
  double r = gRandom->Uniform(-1,1);

  for(unsigned int i = 0 ; i < 100 ; i++){
    double X = gRandom->Uniform(-100,100);
    double Y = gRandom->Uniform(-100,100);
    TVector3 p(X,Y,-3500); 
    //TVector3 d(gRandom->Uniform(-0.05,0.05),gRandom->Uniform(-0.05,0.05),1); 
 //
  //TVector3 p(0,0,-3500); 
  TVector3 d(0,0,1); 
  double b = gRandom->Uniform(3,7);
  std::vector< TVector3 > track = field.Propagate(b,p,d,false);
  // rotate from lab angle to FDC2 frame
  track.back().RotateY(-fdc2angle+angle);
  double br=field.FindBrho(p,d,track.back(),d);
  //double pc = 100*(br-b)/b;

  double pc = br-b;
//  cout <<"  " <<  100*(br-b)/b  << " " << br << " " << b << endl;
  h->Fill(pc);
  hX->Fill(X,pc);
  hY->Fill(Y,pc);
  }
  
  h->Draw();
  new TCanvas();
  hX->Draw("colz");
  new TCanvas();
  hY->Draw("colz");

} 
