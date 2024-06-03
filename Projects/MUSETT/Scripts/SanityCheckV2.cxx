#include "SanityCheckV2.hh"


void FillHit(const TMUSETTPhysics* event, const AnalysisConfig& config)
{
  config.EventMultHist->Fill(event->EventMultiplicity);
  for(int j = 0; j < event->EventMultiplicity; j++)
  {
    if((event->DSSD_E[j]>config.E_min) && (event->DSSD_E[j] < config.E_max))
    {
      if(!config.isDead(event->DetectorNumber[j],event->DSSD_X[j],event->DSSD_Y[j]))
      {
        config.HitPatterns[event->DetectorNumber[j]]->Fill(event->DSSD_X[j], event->DSSD_Y[j]);
        config.Xhits[event->DetectorNumber[j]]->Fill(event->DSSD_X[j]);
        config.Yhits[event->DetectorNumber[j]]->Fill(event->DSSD_Y[j]);
      }
    }
    config.SumEnergy[event->DetectorNumber[j]]->Fill(event->DSSD_E[j]);
  }
}

void FillCorrelation(const TMUSETTPhysics* event, const AnalysisConfig& config)
{
  int mult = event->EventMultiplicity;
  int detA, detB;
  for(int j1 = 0; j1 < mult; j1++)
  {
    detA = event->DetectorNumber[j1];
    for(int j2 = j1+1; j2 < mult; j2++)
    {
      detB = event->DetectorNumber[j2];
      config.E1E2_detAdetB[detA][detB]->Fill(event->DSSD_E[j1],event->DSSD_E[j2]);
      config.E1E2_all->Fill(event->DSSD_E[j1],event->DSSD_E[j2]);
      if(detA!=detB)
      {
        config.E1E2_all_diffDet->Fill(event->DSSD_E[j1],event->DSSD_E[j2]);
      }
    }
  }
}

void AnalyzeData(const AnalysisConfig& config) {

    TFile *inputFile = new TFile(config.inputFileName.c_str(), "READ"); // Open your ROOT file for reading
    TTree *tree = (TTree*)inputFile->Get("PhysicsTree");  // Replace "your_tree_name" with the name of your TTree

    TObjArray* branches = tree->GetListOfBranches();
    TMUSETTPhysics* Mumu = new TMUSETTPhysics();
    TBranch *branch = tree->GetBranch("MUSETT");
    branch->SetAddress(&Mumu);


    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        FillHit(Mumu,config);
        if(Mumu->EventMultiplicity>=2)FillCorrelation(Mumu,config);
      }
}

void SanityCheckV2() {
    AnalysisConfig config;
    config.inputFileName = "/Users/lh270370/Software/nptool/Outputs/Analysis/Simu/Isolde/energy_beam_test/222Ra_fixedGeometry_1e6__energyBeam50_physics.root";
    config.outputFileName = "/Users/lh270370/Software/nptool/Outputs/Analysis/Simu/Isolde/energy_beam_test/histoSanity_222Ra_fixedGeometry_1e6__energyBeam50_physics.root";
    config.E_min = 6.4;
    config.E_max = 6.6;
    config.tau = 100;
    config.T = 1e+09;
    AnalyzeData(config);
    config.WriteHistograms();
    config.ClearHistograms();
}
