#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "TCATSPhysics.h"
#include "TMust2Physics.h"
#include "TModularLeafPhysics.h"

#include "ROOT/RDataFrame.hxx"

void AnalyzeCATS()
{
    ROOT::EnableImplicitMT();
    
    auto* inFile {new TFile("/home/miguel/nptool/Projects/e748/output/analysis/Test_12Be.root")};
    //inFile->ls();
    auto* tree {inFile->Get<TTree>("PhysicsTree")};
    //tree->Print();

    ROOT::RDataFrame df(*tree);
    auto cond {df.Filter("CATS.PositionOnTargetX != -1000 && CATS.PositionOnTargetY != -1000")};
    auto hPos {cond.Define("x", "CATS.PositionOnTargetX").Define("y", "CATS.PositionOnTargetY")
        .Histo2D({"hPos", "Position on CATS;X;Y", 100, -70., 70., 100, -70., 70.}, "x", "y")};
    hPos->DrawClone("colz");

    auto hCATS1Cal {cond.Histo1D("T_CATS1_CAV_Cal")};
    hCATS1Cal->DrawClone();

    auto hCATSBeamEnergy {cond.Filter("T_CATS1_CAV > 0")
        .Histo2D({"hCATSBeamEnergy", "CATS against BeamEnergy", 100, 0., 300., 200, 0., 500.}, "T_CATS1_CAV", "BeamEnergy")};
    hCATSBeamEnergy->DrawClone("colz");
  
}
