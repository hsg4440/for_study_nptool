#include <ROOT/RDataFrame.hxx>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>

#include "TMust2Physics.h"
#include "Utils.cpp"

#include <set>

void AnalyzeBeam()
{
    ROOT::EnableImplicitMT();
    //gStyle->SetOptStat("nmeruoi");

    //ROOT::RDataFrame df("BaseCutTree", "./RootFiles/BaseCutTree_12Be.root");
    auto df {ReadAll12BeData()};

    auto gated {df.Filter(
                          [](const TMust2Physics& must2)
                          {
                              return (must2.TelescopeNumber.size() == 0);
                          },
                          {"MUST2"})};
    auto hCATSchio {gated.Histo2D({"hCATSchio", "Raw TOF vs IC", 700, 0., 300., 900, 0, 900}, "T_CATS1_CAV", "IC_E")};
    auto hCATSplast {gated.Histo2D({"hCATSplast", "Raw TOF vs QPlast", 700, 0., 300., 4000, 0, 16000}, "T_CATS1_CAV", "QPlast")};
    auto hCHIOplast {gated.Histo2D({"hCHIOplast", "CHIO vs QPlast", 900, 0, 900, 4000, 0, 16000}, "IC_E", "QPlast")};
    
    //plotting
    auto* c1 {new TCanvas("c1", "Beam ID")};
    c1->DivideSquare(4);
    c1->cd(1);
    hCATSchio->DrawClone("colz");
    c1->cd(2);
    hCATSplast->DrawClone("colz");
    c1->cd(3);
    hCHIOplast->DrawClone("colz");
}
