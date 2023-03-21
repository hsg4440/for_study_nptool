#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCATSPhysics.h"
#include "TMust2Physics.h"
#include "NPReaction.h"

#include "ROOT/RDataFrame.hxx"
#include <ROOT/RVec.hxx>

#include "Utils.cpp"

#include <Rtypes.h>
#include <TCanvas.h>
#include <algorithm>
#include <iostream>
#include <vector>

void GetSpectrum()
{
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat("nmeruoi");

    ROOT::RDataFrame df("FinalTree", "./RootFiles/FinalTree_12Be.root");
    auto dsc {df.Describe()};
    dsc.Print();

    //--> Read TCutGs
    auto* cutChio {GetGraphCutFromFile("./Cuts/chio_v0.root")};
    auto* cutPID {GetGraphCutFromFile("./../TOF.root", "TOF")};
    
    //--> Get histograms
    // kinematic lines
    auto hKin {df.Histo2D({"hKin", "Kinematic lines;#theta_{Lab} [degree];E_{Vertex} [MeV]", 500, 0., 60., 500, 0., 30.}, "ThetaLab", "ELab")};
    auto hBeam {df.Histo1D({"hBem", "Beam raw TOF", 500, 100., 300.}, "T_CATS1_CAV")};
    // excitation energy
    auto hEex {df.Histo1D("Ex")};
    // beam energy
    auto hBeamE {df.Histo1D("BeamEnergy")};
    
    // Reaction
    NPL::Reaction reaction("12Be(d,3He)11Li@355");
    auto* g {reaction.GetKinematicLine3()};
    
    auto* cRes {new TCanvas("cRes", "Results canvas", 1)};
    cRes->DivideSquare(2);
    cRes->cd(1);
    hKin->DrawClone("colz");
    g->SetLineColor(kRed); g->SetLineWidth(2);
    g->Draw("l same");
    cRes->cd(2);
    hEex->DrawClone();

    auto* cBeam {new TCanvas("cBeam", "Beam canvas", 1)};
    cBeam->DivideSquare(2);
    cBeam->cd(1);
    hBeam->DrawClone();
    cBeam->cd(2);
    hBeamE->DrawClone();
}
