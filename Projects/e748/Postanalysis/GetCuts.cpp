#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCATSPhysics.h"
#include "TMust2Physics.h"
#include "Math/Point3D.h"

#include "ROOT/RDataFrame.hxx"
#include <ROOT/RVec.hxx>

#include "Utils.cpp"

#include <Rtypes.h>
#include <TCanvas.h>
#include <algorithm>
#include <iostream>
#include <vector>

void GetCuts()
{
    ROOT::EnableImplicitMT();
    /////////////////////////////////////////////////////
    ///READ DATA
    gStyle->SetOptStat("nmeruoi");
    
    auto* inFile {new TFile("/home/miguel/nptool/Projects/e748/output/analysis/Physics_12Be.root")};
    //inFile->ls();
    auto* tree {inFile->Get<TTree>("PhysicsTree")};
    //tree->Print();

    //load also chio chain
    auto* chioChain {new TChain("Numexo2")};//runs for 12Be
    std::set<int> runs {315, 316, 317, 318, 320, 321, 323, 325, 326, 327, 328, 329, 330, 331, 339, 341, 342, 346, 347, 348};
    for(const auto& run : runs)
        chioChain->Add(TString::Format("/home/miguel/nptool/Projects/e748/Data/Merged/run_%04d.root", run));
    tree->AddFriend(chioChain);

    ////// DATAFRAME
    ROOT::RDataFrame df(*tree);
    //Apply base cut
    auto baseCut = [](const TCATSPhysics& cats,
                      const std::vector<short>& vtel, const std::vector<double>& vCsIE,
                      unsigned short chioE, double QPlast, double T_CATS1_CAV)
    {
        //0--> Only one hit per event!
        bool condSize {vtel.size() == 1};
        //1--> Position on target should be inside its dimensions
        bool condPosition {cats.PositionOnTargetY > -10 && cats.PositionOnTargetY < 10.
            && cats.PositionOnTargetX > -20. && cats.PositionOnTargetX < 20.};
        //2--> Telescope number lower than 5
        bool condTelescopeN {std::all_of(vtel.begin(), vtel.end(),
                                         [](const auto& t){return t < 5;})};
        //3--> CsI energy == 0!!
        bool condCsIEnergy {std::all_of(vCsIE.begin(), vCsIE.end(),
                                        [](const auto& e){return e < 0;})};
        //4--> Recoil measured in IC and plastic
        bool condRecoil {chioE > 0 && QPlast > 0};
        //5--> TOF for beam, raw cut
        bool condBeam {T_CATS1_CAV > 160. && T_CATS1_CAV < 235.};
        return (condSize && condPosition && condTelescopeN && condCsIEnergy && condRecoil && condBeam);
    };
    
    auto gated {df.Filter(baseCut,
                          {"CATS",
                           "Must2Telescopes", "Must2CsIE",
                           "IC_E", "QPlast", "T_CATS1_CAV"})};
    
    //Build ID for light recoil, based on TOF between CATS and MUST2
    gated = gated.Define("TOFLight",
                         [](const std::vector<short>& vTel, const std::vector<double>& vSiT,
                            const double& timeCorr)
                         {
                             ROOT::RVecD v {};//serialize use of RVec even though we are working with a size=1 vector
                             for(int i = 0; i < vTel.size(); i++)
                             {
                                 auto val1 {vSiT.at(i) + timeCorr};
                                 auto tn {vTel.at(i)};
                                 double opt {};
                                 if(tn == 1)
                                     opt = -2.521;
                                 else if(tn == 2)
                                     opt = 0.148;
                                 else if(tn == 3)
                                     opt = -1.922;
                                 else if(tn == 4)
                                     opt = -7.176;
                                 v.push_back(val1 + opt);
                             }
                             return v;
                         },
                         {"Must2Telescopes", "Must2SiT", "TimeCorr"});

    // //HISTOGRAMS
    auto hCHIO { gated.Histo2D({"hCHIO", "Heavy ID;QPlast;IC_E", 4000, 0., 16000., 900, 0., 900.}, "QPlast", "IC_E")};
    auto hPID {gated.Histo2D({"hPID", "PID;E_{Si} [MeV];TOF_{Heavy}[au]", 1000, 0., 30., 1000, 460., 580.}, "Must2SiE", "TOFLight")};

    //--> Read TCutGs
    auto* cutChio {GetGraphCutFromFile("./Cuts/chio_v0.root")};
    auto* cutPID {GetGraphCutFromFile("./../TOF.root", "TOF")};
    
    //plotting
    auto* cID {new TCanvas("cID", "ID canvas", 1)};
    cID->DivideSquare(2);
    cID->cd(1);
    hCHIO->DrawClone("colz");
    cutChio->SetLineColor(kRed); cutChio->SetLineWidth(2);
    cutChio->Draw("same");
    cID->cd(2);
    hPID->DrawClone("colz");
    cutPID->SetLineColor(kRed); cutPID->SetLineWidth(2);
    cutPID->Draw("same");
    
    //save to disk gated TTree
    gated.Snapshot("BaseCutTree", "./RootFiles/BaseCutTree_12Be.root");
}
