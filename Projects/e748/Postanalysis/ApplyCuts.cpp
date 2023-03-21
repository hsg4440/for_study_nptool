#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCATSPhysics.h"
#include "TMust2Physics.h"

#include "ROOT/RDataFrame.hxx"
#include <ROOT/RVec.hxx>

#include "Utils.cpp"

#include <Rtypes.h>
#include <TCanvas.h>
#include <algorithm>
#include <iostream>
#include <vector>

void ApplyCuts()
{
    ROOT::EnableImplicitMT();
    gStyle->SetOptStat("nmeruoi");
    
    ROOT::RDataFrame df("BaseCutTree", "./RootFiles/BaseCutTree_12Be.root");
    auto dsc {df.Describe()};
    //dsc.Print();

    //--> Read TCutGs
    auto* cutChio {GetGraphCutFromFile("./Cuts/chio_v0.root")};
    auto* cutPID {GetGraphCutFromFile("./../TOF.root", "TOF")};
    //--> Apply them. Must apply cuts to the vectors stored in classes
    auto gated {df.Filter([&](const ROOT::RVecD& vSiE, const ROOT::RVecD& pid,
                              const unsigned short& chioE, const double& Qplast)
    {
        std::vector<bool> vCondPID {};
        for(int i = 0; i < vSiE.size(); i++)
        {
            vCondPID.push_back(cutPID->IsInside(vSiE.at(i), pid.at(i)));
        }
        //all hits in MUST2 should be inside TCutG
        bool condPID {std::all_of(vCondPID.begin(), vCondPID.end(), [](const bool& e){return e;})};
        //same but for heavy ID
        bool condChio {static_cast<bool>(cutChio->IsInside(Qplast, chioE))};
        return (condChio && condPID);
    },
                          {"Must2SiE", "TOFLight", "IC_E", "QPlast"})};
    
   
    //--> Write to file
    gated.Snapshot("FinalTree", "./RootFiles/FinalTree_12Be.root");    
}
