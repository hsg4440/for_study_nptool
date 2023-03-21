#include <ROOT/RDataFrame.hxx>
#include <TCutG.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

#include <stdexcept>
#include <string>
#include <set>

TCutG* GetGraphCutFromFile(const std::string& fileName, const std::string& name = "CUTG")
{
    auto* file {new TFile(fileName.c_str())};
    auto* cut {file->Get<TCutG>(name.c_str())};
    if(!cut)
        throw std::runtime_error("Error loading TCutG: Check file existence or name!");
    return cut;
    
}

ROOT::RDataFrame ReadAll12BeData()
{
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
    return ROOT::RDataFrame(*tree);
}
