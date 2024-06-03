// AnalysisConfig.hh
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <cmath> // For sqrt
#include <map>
#include <algorithm> // For std::min and std::max
#include <TCanvas.h>
#include <iostream>
#include "TGraphErrors.h"
#include "TString.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <string>
#include <set>
#include "TMUSETTPhysics.h"

#ifndef ANALYSISCONFIG_HH
#define ANALYSISCONFIG_HH

const double Emin_222Ra = 6.46;
const double Emax_222Ra = 6.50;

const double Emin_220Ra = 7.40;
const double Emax_220Ra = 7.50;



class AnalysisConfig {
public:
    std::string inputFileName;
    std::string outputFileName;
    double E_min;
    double E_max;
    double tau;
    double T;
    int NDet ;
    int NStrip;

    std::vector<TH2F*> HitPatterns;
    std::vector<TH1F*> Xhits;
    std::vector<TH1F*> Yhits;
    std::vector<TH1F*> SumEnergy;
    std::vector<std::vector<TH2F*>> E1E2_detAdetB;
    TH2F* E1E2_all;
    TH2F* E1E2_all_diffDet;
    TH1F* EventMultHist;
    std::map<std::pair<int, int>, int> detPairToIndexMap;

    std::map<int, std::set<int>> deadStripsX = {
        {0, {10, 26, 28}},
        {1, {}},
        {2, {104}},
        {3, {22}}
    };

    // Declaration for dead strips along the Y axis
    std::map<int, std::set<int>> deadStripsY = {
        {0, {30}},
        {1, {12}},
        {2, {118}},
        {3, {}}
    };
    AnalysisConfig() : E_min(0), E_max(100), tau(1), T(100), NDet(4), NStrip(128) {
        InitializeHistograms();
    }

    ~AnalysisConfig() {
        ClearHistograms();
    }

    void InitializeHistograms() {
        // Initialize the histograms for each detector
        //ClearHistograms();
        EventMultHist = new TH1F("EventMultiplicities","Multiplicity of events", 20,0,20);
        for (int det = 0; det < NDet; ++det) {
            HitPatterns.push_back(new TH2F(Form("HitPattern_Det%d", det),
                                           Form("Spatial Hit Distribution for Detector %d;Strip X;Strip Y", det),
                                           NStrip, 0, NStrip, NStrip, 0, NStrip));
            Xhits.push_back(new TH1F(Form("Xhits_Det%d", det),
                                     Form("X Hits Distribution for Detector %d;Strip X;Hits Count", det),
                                     NStrip, 0, NStrip));
            Yhits.push_back(new TH1F(Form("Yhits_Det%d", det),
                                     Form("Y Hits Distribution for Detector %d;Strip Y;Hits Count", det),
                                     NStrip, 0, NStrip));
            SumEnergy.push_back(new TH1F(Form("SumEnergy_Det%d", det),
                                     Form("Summed Energy for Hits in Detector %d;Energy (MeV);Event Count", det),
                                     400000, 0, 20));
        }

        E1E2_detAdetB.resize(NDet);
        for (int detA = 0; detA < NDet; detA++) {
          E1E2_detAdetB[detA].resize(NDet);
          }
        // Initialize the 2D histograms for energy correlations between different detectors
        for (int detA = 0; detA < NDet; detA++) {
            for (int detB = detA; detB < NDet; detB++) {
                std::string histName = Form("E1E2_det%d_det%d", detA, detB);
                std::string histTitle = Form("Energy Correlation between Detectors %d and %d;Energy E%d (MeV);Energy E%d (MeV)", detA, detB, detA, detB);
                E1E2_detAdetB[detA][detB] = new TH2F(histName.c_str(), histTitle.c_str(), 200, 0, 20, 200, 0, 20);
            }
        }
        E1E2_all = new TH2F("E1E2_all", "Energy correlation betwenn all detectors;Energy E__A( MeV);Energy E_B (MeV)", 200,0,20,200,0,20);
        E1E2_all_diffDet = new TH2F("E1E2_all_diffDet", "Energy correlation betwenn all detectors detA!= detB;Energy E_A (MeV);Energy E_B (MeV)", 200,0,20,200,0,20);
    }

    void ClearHistograms() {
        // Properly delete all histogram objects to avoid memory leaks
        for (auto& hist : HitPatterns) delete hist;
        for (auto& hist : Xhits) delete hist;
        for (auto& hist : Yhits) delete hist;
        for (auto& hist : SumEnergy) delete hist;
        for(int detA = 0; detA < NDet ;  ++detA)
        {
          E1E2_detAdetB[detA].clear();
        }

        HitPatterns.clear();
        Xhits.clear();
        Yhits.clear();
        SumEnergy.clear();
        E1E2_detAdetB.clear();
        //delete E1E2_all;
    }

    void WriteHistograms  () {
        // Open a file with the provided output file name
        TFile outputFile(outputFileName.c_str(), "RECREATE");

        if (!outputFile.IsOpen()) {
            std::cerr << "Error: File could not be opened for writing: " << outputFileName << std::endl;
            return;
        }

        EventMultHist->Write();
        // Create the main groups (directories)
        TDirectory* geometry = outputFile.mkdir("Geometry");
        TDirectory* allEnergy = outputFile.mkdir("All Energy");
        TDirectory* energyCorrelation = outputFile.mkdir("Energy correlation");

        // Create subgroups for "Geometry"
        TDirectory* hitPatternDir = geometry->mkdir("HitPattern");
        TDirectory* xhitsDir = geometry->mkdir("Xhits");
        TDirectory* yhitsDir = geometry->mkdir("Yhits");

        // Create subgroups for "Energy correlation"
        TDirectory* twoByTwoDir = energyCorrelation->mkdir("2by2");



        // Writing "Geometry" histograms
        hitPatternDir->cd();
        for (auto& hist : HitPatterns) hist->Write();

        xhitsDir->cd();
        for (auto& hist : Xhits) hist->Write();

        yhitsDir->cd();
        for (auto& hist : Yhits) hist->Write();

        // Writing "All Energy" histograms
        allEnergy->cd();
        for (auto& hist : SumEnergy) hist->Write();

        // Writing "Energy correlation" histograms
        twoByTwoDir->cd();
        for(int detA = 0; detA < NDet ;  ++detA)
        {
          for(int detB = detA; detB < NDet ;  ++detB)
          {
            E1E2_detAdetB[detA][detB]->Write();
          }
        }

        energyCorrelation->cd();
        E1E2_all->Write();
        E1E2_all_diffDet->Write();


        outputFile.Close();
        std::cout << "Histograms successfully written to " << outputFileName << std::endl;
    }

    bool isDead(int detector, int stripX, int stripY) const {
        bool isXDead = false, isYDead = false;

        // Check if the X strip is dead
        auto itX = deadStripsX.find(detector);
        if (itX != deadStripsX.end()) {
            isXDead = itX->second.find(stripX) != itX->second.end();
        }

        // Check if the Y strip is dead
        auto itY = deadStripsY.find(detector);
        if (itY != deadStripsY.end()) {
            isYDead = itY->second.find(stripY) != itY->second.end();
        }

        // Return true if either X or Y strip is dead
        return isXDead || isYDead;
    }


};

#endif // ANALYSISCONFIG_HH
