#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h"
#include "TCollection.h"
#include "TSpectrum.h"
#include "TGraphErrors.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "TMUSETTPhysics.h"


void Moulinette(const TString runNr){


    TChain *Chain ;
    TString path = "../../Outputs/Analysis/";
    Chain = new TChain("PhysicsTree");
    std::cout << "Initial file " << path+runNr+".root" << std::endl;
    Chain->Add(path+runNr+"*.root");
    ULong64_t GetEntries = Chain->GetEntries();
    if (!GetEntries){
        printf("ERROR in the Chain !! \n");
        return;
    }
    TTreeReader TreeReader(Chain);
    
    
    TTreeReaderValue<TMUSETTPhysics> *MUSETT_= new TTreeReaderValue<TMUSETTPhysics>(TreeReader,"MUSETT");

    TMUSETTPhysics MUSETT;

    clock_t start = clock(), current, end;
    while (TreeReader.Next())
    {
        //Long64_t cEntry = TreeReader.GetCurrentEntry();;
        //if (cEntry%10000 == 0){
        //    current = clock();
        //    Double_t Frac = 1.0*cEntry/GetEntries;
        //    Double_t Time = ((double) (current-start)/CLOCKS_PER_SEC);
        //    Double_t TimeLeft = Time*(1/Frac - 1.);

        //    std::cout << "\rEntry : " << cEntry
        //         << "/" << GetEntries
        //         << " --- "
        //         << Form("%.2f",100.*cEntry/GetEntries) <<" %"
        //         << " --- "
        //         <<Form("%.00f RunEvt/sec",cEntry/Time)
        //        <<" --- "
        //       << " Time Left : "<< Form("%d min ",(int)TimeLeft/60)
        //       << Form("%01.00d sec",(int)TimeLeft%60)
        //       << std::flush;
        //}

    
    MUSETT = **MUSETT_;
    if(MUSETT.DSSD_E.size() == 2){
        if((MUSETT.DSSD_E[0] > 9 && MUSETT.DSSD_E[1] < 7) ||(MUSETT.DSSD_E[1] > 9 && MUSETT.DSSD_E[0] < 7)  ){
            std::cout << "Energy alpha 1 : " << MUSETT.DSSD_E[0] << std::endl;
            std::cout << "Energy alpha 2 : " << MUSETT.DSSD_E[1] << std::endl;
            std::cout << "alpha 1 + alpha 2 : " << MUSETT.DSSD_E[1]+MUSETT.DSSD_E[0] << std::endl;
            std::cout << "Angle : " << MUSETT.RelativeAngle << std::endl;
            std::cout << std::endl;
        }
    }
    }

}