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
#include "TMUSETTData.h"


void NPConvert(const TString runNr, const unsigned short int NDet, const unsigned short int FirstDet, const unsigned short int LastDet){

    if(LastDet - FirstDet + 1 != NDet){
        std::cout << "Error, NDet do not fit with first and last detectors" << std::endl;
        return;
    }
    const unsigned int FIRSTMUSETTDAT = 6;
    const unsigned int NBDET = 4;
    const unsigned int NBSTRIP = 128;
    bool val_X;
    bool val_E;
    unsigned int det_numb;
    unsigned int strip_numb;

    const bool SetMap = false;

    TChain *Chain ;
    TString path = "RootRaw/";
    Chain = new TChain("T");
    std::cout << "Initial file " << path+ "run_0" +runNr+"_label.root" << std::endl;
    Chain->Add(path + "run_0" +runNr+"_label.root");
    ULong64_t GetEntries = Chain->GetEntries();
    if (!GetEntries){
        printf("ERROR in the Chain !! \n");
        return;
    }
    TTreeReader TreeReader(Chain);

    UShort_t MUMU_D;
    TTreeReaderArray<UShort_t> *MUMU_C_;
    //for(unsigned int i = 0; i < NDet; i++){
    //    for(unsigned int j = 0; j < 512; j++){
    //        //std::cout << "test " << i << " " << j << std::endl;
    //        MUMU_C_[i][j] = new TTreeReaderArray<UShort_t>(TreeReader,Form("Event"));
    //    }
    //}
    TTreeReaderValue<std::vector<UShort_t>> *label_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"label");
    TTreeReaderValue<std::vector<UShort_t>> *value_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"value");
    //TTreeReaderValue<std::vector<UShort_t>> *Det_X_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Det_X");
    //TTreeReaderValue<std::vector<UShort_t>> *Det_Y_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Det_Y");
    //TTreeReaderValue<std::vector<UShort_t>> *StripX_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Strip_X");
    //TTreeReaderValue<std::vector<UShort_t>> *StripY_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Strip_Y");
    //TTreeReaderValue<std::vector<UShort_t>> *EnergyY_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Energy_Y");
    //TTreeReaderValue<std::vector<UShort_t>> *EnergyX_= new TTreeReaderValue<std::vector<UShort_t>>(TreeReader,"Energy_X");
    //TTreeReaderValue<std::vector<ULong64_t>> *TimeY_= new TTreeReaderValue<std::vector<ULong64_t>>(TreeReader,"Time_Y");
    //TTreeReaderValue<std::vector<ULong64_t>> *TimeX_= new TTreeReaderValue<std::vector<ULong64_t>>(TreeR//////////eader,"Time_X");
    ////TTreeReaderValue<unsigned long> *Timestamp_= new TTreeReaderValue<unsigned long>(TreeReader,"Timestamp");

    std::vector<UShort_t> label;
    std::vector<UShort_t> value;
    //std::vector<UShort_t> DetX;
    //std::vector<UShort_t> DetY;
    //std::vector<UShort_t> StripX;
    //std::vector<UShort_t> StripY;
    //std::vector<UShort_t> EnergyY;
    //std::vector<UShort_t> EnergyX;
    //std::vector<ULong64_t> TimeY;
    //std::vector<ULong64_t> TimeX;
    //unsigned long Timestamp;




    TMUSETTData* MUMUData = new TMUSETTData;
    TTree* outtree = new TTree("RD","RD");
    //outtree->Branch("TH_UP",&TH_UP,"TH_UP/s");
    //outtree->Branch("TH",&TH,"TH/s");
    //outtree->Branch("TML_UP",&TML_UP,"TML_UP/s");
    //outtree->Branch("TML",&TML,"TML/s");
    //outtree->Branch("DATATRIG",&DATATRIG,"DATATRIG/s");

    outtree->Branch("MUSETT", "TMUSETTData", &MUMUData);

    clock_t start = clock(), current, end;
    while (TreeReader.Next())
    {
        Long64_t cEntry = TreeReader.GetCurrentEntry();;
        if (cEntry%10000 == 0){
            current = clock();
            Double_t Frac = 1.0*cEntry/GetEntries;
            Double_t Time = ((double) (current-start)/CLOCKS_PER_SEC);
            Double_t TimeLeft = Time*(1/Frac - 1.);

            std::cout << "\rEntry : " << cEntry
                 << "/" << GetEntries
                 << " --- "
                 << Form("%.2f",100.*cEntry/GetEntries) <<" %"
                 << " --- "
                 <<Form("%.00f RunEvt/sec",cEntry/Time)
                <<" --- "
               << " Time Left : "<< Form("%d min ",(int)TimeLeft/60)
               << Form("%01.00d sec",(int)TimeLeft%60)
               << std::flush;
        }
        // std::cout << "test 1" << std::endl;
        //DetX = **Det_X_;
        //DetY = **Det_Y_;
        //StripX = **StripX_;
        //StripY = **StripY_;
        //EnergyY = **EnergyY_;
        //EnergyX = **EnergyX_;
        //TimeY = **TimeY_;
        //TimeX = **TimeX_;
        // std::cout << "test 2 " << StripX.size() << " " << EnergyX.size() << " " << TimeX.size() << std::endl;
        // std::cout << "test 2 " << StripY.size() << " " << EnergyY.size() << " " << TimeY.size() << std::endl;
        //unsigned short multX = DetX.size();
        //unsigned short multY = DetY.size();
        // std::cout << "test 3" << std::endl;
        /*for(unsigned short i = 0; i < multX; i++){
            // std::cout << "test 4 " << i << std::endl;
            if(i%2 == 0)
                MUMUData->SetDSSDXE(SetMap,DetX[i], StripX[i], EnergyX[i/2]);
            else
                MUMUData->SetDSSDXT(SetMap,DetX[i], StripX[i], TimeX[i/2]);
        }
        for(unsigned short i = 0; i < multY; i++){
            if(i%2 == 0)
                MUMUData->SetDSSDYE(SetMap,DetY[i], StripY[i], EnergyY[i/2]);
            else
                MUMUData->SetDSSDYT(SetMap,DetY[i], StripY[i], TimeY[i/2]);
        }
                    //    else
                    //        MUMUData->SetDSSDYE(SetMap,i+FirstDet,(j-256)/2 + 1,MUMU_D);
                    //else{
                    //    if(j < 256)
                    //        MUMUData->SetDSSDXT(SetMap,i+FirstDet,(j-1)/2 + 1,MUMU_D);
                    //    else
                    //        MUMUData->SetDSSDYT(SetMap,i+FirstDet,(j-256-1)/2 + 1,MUMU_D);

    */

    label = **label_;
    int TRIG = 0;
    unsigned int k = 0;
    while (k< label.size() && label[k] != 5)
    {
      k++;
      if( k< label.size() && label[k] == 5)
      {
        value = **value_;
        TRIG = value[k];
      }
    }
    MUMUData->SetTRIGGER(TRIG);
    for(unsigned int i = 0; i < label.size(); i++)
        if(label[i] >= FIRSTMUSETTDAT)
            {
                label[i] -= FIRSTMUSETTDAT;
                value = **value_;
                val_X = (label[i]/(2*NBSTRIP)%2 == 0);
                // std::cout << "test " << label[i]/(2*NBSTRIP) << std::endl;
                val_E = (label[i]%2 == 0);
                det_numb = (label[i]/(4*NBSTRIP));
                strip_numb = ((label[i] - det_numb*4*NBSTRIP)/2);
//                std::cout << label[i] << " " << value[i] << std::endl;
//                std::cout << val_X << " " << val_E << std::endl;
//                std::cout << det_numb << " " << strip_numb << std::endl;
                if(val_X){
                    if(val_E){
                        MUMUData->SetDSSDXE(SetMap,det_numb, strip_numb,value[i]);
                    }
                    else{
                        MUMUData->SetDSSDXT(SetMap,det_numb, strip_numb,value[i]);
                    }
                }
                else{
                    if(val_E){
                        MUMUData->SetDSSDYE(SetMap,det_numb, strip_numb - NBSTRIP,value[i]);
                    }
                    else{
                        MUMUData->SetDSSDYT(SetMap,det_numb, strip_numb - NBSTRIP,value[i]);
                    }
                }

            }
    outtree->Fill();
    MUMUData->Clear();
    }
    TFile *fout = new TFile(Form("./RootR/run_0%s.root",runNr.Data()),"RECREATE");
    //TFile *fout = new TFile(Form("./RootR/run_0273_0276.root"),"RECREATE");
    std::cout << ">>File Openned : " << "\n"<< fout->GetName() << "\n" << std::endl;
    outtree->Write();
    fout->Close();

}


void NPConvertSeveral(int RunMin, int RunMax) {
    for (int k = RunMin; k <= RunMax; ++k) {
        std::string filename =  std::to_string(k);
        NPConvert(filename.c_str(), 4, 0, 3);
    }
}
