//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun May 28 09:29:12 2023 by ROOT version 6.26/10
// from TTree AD/Analysed Data Tree
// found on file: RootA/r0514_000a.root
//////////////////////////////////////////////////////////

#ifndef Tree_e843_h
#define Tree_e843_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TMust2Data.h"

#include "TMugastData.h"



class Tree_e843 : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UShort_t> GATCONFMASTER = {fReader, "GATCONFMASTER"};
   TTreeReaderValue<ULong64_t> GATCONFMASTERTS = {fReader, "GATCONFMASTERTS"};
   TTreeReaderValue<Float_t> GATCONFMASTER_C = {fReader, "GATCONFMASTER_C"};
   TTreeReaderValue<UShort_t> DATATRIG1 = {fReader, "DATATRIG1"};
   TTreeReaderValue<ULong64_t> DATATRIG1TS = {fReader, "DATATRIG1TS"};
   TTreeReaderValue<Float_t> DATATRIG1_C = {fReader, "DATATRIG1_C"};
   TTreeReaderValue<UShort_t> DATATRIG_CATS = {fReader, "DATATRIG_CATS"};
   TTreeReaderValue<ULong64_t> DATATRIG_CATSTS = {fReader, "DATATRIG_CATSTS"};
   TTreeReaderValue<Float_t> DATATRIG_CATS_C = {fReader, "DATATRIG_CATS_C"};
   TTreeReaderValue<TMust2Data> MUST2 = {fReader, "MUST2"};
   TTreeReaderValue<TMugastData> Mugast = {fReader, "Mugast"};
   TTreeReaderValue<ULong64_t> MUGAST_TS = {fReader, "MUGAST_TS"};
   TTreeReaderValue<UShort_t> DATATRIG_MG1 = {fReader, "DATATRIG_MG1"};
   TTreeReaderValue<ULong64_t> DATATRIG_MG1TS = {fReader, "DATATRIG_MG1TS"};
   TTreeReaderValue<Float_t> DATATRIG_MG1_C = {fReader, "DATATRIG_MG1_C"};
   TTreeReaderValue<Float_t> CATS1_X = {fReader, "CATS1_X"};
   TTreeReaderValue<Float_t> CATS1_Y = {fReader, "CATS1_Y"};
   TTreeReaderValue<Float_t> CATS1_XWA = {fReader, "CATS1_XWA"};
   TTreeReaderValue<Float_t> CATS1_YWA = {fReader, "CATS1_YWA"};
   TTreeReaderValue<Int_t> CATS1XVM = {fReader, "CATS1XVM"};
   TTreeReaderArray<Float_t> CATS1XV = {fReader, "CATS1XV"};
   TTreeReaderArray<UShort_t> CATS1XVN = {fReader, "CATS1XVN"};
   TTreeReaderValue<Int_t> CATS1YVM = {fReader, "CATS1YVM"};
   TTreeReaderArray<Float_t> CATS1YV = {fReader, "CATS1YV"};
   TTreeReaderArray<UShort_t> CATS1YVN = {fReader, "CATS1YVN"};
   TTreeReaderValue<Float_t> CATS2_X = {fReader, "CATS2_X"};
   TTreeReaderValue<Float_t> CATS2_Y = {fReader, "CATS2_Y"};
   TTreeReaderValue<Float_t> CATS2_XWA = {fReader, "CATS2_XWA"};
   TTreeReaderValue<Float_t> CATS2_YWA = {fReader, "CATS2_YWA"};
   TTreeReaderValue<Int_t> CATS2XVM = {fReader, "CATS2XVM"};
   TTreeReaderArray<Float_t> CATS2XV = {fReader, "CATS2XV"};
   TTreeReaderArray<UShort_t> CATS2XVN = {fReader, "CATS2XVN"};
   TTreeReaderValue<Int_t> CATS2YVM = {fReader, "CATS2YVM"};
   TTreeReaderArray<Float_t> CATS2YV = {fReader, "CATS2YV"};
   TTreeReaderArray<UShort_t> CATS2YVN = {fReader, "CATS2YVN"};
   TTreeReaderValue<Float_t> Xf = {fReader, "Xf"};
   TTreeReaderValue<Float_t> Tf = {fReader, "Tf"};
   TTreeReaderValue<Float_t> Yf = {fReader, "Yf"};
   TTreeReaderValue<Float_t> Pf = {fReader, "Pf"};
   TTreeReaderValue<Float_t> XfWa = {fReader, "XfWa"};
   TTreeReaderValue<Float_t> TfWa = {fReader, "TfWa"};
   TTreeReaderValue<Float_t> YfWa = {fReader, "YfWa"};
   TTreeReaderValue<Float_t> PfWa = {fReader, "PfWa"};
   TTreeReaderArray<UShort_t> PlasticRaw = {fReader, "PlasticRaw"};
   TTreeReaderArray<ULong64_t> PlasticRawTS = {fReader, "PlasticRawTS"};
   TTreeReaderValue<UShort_t> TAC_CATS_PL = {fReader, "TAC_CATS_PL"};
   TTreeReaderValue<ULong64_t> TAC_CATS_PLTS = {fReader, "TAC_CATS_PLTS"};
   TTreeReaderValue<UShort_t> TAC_CATS_HF = {fReader, "TAC_CATS_HF"};
   TTreeReaderValue<ULong64_t> TAC_CATS_HFTS = {fReader, "TAC_CATS_HFTS"};
   TTreeReaderValue<UShort_t> TAC_CATS_EXOGAM = {fReader, "TAC_CATS_EXOGAM"};
   TTreeReaderValue<ULong64_t> TAC_CATS_EXOGAMTS = {fReader, "TAC_CATS_EXOGAMTS"};
   TTreeReaderValue<UShort_t> TAC_MMG_CATS2 = {fReader, "TAC_MMG_CATS2"};
   TTreeReaderValue<ULong64_t> TAC_MMG_CATS2TS = {fReader, "TAC_MMG_CATS2TS"};
   TTreeReaderValue<UShort_t> TAC_MMG_CATS1 = {fReader, "TAC_MMG_CATS1"};
   TTreeReaderValue<ULong64_t> TAC_MMG_CATS1TS = {fReader, "TAC_MMG_CATS1TS"};
   TTreeReaderValue<UShort_t> TAC_MMG_EXOGAM = {fReader, "TAC_MMG_EXOGAM"};
   TTreeReaderValue<ULong64_t> TAC_MMG_EXOGAMTS = {fReader, "TAC_MMG_EXOGAMTS"};
   TTreeReaderValue<UShort_t> TAC_CATS2_CATS1 = {fReader, "TAC_CATS2_CATS1"};
   TTreeReaderValue<ULong64_t> TAC_CATS2_CATS1TS = {fReader, "TAC_CATS2_CATS1TS"};
   TTreeReaderValue<UShort_t> TAC_D4_CATS1 = {fReader, "TAC_D4_CATS1"};
   TTreeReaderValue<ULong64_t> TAC_D4_CATS1TS = {fReader, "TAC_D4_CATS1TS"};
   TTreeReaderValue<UShort_t> TAC_PL_1 = {fReader, "TAC_PL_1"};
   TTreeReaderValue<ULong64_t> TAC_PL_1TS = {fReader, "TAC_PL_1TS"};
   TTreeReaderValue<UShort_t> TAC_PL_2 = {fReader, "TAC_PL_2"};
   TTreeReaderValue<ULong64_t> TAC_PL_2TS = {fReader, "TAC_PL_2TS"};
   TTreeReaderValue<UShort_t> TAC_PL_3 = {fReader, "TAC_PL_3"};
   TTreeReaderValue<ULong64_t> TAC_PL_3TS = {fReader, "TAC_PL_3TS"};
   TTreeReaderValue<UShort_t> TAC_PL_4 = {fReader, "TAC_PL_4"};
   TTreeReaderValue<ULong64_t> TAC_PL_4TS = {fReader, "TAC_PL_4TS"};
   TTreeReaderValue<UShort_t> TAC_PL_5 = {fReader, "TAC_PL_5"};
   TTreeReaderValue<ULong64_t> TAC_PL_5TS = {fReader, "TAC_PL_5TS"};
   TTreeReaderValue<UShort_t> T_DD6_DD4 = {fReader, "T_DD6_DD4"};
   TTreeReaderValue<ULong64_t> T_DD6_DD4TS = {fReader, "T_DD6_DD4TS"};
   TTreeReaderValue<UShort_t> T_CATSD6_1_DD4 = {fReader, "T_CATSD6_1_DD4"};
   TTreeReaderValue<ULong64_t> T_CATSD6_1_DD4TS = {fReader, "T_CATSD6_1_DD4TS"};
   TTreeReaderArray<UShort_t> IC_ZDDRaw = {fReader, "IC_ZDDRaw"};
   TTreeReaderArray<ULong64_t> IC_ZDDRawTS = {fReader, "IC_ZDDRawTS"};
   TTreeReaderValue<Int_t> DCRawM = {fReader, "DCRawM"};
   TTreeReaderArray<UShort_t> DCRawNr = {fReader, "DCRawNr"};
   TTreeReaderArray<UShort_t> DCRaw = {fReader, "DCRaw"};
   TTreeReaderArray<ULong64_t> DCRawTS = {fReader, "DCRawTS"};
   TTreeReaderValue<Float_t> DC_00_C = {fReader, "DC_00_C"};
   TTreeReaderValue<Float_t> DC_01_C = {fReader, "DC_01_C"};
   TTreeReaderValue<Float_t> DC_02_C = {fReader, "DC_02_C"};
   TTreeReaderValue<Float_t> DC_03_C = {fReader, "DC_03_C"};
   TTreeReaderValue<Int_t> EXO_ZDDRawM = {fReader, "EXO_ZDDRawM"};
   TTreeReaderArray<UShort_t> EXO_ZDDRawNr = {fReader, "EXO_ZDDRawNr"};
   TTreeReaderArray<UShort_t> EXO_ZDDRaw = {fReader, "EXO_ZDDRaw"};
   TTreeReaderArray<ULong64_t> EXO_ZDDRawTS = {fReader, "EXO_ZDDRawTS"};
   TTreeReaderValue<Float_t> EXO_ZDD_A_C = {fReader, "EXO_ZDD_A_C"};
   TTreeReaderValue<Float_t> EXO_ZDD_B_C = {fReader, "EXO_ZDD_B_C"};
   TTreeReaderValue<Float_t> EXO_ZDD_C_C = {fReader, "EXO_ZDD_C_C"};
   TTreeReaderValue<Float_t> EXO_ZDD_D_C = {fReader, "EXO_ZDD_D_C"};
   TTreeReaderValue<Int_t> Inner6MVM = {fReader, "Inner6MVM"};
   TTreeReaderArray<Float_t> Inner6MV = {fReader, "Inner6MV"};
   TTreeReaderArray<UShort_t> Inner6MVN = {fReader, "Inner6MVN"};
   TTreeReaderArray<ULong64_t> Inner6MVTS = {fReader, "Inner6MVTS"};
   TTreeReaderValue<Int_t> Inner20MVM = {fReader, "Inner20MVM"};
   TTreeReaderArray<Float_t> Inner20MV = {fReader, "Inner20MV"};
   TTreeReaderArray<UShort_t> Inner20MVN = {fReader, "Inner20MVN"};
   TTreeReaderArray<ULong64_t> Inner20MVTS = {fReader, "Inner20MVTS"};
   TTreeReaderValue<Int_t> DeltaTVM = {fReader, "DeltaTVM"};
   TTreeReaderArray<Float_t> DeltaTV = {fReader, "DeltaTV"};
   TTreeReaderArray<UShort_t> DeltaTVN = {fReader, "DeltaTVN"};
   TTreeReaderArray<ULong64_t> DeltaTVTS = {fReader, "DeltaTVTS"};
   TTreeReaderValue<Int_t> OutersVM = {fReader, "OutersVM"};
   TTreeReaderArray<Float_t> OutersV = {fReader, "OutersV"};
   TTreeReaderArray<UShort_t> OutersVN = {fReader, "OutersVN"};
   TTreeReaderValue<Int_t> BGOVM = {fReader, "BGOVM"};
   TTreeReaderArray<Float_t> BGOV = {fReader, "BGOV"};
   TTreeReaderArray<UShort_t> BGOVN = {fReader, "BGOVN"};
   TTreeReaderValue<Int_t> CSIVM = {fReader, "CSIVM"};
   TTreeReaderArray<Float_t> CSIV = {fReader, "CSIV"};
   TTreeReaderArray<UShort_t> CSIVN = {fReader, "CSIVN"};


   Tree_e843(TTree * /*tree*/ =0) { }
   virtual ~Tree_e843() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Tree_e843,0);

};

#endif

#ifdef Tree_e843_cxx
void Tree_e843::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t Tree_e843::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef Tree_e843_cxx
