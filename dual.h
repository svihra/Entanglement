#ifndef DUAL_H
#define DUAL_H

#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>

#include <iostream>

#define MAX_HITS  65536      // minimal value for sorting 65536
#define MAX_FILES 256
#define DUAL_ENTRY_LOOP 100
#define DUAL_MAX_DIFF    1562.5//15625//156.25//15600 //936        //ns
#define CUT_TOT         0
#define CUT_SIZE        1

// lab I
#define LAB_X_SIZE  11
#define LAB_X_LOW   116
#define LAB_X_HIGH  (LAB_X_LOW + LAB_X_SIZE)
#define LAB_Y_SIZE  11
#define LAB_Y_LOW   108
#define LAB_Y_HIGH  (LAB_Y_LOW + LAB_Y_SIZE)

// tower II
#define TOWER_X_SIZE  8
#define TOWER_X_LOW   107
#define TOWER_X_HIGH  (TOWER_X_LOW + TOWER_X_SIZE)
#define TOWER_Y_SIZE  8
#define TOWER_Y_LOW   89
#define TOWER_Y_HIGH  (TOWER_Y_LOW + TOWER_Y_SIZE)

class Dual
{
public:
    Dual(TString file, TString file2, UInt_t start, Int_t time, Int_t time2, TString tree, UInt_t maxEntries, TString name = "");
    void Process();

private:
    void Init(TString file, TString file2, UInt_t start, Int_t time, Int_t time2, TString tree, UInt_t maxEntries, TString name);

    Bool_t PositionCheck(UInt_t area[4]);
    Bool_t PositionCheck2(UInt_t area[4]);

    Double_t FindDelta(ULong64_t &diff);
    Long64_t FindPairs(UInt_t area[4], Long64_t &entry);
    bool FindEntry(Long64_t &entry2);
    void FindTrig(Long64_t &entry2);
    void ScanEntry(Long64_t &entry, Long64_t &entry2);
    void PrintCsv();

private:
    Long64_t maxEntries_;
    TString inputName_;
    TString inputName2_;

    TFile* fileRoot_;
    TTree* tree_;
    TTree* timeTree_;

    TFile* fileRoot2_;
    TTree* tree2_;
    TTree* timeTree2_;

    TFile* outputRoot_;
    TString outputName_;
    TTree* entTree_;

    TCanvas* can_;
    TLegend* leg_;
    Int_t colors_[5] = {kBlack, kAzure+2, kOrange+1, kRed+1, kGreen+3 };

    Long64_t entry_;
    Long64_t entry2_;

    UInt_t      id_;
    Long64_t    mainEntry_;
    Long64_t    nextEntry_;
    ULong64_t   ToAdiff_;
    ULong64_t   ToAzero_;
    ULong64_t   ToAzero2_;

    Long64_t    Entries_;
    UInt_t      Size_;
    UInt_t      Cols_[MAX_HITS];
    UInt_t      Rows_[MAX_HITS];
    ULong64_t   ToAs_[MAX_HITS];
    ULong64_t   ToATrigs_[MAX_HITS];
    UInt_t      ToTs_[MAX_HITS];

    UInt_t      TrigStart_;

    UInt_t      TrigId_;
    UInt_t      TrigIdNext_;
    ULong64_t   TrigWalk_;
    ULong64_t   TrigTime_;
    ULong64_t   TrigTimeNext_;
    ULong64_t   TrigDiff_;

    Long64_t    Entries2_;
    UInt_t      Size2_;
    UInt_t      Cols2_[MAX_HITS];
    UInt_t      Rows2_[MAX_HITS];
    ULong64_t   ToAs2_[MAX_HITS];
    ULong64_t   ToATrigs2_[MAX_HITS];
    UInt_t      ToTs2_[MAX_HITS];

    UInt_t      TrigId2_;
    UInt_t      TrigIdNext2_;
    ULong64_t   TrigWalk2_;
    ULong64_t   TrigTime2_;
    ULong64_t   TrigTimeNext2_;
    ULong64_t   TrigDiff2_;


    TString     inputFiles_[MAX_FILES];
    TString     inputDir_;

    TString     inputFiles2_[MAX_FILES];
    TString     inputDir2_;

    Int_t       numberOfFiles_;

    UInt_t area1All_[4] = {LAB_X_LOW, LAB_X_HIGH, LAB_Y_LOW, LAB_Y_HIGH};
    UInt_t area2All_[4] = {TOWER_X_LOW, TOWER_X_HIGH, TOWER_Y_LOW, TOWER_Y_HIGH};

};

#endif // DUAL_H
