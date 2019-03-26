#ifndef ENTANGLED_H
#define ENTANGLED_H

#include <TString.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

#include <TSystemDirectory.h>
#include <TChain.h>

#define MAX_HITS  65536      // minimal value for sorting 65536
#define MAX_FILES 256
#define ENTRY_LOOP 150
#define MAX_DIFF    156.25//15600 //936        //ns
#define CUT_TOT         0
#define CUT_SIZE        1

// first
#define X1_SIZE  20
#define X1_LOW   187
#define X1_HIGH  (X1_LOW + X1_SIZE)
#define Y1_SIZE  20
#define Y1_LOW   87
#define Y1_HIGH  (Y1_LOW + Y1_SIZE)

#define X1_SIZE_WHOLE  20
#define X1_LOW_WHOLE   187
#define X1_HIGH_WHOLE  (X1_LOW_WHOLE + X1_SIZE_WHOLE)
#define Y1_SIZE_WHOLE  20
#define Y1_LOW_WHOLE   87
#define Y1_HIGH_WHOLE  (Y1_LOW_WHOLE + Y1_SIZE_WHOLE)

// second
#define X2_SIZE  20
#define X2_LOW   115
#define X2_HIGH  (X2_LOW + X2_SIZE)
#define Y2_SIZE  20
#define Y2_LOW   80
#define Y2_HIGH  (Y2_LOW + Y2_SIZE)

#define X2_SIZE_WHOLE  20
#define X2_LOW_WHOLE   115
#define X2_HIGH_WHOLE  (X2_LOW_WHOLE + X2_SIZE_WHOLE)
#define Y2_SIZE_WHOLE  20
#define Y2_LOW_WHOLE   80
#define Y2_HIGH_WHOLE  (Y2_LOW_WHOLE + Y2_SIZE_WHOLE)

#define X1_CUT  1
#define Y1_CUT  1
#define X2_CUT  1
#define Y2_CUT  1

class Entangled
{
public:
    Entangled(TString fileName, TString tree = "rawtree", UInt_t maxEntries = 0, Int_t startEntryPart = -1, Int_t parts = 25);
    ~Entangled();
    void Process();

private:
    Bool_t Init(TString file, TString tree, UInt_t maxEntries, Int_t startEntryPart, Int_t parts);
    Bool_t AddFiles(TSystemFile* dir, TChain* chainDat, TChain* chainTime);

    Bool_t PositionCheck(UInt_t area[4]);
    Long64_t FindPairs(UInt_t area[4], Long64_t &entry, bool inverse = false);
    void ScanEntry(Long64_t &entry);
    void PrintCsv();

private:
    TChain* treeChain_;
    TChain* timeChain_;

    TSystemFile* dir_;

    TString inputName_, outputName_;
    Int_t nThreads_;
    UInt_t maxEntries_;

    TFile* fileRoot_;
//    TTree* rawTree_;

    TFile* outputRoot_;
    TDirectory* combinationDir_[X1_CUT][Y1_CUT][X2_CUT][Y2_CUT];
    TTree* entTree_[X1_CUT][Y1_CUT][X2_CUT][Y2_CUT];
    TTree* entTreeAll_;

    UInt_t      id_;
    Long64_t      mainEntry_;
    Long64_t      nextEntry_;

    UInt_t      Size_;
    UInt_t      Cols_[MAX_HITS];
    UInt_t      Rows_[MAX_HITS];
    ULong64_t   ToAs_[MAX_HITS];
    UInt_t      ToTs_[MAX_HITS];

    UInt_t      Size2_;
    UInt_t      Cols2_[MAX_HITS];
    UInt_t      Rows2_[MAX_HITS];
    ULong64_t   ToAs2_[MAX_HITS];
    UInt_t      ToTs2_[MAX_HITS];

    Long64_t       Entries_;
    Int_t       Part_;
    Int_t       Parts_;

    TString     inputFiles_[MAX_FILES];
    TString     inputDir_;
    Int_t       numberOfFiles_;

    UInt_t area1All_[4] = {X1_LOW_WHOLE, X1_HIGH_WHOLE, Y1_LOW_WHOLE, Y1_HIGH_WHOLE};
    UInt_t area2All_[4] = {X2_LOW_WHOLE, X2_HIGH_WHOLE, Y2_LOW_WHOLE, Y2_HIGH_WHOLE};
    UInt_t area1_[X1_CUT][Y1_CUT][4];// = {X1_LOW, X1_HIGH, Y1_LOW, Y1_HIGH};
    UInt_t area2_[X2_CUT][Y2_CUT][4];// = {X2_LOW, X2_HIGH, Y2_LOW, Y2_HIGH};

};

#endif // ENTANGLED_H
