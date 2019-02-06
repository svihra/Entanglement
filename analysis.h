#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <TString.h>
#include <TStopwatch.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TColor.h>

#include <TGFrame.h>

#include <deque>

#define MAX_HITS  65536      // minimal value for sorting 65536
#define MAX_FILES 256
#define ENTRY_LOOP 150

struct data
{
    UInt_t      Entries;
    UInt_t      Size;
    UInt_t      Col;
    UInt_t      Row;
    UInt_t      ToT;
    ULong64_t   ToA;
    UInt_t      Sizes[ENTRY_LOOP*2];
    UInt_t      Cols[ENTRY_LOOP*2];
    UInt_t      Rows[ENTRY_LOOP*2];
    ULong64_t   ToAs[ENTRY_LOOP*2];
    UInt_t      ToTs[ENTRY_LOOP*2];
};

struct LTnames
{
    TString time;
    TString name;
};

class Analysis : public TGMainFrame
{
public:
    Analysis();

    void Processing();
    void Entanglement();
    void Correlation();
    void Mapper();
    void TimeScan();
    void Comparison(TString file);
    void Spatial(TString file = "/home/svihra/Data/TPX3/ENT/long/spatial.txt");
    void Fitter(TString file = "/home/svihra/Data/TPX3/ENT/laser/PolA000PolB000_DiodeCurrent_35,2mAW0028_H11-180711-164524-1_proctree_processed.root");
    void Plotter(TString file, Int_t index);

private:
    void Init(TString file = "/home/svihra/Documents/Timepix3/Data_Acquired/ENT/ent_THL230_gain00_polA000B000_60s_W0028_H11-171114-175414-1.root");
    bool Browse();
    Int_t Degrees(int x, int y);
    Bool_t PositionCheck(UInt_t &x, UInt_t &xLow, UInt_t &xHigh, UInt_t &y, UInt_t &yLow, UInt_t &yHigh, UInt_t tot = 0, UInt_t size = 0);
    Bool_t PositionCheck(UInt_t &x, UInt_t &y,  UInt_t area[4], UInt_t tot, UInt_t size);
    UInt_t FindPairs(UInt_t area[4], Long64_t &entry, data &fiber);
    void Randomizer();

    void GetWindow(data &out, UInt_t &entry, UInt_t &size, UInt_t cols[MAX_HITS], UInt_t rows[MAX_HITS], ULong64_t toas[MAX_HITS], UInt_t tots[MAX_HITS]);
    void SaveWindow(data &src, UInt_t &size, UInt_t sizes[MAX_HITS], UInt_t cols[MAX_HITS], UInt_t rows[MAX_HITS], ULong64_t toas[MAX_HITS], UInt_t tots[MAX_HITS]);

private:


    TString fileName_, pdfName_, outputName_;

    TFile* fileRoot_;
    TTree* rawTree_;
    TTree* procTree_;

    TFile* outputRoot_;

    TCanvas* can_;
    TLegend* leg_;
    Color_t colors_[4] = {kAzure+2, kOrange+1, kRed+1, kGreen+3 };

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

    Long64_t    Entries_;

    TString     inputFiles_[MAX_FILES];
    TString     inputDir_;
    Int_t       numberOfFiles_;

};

#endif // ANALYSIS_H
