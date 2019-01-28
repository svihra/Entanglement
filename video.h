#ifndef VIDEO_H
#define VIDEO_H

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>

#define MAX_HITS  65536      // minimal value for sorting 65536

class Video
{
public:
    Video();

private:
    void LoadFiles(TString folder = "/home/svihra/Data/TPX3/ENT/5min/");
    void RunLoop(TTree *tree, TString name, bool dual = false);

    void Plot(TH2F* map, TCanvas* can, TString name);

private:
    UInt_t      Id_;
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

    TFile* procFile_;
    TFile* entFile_;

    TTree* procTree_;
    TTree* entExclTree_;
    TTree* entTree_;

    ULong64_t timeStart_;
    ULong64_t timeStep_;

    UInt_t nFrames_;

    TString folder_;
};

#endif // VIDEO_H
