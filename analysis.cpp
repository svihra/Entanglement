#include "analysis.h"
#include <algorithm>
#include <iostream>

#include <fstream>

#include "TStopwatch.h"

#include "TString.h"
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <TLegend.h>

#include <TTree.h>
#include <TH1.h>
#include <TH2.h>

#include <TGFileDialog.h>
#include <TRandom.h>

#define MAX_DIFF    156.25//15600 //936        //ns
#define TIME_CORR   50 //ns
#define CUT_TOT         0
#define CUT_SIZE        1
#define TIME         300 // seconds

#define NEVENTS 0

#define TOA_CHECK   100 //ns

// y = 35 - 66
// x = 171 - 200

// y = 168 - 211
// x = 114 - 156


//UInt_t X1_LOW    = 136;
//UInt_t X1_HIGH   = 146;
//UInt_t Y1_LOW    =  97;
//UInt_t Y1_HIGH   = 107;

//UInt_t X2_LOW    = 125;
//UInt_t X2_HIGH   = 135;
//UInt_t Y2_LOW    = 112;
//UInt_t Y2_HIGH   = 122;
UInt_t X1_LOW    = 171;
UInt_t X1_HIGH   = 200;
UInt_t Y1_LOW    =  35;
UInt_t Y1_HIGH   =  66;

UInt_t X2_LOW    = 114;
UInt_t X2_HIGH   = 156;
UInt_t Y2_LOW    = 168;
UInt_t Y2_HIGH   = 211;

//UInt_t X1_LOW    = 101;
//UInt_t X1_HIGH   = 130;
//UInt_t Y1_LOW    =  35;
//UInt_t Y1_HIGH   =  66;

//UInt_t X2_LOW    = 44;
//UInt_t X2_HIGH   = 86;
//UInt_t Y2_LOW    = 168;
//UInt_t Y2_HIGH   = 211;


const char * fileTypes[] = {
        "ROOT files", "*.root",
        0, 0
};

Analysis::Analysis()
{
//    std::cout << "Starting Time Scan: " << std::endl;
//    TimeScan();

//    Randomizer();

    Fitter();
    fileRoot_->Close();
//    if (!Browse())
//    {
//        std::cout << "Failed to load data" << std::endl;
//    }
//    else
//    {
//        for (Int_t file = 0; file < numberOfFiles_; file++)
//        {
//            if (inputFiles_[file].EndsWith(".root"))
//            {
//                std::cout << "Starting Init: " << std::endl;
////                Comparison(inputFiles_[file]);
////                Init(inputFiles_[file]);
////    //            std::cout << "Starting Entanglement for file: " << inputFiles_[file] << std::endl;
////                Processing();
////    //            Entanglement();
//            }
//            else
//            {
//                std::cout << "Invalid file name: " << inputFiles_[file] << std::endl;
//            }
//        }

//        fileRoot_->Close();
//    //    Mapper();
//    //    Correlation();

//    }
}

void Analysis::Init(TString file)
{
    std::cout << "Reading file" << std::endl;
    fileName_ = file;
    pdfName_ = fileName_;
    pdfName_.ReplaceAll(".root","_coin.pdf");

    std::cout << "Reading tree" << std::endl;
    fileRoot_   = new TFile(fileName_, "UPDATE");
    rawTree_ = (TTree *) fileRoot_->Get("proctree");
    std::cout << " - setting branches" << std::endl;
    rawTree_->SetBranchAddress("Size", &Size_);
    rawTree_->SetBranchAddress("Col",  Cols_);
    rawTree_->SetBranchAddress("Row",  Rows_);
    rawTree_->SetBranchAddress("ToT",  ToTs_);
    rawTree_->SetBranchAddress("ToA",  ToAs_);

    Entries_ = rawTree_->GetEntries();

    std::cout << "Create writing file" << std::endl;
    outputName_ = fileName_;
    outputName_.ReplaceAll(".root","_processed.root");
    outputRoot_ = new TFile(outputName_,"RECREATE");
}

void Analysis::Fitter(TString file)
{
    fileRoot_   = new TFile(file, "UPDATE");
    TDirectory* dir = fileRoot_->GetDirectory("entry_0x0_0x0");

    TFile* oldRoot = new TFile("/home/svihra/Data/TPX3/ENT/laser/PolA000PolB000_DiodeCurrent_35,2mAW0028_H11-180711-164524-1_rawtree_processed.root","READ");
    TDirectory* ddir = oldRoot->GetDirectory("entry_0x0_0x0");

    TString ltFile = file;
    ltFile.Replace(ltFile.Last('/'),200,"/LT.csv");

    TCanvas* can = new TCanvas("can","",800,800);
    can->cd();

    TTree* tmpTree = (TTree * ) dir->Get("entTree");
    TTree* ttmpTree = (TTree * ) ddir->Get("entTree");

    Float_t x_min, x_max;
    Float_t y_min, y_max;
    x_min = -12.5;
    x_max = 2512.5;
    y_min = 0;
    y_max = 1;
    tmpTree->Draw("ToT[0]>>tot1(101,-12.5,2512.5)","ID==0","goff");
    tmpTree->Draw("ToT2[0]>>tot2(101,-12.5,2512.5)","ID==0","goff");

    // sigma fit values
    Int_t sig_entries = 100;
    Float_t x[sig_entries];
    Float_t xerr[sig_entries];
    Float_t y[sig_entries];
    Float_t yerr[sig_entries];

    Float_t yO[sig_entries];
    Float_t yerrO[sig_entries];

    // new
    for (int i = 1; i <= sig_entries; i++)
    {
        x[i] = i*25;
        xerr[i] = 12.5;
        TF1* func = new TF1("fit","[0] + [1]*exp((-0.5)*pow((x-[2])/[3],2)) + [4]*exp(-0.5*pow((x-[5])/[6],2))");
        TF1* bcg  = new TF1("bcg","[0] + [1]*exp((-0.5)*pow((x-[2])/[3],2))");
        bcg->SetRange(-150,150);
        func->SetParLimits(0,0,100);
        func->SetParLimits(1,0,1e4);
        func->SetParLimits(2,-5,15);
        func->SetParLimits(3,5,100);
        func->SetParLimits(4,0,1e9);
        func->SetParLimits(5,-5,15);
        func->SetParLimits(6,1,10);
        func->SetParameters(5,10,8,70,20000/(i+2),3,3);
        TString range, name, result;
        range.Form("ID==0 && ToT[0] >= %d && ToT2[0] >= %d",(int) x[i], (int) x[i]);
        name.Form("((ToA[0]-ToA2[0])*25/4096)>>coin%d(65,-50-0.78125,50+0.78125)",(int) i);
        int count = tmpTree->Draw(name,range,"");
        result.Form("coin%d",(int) i);

        TH1D *coin = (TH1D*) gDirectory->Get(result);
        if (count > 200)
        {
            coin->Fit(func);

            y[i] = func->GetParameter(6);
            yerr[i] = func->GetParError(6);

            bcg->SetParameter(0, func->GetParameter(0));
            bcg->SetParameter(1, func->GetParameter(1));
            bcg->SetParameter(2, func->GetParameter(2));
            bcg->SetParameter(3, func->GetParameter(3));
            bcg->Draw("same");

            if (y[i] == 2)
            {
                y[i] = 0;
                yerr[i] = 0;
            }
        }
        else
        {
            y[i] = 0;
            yerr[i] = 0;
        }
        TString ttt = file;
        ttt.Replace(ttt.Last('/'),300,"/fits/"+result+".png");
        can->Print(ttt);
        delete func;
        delete bcg;
    }

    // old
    for (int i = 1; i <= sig_entries; i++)
    {
        TF1* func = new TF1("fit","[0] + [1]*exp((-0.5)*pow((x-[2])/[3],2)) + [4]*exp(-0.5*pow((x-[5])/[6],2))");
        TF1* bcg  = new TF1("bcg","[0] + [1]*exp((-0.5)*pow((x-[2])/[3],2))");
        bcg->SetRange(-150,150);
        func->SetParLimits(0,0,200);
        func->SetParLimits(1,0,1e4);
        func->SetParLimits(2,0,15);
        func->SetParLimits(3,6,200);
        func->SetParLimits(4,0,1e9);
        func->SetParLimits(5,-20,30);
        func->SetParLimits(6,2,20);
        func->SetParameters(20,150,8,100,6000/(0.75*i),5,7);
        TString range, name, result;
        range.Form("ID==0 && ToT[0] >= %d && ToT2[0] >= %d",(int) x[i], (int) x[i]);
        name.Form("((ToA[0]-ToA2[0])*25/4096)>>coin%d(201,-156.25-0.78125,156.25+0.78125)",(int) i);
        int count = ttmpTree->Draw(name,range,"");
        result.Form("coin%d",(int) i);

        TH1D *coin = (TH1D*) gDirectory->Get(result);
        if (count > 200)
        {
            coin->Fit(func);

            yO[i] = func->GetParameter(6);
            yerrO[i] = func->GetParError(6);

            bcg->SetParameter(0, func->GetParameter(0));
            bcg->SetParameter(1, func->GetParameter(1));
            bcg->SetParameter(2, func->GetParameter(2));
            bcg->SetParameter(3, func->GetParameter(3));
            bcg->Draw("same");

            if (y[i] == 2)
            {
                yO[i] = 0;
                yerrO[i] = 0;
            }
        }
        else
        {
            yO[i] = 0;
            yerrO[i] = 0;
        }
        TString ttt = file;
        ttt.Replace(ttt.Last('/'),300,"/fits/old_"+result+".png");
        can->Print(ttt);
        delete func;
        delete bcg;
    }

    // first ToT dist
    TH1D *tot1 = (TH1D*) gDirectory->Get("tot1");
    tot1->SetTitle("");
    tot1->SetStats(kFALSE);
    y_max = tot1->GetMaximum()*1.1;
    tot1->GetXaxis()->SetTitle("ToT [ns]");
    tot1->GetYaxis()->SetTitle("entries [a.u.]");
    tot1->GetYaxis()->SetTitleOffset(1.55);
    tot1->GetYaxis()->SetRangeUser(y_min, y_max);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.20);
    tot1->Draw();

    // second ToT dist
    TH1D *tot2 = (TH1D*) gDirectory->Get("tot2");
    tot2->SetLineColor(kBlack);
    tot2->SetTitle("");
    tot2->SetStats(kFALSE);
    tot2->Draw("same");

    // LT table data
    TGraph* graph = new TGraph(ltFile);
    graph->SetTitle("");

    Float_t rightmax = 1.1*120;
    Float_t scale = y_max/rightmax;
    for (int i=0;i<graph->GetN();i++) graph->GetY()[i] *= -1000*scale;
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);
    graph->Draw("same*");

    TGaxis *axis = new TGaxis(x_max,0, x_max, y_max,0,rightmax,510,"+L");
    axis->SetTitle("dToA [ns]");
    axis->SetTextFont(1);
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();


    Float_t rightmax_sig = 1.1*12;
    Float_t scale_sig = y_max/rightmax_sig;
    for (int i=0;i<sig_entries;i++) {y[i] *= scale_sig; yerr[i] *= scale_sig;}

    TGraph* sigma = new TGraphErrors(sig_entries,x,y,xerr,yerr);
    sigma->SetLineColor(kGreen+3);
    sigma->SetMarkerColor(kGreen+3);
    sigma->SetFillColor(kGreen+3);
    sigma->SetFillStyle(3010);
    sigma->Draw("samep");

    for (int i=0;i<sig_entries;i++) {yO[i] *= scale_sig; yerrO[i] *= scale_sig;}

    TGraph* sigmaO = new TGraphErrors(sig_entries,x,yO,xerr,yerrO);
    sigmaO->SetLineColor(kSpring);
    sigmaO->SetMarkerColor(kSpring);
    sigmaO->SetFillColor(kSpring);
    sigmaO->SetFillStyle(3010);
    sigmaO->Draw("samep");


    TGaxis *axis_sig = new TGaxis(1.15*x_max, 0, 1.15*x_max, y_max, 0, rightmax_sig,510,"+L");
    axis_sig->SetTitle("sigma [ns]");
    axis_sig->SetTextFont(1);
    axis_sig->SetLineColor(kGreen+3);
    axis_sig->SetLabelColor(kGreen+3);
    axis_sig->Draw();

    TLegend* legend = new TLegend(0.5, 0.65, 0.75, 0.85);
    legend->AddEntry(tot1  , "ToT1 distribution", "l");
    legend->AddEntry(tot2  , "ToT2 distribution", "l");
    legend->AddEntry(graph , "ToT correction"   , "l");
    legend->AddEntry(sigmaO, "#sigma entries raw"   , "l");
    legend->AddEntry(sigma , "#sigma entries corr"   , "l");
    legend->Draw("same");

    // saving file
    TString tmpFile = file;
    tmpFile.Replace(tmpFile.Last('.'),200,"_plot.png");
    can->Print(tmpFile);
}

void Analysis::Comparison(TString file)
{
    UInt_t area1[4] = {X1_LOW, X1_HIGH, Y1_LOW, Y1_HIGH};
    UInt_t area2[4] = {X2_LOW, X2_HIGH, Y2_LOW, Y2_HIGH};

    std::cout << "Reading tree" << std::endl;
    fileRoot_   = new TFile(file, "UPDATE");
    rawTree_  = (TTree *) fileRoot_->Get("rawtree");
    procTree_ = (TTree *) fileRoot_->Get("proctree");
    std::cout << " - setting branches" << std::endl;
    rawTree_->SetBranchAddress("Size", &Size_);
    rawTree_->SetBranchAddress("Col",  Cols_);
    rawTree_->SetBranchAddress("Row",  Rows_);
    rawTree_->SetBranchAddress("ToT",  ToTs_);
    rawTree_->SetBranchAddress("ToA",  ToAs_);
    procTree_->SetBranchAddress("Size", &Size2_);
    procTree_->SetBranchAddress("Col",  Cols2_);
    procTree_->SetBranchAddress("Row",  Rows2_);
    procTree_->SetBranchAddress("ToT",  ToTs2_);
    procTree_->SetBranchAddress("ToA",  ToAs2_);

    Entries_ = rawTree_->GetEntries();

    fileRoot_->cd();
    TH1F* subtract    = new TH1F("subtract", "subtract", 301, 1000 -157.03125, 1000 + 157.03125);
    TH2F* subtractMap = new TH2F("subtractMap", "subtractMap", 301, 1000 - 157.03125, 1000+ 157.03125, 400,0,10000);

    Float_t dToA;
    for (int i = 0; i < Entries_; i++)
    {
        if (i % 100000 == 0)
        {
            std::cout << i << " out of " << Entries_ << std::endl;
        }

        rawTree_->GetEntry(i);
        procTree_->GetEntry(i);

        if (PositionCheck(Cols_[0], Rows_[0], area1, ToTs_[0], Size_) || PositionCheck(Cols_[0], Rows_[0], area2, ToTs_[0], Size_))
        {
            if (ToAs_[0] >= ToAs2_[0])
                dToA = ((Float_t) (ToAs_[0] - ToAs2_[0]) * (25.0/4096));
            else
                dToA = ((Float_t) (ToAs2_[0] - ToAs_[0]) * (25.0/4096));

            subtract->Fill(dToA);
            subtractMap->Fill(dToA,ToTs_[0]);
        }
    }

    subtract->Write();
    subtractMap->Write();
    fileRoot_->Write();
}

Bool_t Analysis::PositionCheck(UInt_t &x, UInt_t &y,  UInt_t area[4], UInt_t tot, UInt_t size)
{
    if ( !(tot <= CUT_TOT) && !(size < CUT_SIZE) && x > area[0] && x < area[1] && y > area[2] && y < area[3] )
        return kTRUE;
    else
        return kFALSE;
}

Bool_t Analysis::PositionCheck(UInt_t &x, UInt_t &xLow, UInt_t &xHigh,UInt_t &y,  UInt_t &yLow, UInt_t &yHigh, UInt_t tot, UInt_t size)
{
    if ( !(tot <= CUT_TOT) && !(size < CUT_SIZE) && x > xLow && x < xHigh && y > yLow && y < yHigh )
        return kTRUE;
    else
        return kFALSE;
}

void Analysis::GetWindow(data &out, UInt_t &entry, UInt_t &size, UInt_t cols[], UInt_t rows[], ULong64_t toas[], UInt_t tots[])
{
    out.Sizes[entry] = size;
    out.Cols[entry] = cols[0];
    out.Rows[entry] = rows[0];
    out.ToAs[entry] = toas[0];
    out.ToTs[entry] = tots[0];
}

void Analysis::SaveWindow(data &src, UInt_t &size, UInt_t sizes[],  UInt_t cols[], UInt_t rows[], ULong64_t toas[], UInt_t tots[])
{
    size = src.Entries;
    for (UInt_t entry = 0; entry < src.Entries; entry++)
    {
        sizes[entry] = src.Sizes[entry];
        cols[entry]  = src.Cols[entry];
        rows[entry]  = src.Rows[entry];
        toas[entry]  = src.ToAs[entry];
        tots[entry]  = src.ToTs[entry];
    }
}

bool Analysis::Browse()
{
    for (Int_t file = 0; file < MAX_FILES; file++)
    {
        inputFiles_[file].Clear();
    }
    numberOfFiles_ = 0;

    static TString dir(".");
    TGFileInfo fileInfo;
    fileInfo.SetMultipleSelection(kTRUE);
    fileInfo.fFileTypes = fileTypes;
    fileInfo.fIniDir = StrDup(dir);

    //
    // browsing dialog - multiple files or one file only
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fileInfo);
    if (fileInfo.fFileNamesList)
    {
        TObjString *el;
        TIter next(fileInfo.fFileNamesList);
        while ((el = (TObjString *) next()))
        {
            inputFiles_[numberOfFiles_++] = el->GetString();
        }
    }
    else
        return false;
    inputDir_= "/home/svihra/Data/TPX3/ENT/long/";
    inputDir_ = fileInfo.fIniDir;
    return true;
}

void Analysis::TimeScan()
{
    for (Int_t entry = 0; entry < Entries_; entry++)
    {
        rawTree_->GetEntry(entry);
        Float_t ToA = ((Float_t) (ToAs_[0]))*25.0/4096e9;

        if (ToA > 2.898434 && ToA < 2.898440)
        {
            std::cout << "### Main ToA(ns): "<< ((ToAs_[0]*25.0/4096e9)-2.898434)*1e9 <<", for entry N:" << entry << std::endl;
            std::cout << "### List of all parts in cluster:" << entry << std::endl;
            for (Int_t cluster = 0; cluster < (Int_t) Size_; cluster++)
            {
                std::cout << ((ToAs_[cluster]*25.0/4096e9)-2.898434)*1e9 << ", " << Cols_[cluster] << ", " << Rows_[cluster] << std::endl;
            }
            std::cout << "Size of the cluster: " << Size_ << std::endl;
            std::cout << "###############################################" << std::endl;
        }
    }
}

void Analysis::Randomizer()
{
    TH1F* ent = new TH1F("ent", "Background", 1 + ((2.0*15625)/1.5625), -15625-0.78125, 15625+0.78125);
    TCanvas* canvas = new TCanvas("canvas", "plot", 400, 400);

    Int_t laserTime1,laserTime2;
    Int_t bcgTime1,bcgTime2;

    Int_t timeDiff;

    Int_t actTime;

    TRandom* source = new TRandom();

    for (UInt_t entry = 0; entry < 1e5; entry++)
    {
        actTime    = (Int_t) source->Uniform(5e5);
        bcgTime1   = (Int_t) source->Uniform(5e5) - 5e5;
        bcgTime2   = (Int_t) source->Uniform(5e5) + 5e5;

        timeDiff = actTime + bcgTime1;
        if ( std::abs(bcgTime2 - actTime) < std::abs(timeDiff))
            timeDiff = bcgTime2 - actTime;

        laserTime2 = (Int_t) source->Uniform(5e4) - 5e4;
        if ( std::abs(actTime + laserTime2) < std::abs(timeDiff))
            timeDiff = actTime + laserTime2;

        for (UInt_t shot = 0; shot < 11; shot++)
        {
            laserTime1 = ((Int_t) source->Uniform(5e4)) + ((shot)*5e4);

            if (std::abs(actTime - laserTime1) < std::abs(timeDiff))
                timeDiff = actTime - laserTime1;

            if ( shot != 0 && std::abs(laserTime1 - laserTime2) < std::abs(timeDiff))
                timeDiff = laserTime1 - laserTime2;
            else if (std::abs(laserTime1 + laserTime2) < std::abs(timeDiff))
                timeDiff = laserTime1 + laserTime2;

            laserTime2 = laserTime1;
        }
        ent->Fill(timeDiff);
    }
    canvas->cd();
    ent->Draw();
    canvas->Print("random.png");
    canvas->Close();
    delete source;
}

UInt_t Analysis::FindPairs(UInt_t area[4], Int_t &entry, data &fiber)
{
    ULong64_t diffGlobal = (ULong64_t) (163.84 * MAX_DIFF);
    ULong64_t diffToA = (ULong64_t) (163.84 * MAX_DIFF);
    Bool_t bFound   = kFALSE;
    Bool_t bSmaller = kTRUE;
    UInt_t nextEntry;

    UInt_t count = 0;

    // find pairs
    for (Int_t pair = std::max(entry - ENTRY_LOOP, 0); pair < std::min(entry + ENTRY_LOOP, Entries_); pair++)
    {
        rawTree_->GetEntry(pair);
        if (pair == entry)
        {
            bSmaller = kFALSE;
            continue;
        }

        // area of incoming photons
        if ( PositionCheck(Cols_[0], Rows_[0], area, ToTs_[0], Size_) )
        {
            if ( bSmaller && (ToAs2_[0] - ToAs_[0]) < diffGlobal )
            {
                if ( (ToAs2_[0] - ToAs_[0]) < diffToA )
                {
                    diffToA = ToAs2_[0] - ToAs_[0];
                    nextEntry = pair;
                }
                bFound = kTRUE;
                GetWindow(fiber, count, Size_, Cols_, Rows_, ToAs_, ToTs_);
                count++;

            }
            else if ( !bSmaller && (ToAs_[0] - ToAs2_[0]) < diffGlobal )
            {
                if ( (ToAs_[0] - ToAs2_[0]) < diffToA )
                {
                    diffToA = ToAs_[0] - ToAs2_[0];
                    nextEntry = pair;
                }
                bFound = kTRUE;
                GetWindow(fiber, count, Size_, Cols_, Rows_, ToAs_, ToTs_);
                count++;
            }
        }
    }

    if (bFound)
    {
        fiber.Entries = count;
        return nextEntry;
    }
    else
    {
        fiber.Entries = 0;
    }
    return 0;
}

void Analysis::Processing()
{
    UInt_t mainEntry, nextEntry;
    UInt_t area1[4] = {X1_LOW, X1_HIGH, Y1_LOW, Y1_HIGH};
    UInt_t area2[4] = {X2_LOW, X2_HIGH, Y2_LOW, Y2_HIGH};

    outputRoot_->cd();

    UInt_t id = 0;

    TTree* entTree = new TTree("entTree","entangled data");
    entTree->Branch("ID"  ,       &id,"ID/i");
    entTree->Branch("mE"  ,&mainEntry,"mE/i");
    entTree->Branch("nE"  ,&nextEntry,"nE/i");
    entTree->Branch("Size",    &Size_,"Size/i");
    entTree->Branch("Col" ,     Cols_, "Col[Size]/i");
    entTree->Branch("Row" ,     Rows_, "Row[Size]/i");
    entTree->Branch("ToT" ,     ToTs_, "ToT[Size]/i");
    entTree->Branch("ToA" ,     ToAs_, "ToA[Size]/l");
    entTree->Branch("Size2",   &Size2_,"Size2/i");
    entTree->Branch("Col2" ,    Cols2_, "Col2[Size2]/i");
    entTree->Branch("Row2" ,    Rows2_, "Row2[Size2]/i");
    entTree->Branch("ToT2" ,    ToTs2_, "ToT2[Size2]/i");
    entTree->Branch("ToA2" ,    ToAs2_, "ToA2[Size2]/l");

    ULong64_t ToA;
    UInt_t ToT, Col, Row, Size;
    UInt_t Sizes1[MAX_HITS], Sizes2[MAX_HITS];
    data f1, f2;

    TTree* suppTree = new TTree("suppTree","support data");
    suppTree->Branch("Size"   ,    &Size , "Size/i");
    suppTree->Branch("Col"    ,    &Col  , "Col/i");
    suppTree->Branch("Row"    ,    &Row  , "Row/i");
    suppTree->Branch("ToT"    ,    &ToT  , "ToT/i");
    suppTree->Branch("ToA"    ,    &ToA  , "ToA/l");
    suppTree->Branch("f1"     ,    &Size_, "f1/i");
    suppTree->Branch("f1Size" ,    Sizes1, "f1Size[f1]/i");
    suppTree->Branch("f1Col"  ,     Cols_, "f1Col[f1]/i");
    suppTree->Branch("f1Row"  ,     Rows_, "f1Row[f1]/i");
    suppTree->Branch("f1ToT"  ,     ToTs_, "f1ToT[f1]/i");
    suppTree->Branch("f1ToA"  ,     ToAs_, "f1ToA[f1]/l");
    suppTree->Branch("f2"     ,   &Size2_, "f2/i");
    suppTree->Branch("f2Size" ,    Sizes2, "f2Size[f2]/i");
    suppTree->Branch("f2Col2" ,    Cols2_, "f2Col[f2]/i");
    suppTree->Branch("f2Row2" ,    Rows2_, "f2Row[f2]/i");
    suppTree->Branch("f2ToT2" ,    ToTs2_, "f2ToT[f2]/i");
    suppTree->Branch("f2ToA2" ,    ToAs2_, "f2ToA[f2]/l");


    for (Int_t entry = 0; entry < Entries_; entry++)
    {
        if (NEVENTS != 0 && entry >= NEVENTS)
            break;

        rawTree_->GetEntry(entry);
        if (entry % 1000 == 0) std::cout << "Entry " << entry << " of " << Entries_ << " done!" << std::endl;

        // area of incoming photons
        if ( PositionCheck(Cols_[0], Rows_[0], area1, ToTs_[0], Size_) || PositionCheck(Cols_[0], Rows_[0], area2, ToTs_[0], Size_))
        {
            mainEntry = entry;
            Size2_ = Size_;
            for (UInt_t subentry = 0; subentry < Size_; subentry++)
            {
                Cols2_[subentry] = Cols_[subentry];
                Rows2_[subentry] = Rows_[subentry];
                ToAs2_[subentry] = ToAs_[subentry];
                ToTs2_[subentry] = ToTs_[subentry];
            }

            Size = Size2_;
            Col  = Cols2_[0];
            Row  = Rows2_[0];
            ToA  = ToAs2_[0];
            ToT  = ToTs2_[0];

            if (PositionCheck(Cols_[0], Rows_[0], area1, ToTs_[0], Size_))
            {
                nextEntry = FindPairs(area2,entry, f1);
                if (nextEntry != 0)
                {
                    id = 0;
                    rawTree_->GetEntry(nextEntry);
                    entTree->Fill();
                }

                nextEntry = FindPairs(area1,entry, f2);
                if (nextEntry != 0)
                {
                    id = 1;
                    rawTree_->GetEntry(nextEntry);
                    entTree->Fill();
                }

                if (f1.Entries != 0 || f2.Entries != 0)
                {
                    SaveWindow(f1, Size_, Sizes1, Cols_, Rows_, ToAs_, ToTs_);
                    SaveWindow(f2, Size2_, Sizes2, Cols2_, Rows2_, ToAs2_, ToTs2_);
                    suppTree->Fill();
                }
            }
            else
            {
                nextEntry = FindPairs(area1,entry, f2);

                nextEntry = FindPairs(area2,entry, f1);
                if (nextEntry != 0)
                {
                    id = 2;
                    rawTree_->GetEntry(nextEntry);
                    entTree->Fill();
                }

                if (f1.Entries != 0 || f2.Entries != 0)
                {
                    SaveWindow(f1, Size_, Sizes1, Cols_, Rows_, ToAs_, ToTs_);
                    SaveWindow(f2, Size2_, Sizes2, Cols2_, Rows2_, ToAs2_, ToTs2_);
//                    suppTree->Fill();
                }
            }
        }
    }
    entTree->Write();
    suppTree->Write();
    outputRoot_->Close();
}

void Analysis::Entanglement()
{
    ULong64_t   ToA;

    ULong64_t   diffToA;
    Bool_t      bSmaller;
    Bool_t      bFound;
    Bool_t      bSign;

    UInt_t      colMain, col, colNext;
    UInt_t      rowMain, row, rowNext;
    UInt_t      index;

    Int_t       pairs = 0;

    outputRoot_->cd();



    TH1F* histRate      = new TH1F("histRate", "Rate", 7e7, 0, TIME);
    TH1F* histRateCent  = new TH1F("histRateCent", "Rate of centroids", 7e7, 0, TIME);

    TH1F* histEnt      = new TH1F("histEnt", "Entanglement", 1 + ((2.0*MAX_DIFF)/1.56), -MAX_DIFF-0.78, MAX_DIFF+0.78);
    TH1F* histEntCheck = new TH1F("histEntCheck", "Entanglement check", 1 + ((2.0*MAX_DIFF)/1.56), -MAX_DIFF-0.78, MAX_DIFF+0.78);
    TH1F* histSingle1  = new TH1F("histSingle1", "Fiber 1.", 1 + ((MAX_DIFF)/1.56), -MAX_DIFF-0.78, 0.78);
    TH1F* histSingle2  = new TH1F("histSingle2", "Fiber 2.", 1 + ((MAX_DIFF)/1.56), -0.78, 0.78+MAX_DIFF);

    TH1I* histCol  = new TH1I("histCol", "centroid check cols", 250, 0, 250);
    TH1I* histRow  = new TH1I("histRow", "centroid check rows", 250, 0, 250);

    TH2I* mapColRow = new TH2I("mapColvsRow", "difference comparison", ((X1_HIGH-X1_LOW))+1, -0.5, X1_HIGH-X1_LOW+0.5, ((Y1_HIGH-Y1_LOW))+1, -0.5, Y1_HIGH-Y1_LOW+0.5);
    TH2I* mapColRow2 = new TH2I("mapColvsRow2", "difference comparison", ((X2_HIGH-X2_LOW))+1, -0.5, X2_HIGH-X2_LOW+0.5, ((Y2_HIGH-Y2_LOW))+1, -0.5, Y2_HIGH-Y2_LOW+0.5);

    TH2F* mapDiff = new TH2F("mapDiff", "ToT vs ToA diff", 200, 0, 1, 260, 0, 405.6);

    TH1F* histSpatialX     = new TH1F("histSpatialX", "deltaX",X1_HIGH - X2_LOW,0,X1_HIGH - X2_LOW);
    TH1F* histSpatialY     = new TH1F("histSpatialY", "deltaY",Y2_HIGH - Y1_LOW,0,Y2_HIGH - Y1_LOW);
    TH1F* histSpatialR     = new TH1F("histSpatialR", "deltaR",(Int_t) (std::pow(std::pow(X1_HIGH - X2_LOW,2) + std::pow(Y2_HIGH - Y1_LOW,2),0.5) + 1), 0, (Int_t) (std::pow(std::pow(X1_HIGH - X2_LOW,2) + std::pow(Y2_HIGH - Y1_LOW,2),0.5) + 1));

    TH2I* mapSpatialDiff  = new TH2I("mapSpatialDiff" , "deltaX vs deltaY", X1_HIGH - X2_LOW, 0, X1_HIGH - X2_LOW, Y2_HIGH - Y1_LOW, 0, Y2_HIGH - Y1_LOW);
    TH2I* mapSpatialAdd   = new TH2I("mapSpatialAdd"  , "sumX vs sumY", X1_HIGH + X2_HIGH - X2_LOW, 0, X1_HIGH + X2_HIGH - X2_LOW, Y2_HIGH + Y1_HIGH - Y1_LOW, 0, Y2_HIGH + Y1_HIGH - Y1_LOW);
    TH2I* mapSpatial1     = new TH2I("mapSpatial1"    , "X1 vs Y1", X1_HIGH - X1_LOW, 0, X1_HIGH - X1_LOW, Y1_HIGH - Y1_LOW, 0, Y1_HIGH - Y1_LOW);
    TH2I* mapSpatial2     = new TH2I("mapSpatial2"    , "X2 vs Y2", X2_HIGH - X2_LOW, 0, X2_HIGH - X2_LOW, Y2_HIGH - Y2_LOW, 0, Y2_HIGH - Y2_LOW);

    for (Int_t entry = 0; entry < Entries_; entry++)
    {
        if (NEVENTS != 0 && entry >= NEVENTS)
            break;

        rawTree_->GetEntry(entry);
        if (entry % 1000 == 0) std::cout << "Entry " << entry << " of " << Entries_ << " done!" << std::endl;

        // choose fiber 1.
        // area of incoming photons
        if ( PositionCheck(Cols_[0], X1_LOW, X1_HIGH, Rows_[0], Y1_LOW, Y1_HIGH, ToTs_[0], Size_) )
        {
            // get rates
            histRateCent->Fill(ToAs_[0]*25.0/4096e9);
            for (UInt_t subentry = 0; subentry < Size_; subentry++)
            {
                histRate->Fill(ToAs_[subentry]*25.0/4096e9);
            }

            ToA = ToAs_[0];
            diffToA = (ULong64_t) (163840) * MAX_DIFF;
            bSmaller = kTRUE;
            bFound   = kFALSE;
            bSign    = kTRUE;

            colMain = Cols_[0];
            rowMain = Rows_[0];

            // find pairs from fiber 2.
            for (Int_t entangled = std::max(entry - ENTRY_LOOP,0); entangled < std::min(entry + ENTRY_LOOP, Entries_); entangled++)
                //            for (Int_t entangled = entry - ENTRY_LOOP; entangled < entry + ENTRY_LOOP; entangled++)
            {
                rawTree_->GetEntry(entangled);
                if (entangled == entry)
                    bSmaller = kFALSE;

                // area of incoming photons
                if ( PositionCheck(Cols_[0], X2_LOW, X2_HIGH, Rows_[0], Y2_LOW, Y2_HIGH, ToTs_[0], Size_) )
                {
                    if ( bSmaller && (ToA - ToAs_[0]) < diffToA )
                    {
                        diffToA = ToA - ToAs_[0];
                        bFound = kTRUE;
                        colNext = Cols_[0];
                        rowNext = Rows_[0];
                    }
                    else if ( !bSmaller && (ToAs_[0] - ToA) < diffToA )
                    {
                        diffToA = ToAs_[0] - ToA;
                        bFound = kTRUE;
                        bSign  = kFALSE;
                        colNext = Cols_[0];
                        rowNext = Rows_[0];
                    }
                }
            }
            if (bFound)
            {
                if (bSign)
                {
                    histEnt->Fill( ((Float_t) (diffToA) * 25.0 ) / 4096);
                }
                else
                {
                    histEnt->Fill(-((Float_t) (diffToA) * 25.0 ) / 4096);
                }

                histSpatialX->Fill(std::abs((int)colNext - (int)colMain));
                histSpatialY->Fill(std::abs((int)rowMain - (int)rowNext));
                histSpatialR->Fill(std::pow(std::pow((int)rowMain - (int)rowNext, 2) + std::pow((int)rowMain - (int)rowNext, 2),0.5));

                mapSpatialDiff->Fill( std::abs((int)colNext - (int)colMain), std::abs((int)rowMain - (int)rowNext));
                mapSpatialAdd ->Fill( std::abs((int)colNext + (int)colMain) - X2_LOW, std::abs((int)rowMain + (int)rowNext) - Y1_LOW);
                mapSpatial1   ->Fill( std::abs((int)colMain - X1_LOW), std::abs((int)rowMain - Y1_LOW));
                mapSpatial2   ->Fill( std::abs((int)colNext - X2_LOW), std::abs((int)rowNext - Y2_LOW));

//                std::cout << "deltaY: "<<std::abs((int)rowMain - (int)rowNext) << std::endl;

            }

            // reset values for checking
            diffToA = (ULong64_t) (163840) * MAX_DIFF;
            bSmaller = kTRUE;
            bFound   = kFALSE;
            bSign    = kTRUE;
            index = 0;
            // sanity check (compare fiber 1. to itself)
            //            for (Int_t single = entry - ENTRY_LOOP; single < entry + ENTRY_LOOP; single++)
//            for (Int_t single = std::max(entry - ENTRY_LOOP,0); single < std::min(entry + ENTRY_LOOP, Entries_); single++)
            for (Int_t single = entry; single < std::min(entry + ENTRY_LOOP, Entries_); single++)
            {
                rawTree_->GetEntry(single);
                if (single == entry)
                {
                    bSmaller = kFALSE;
                    continue;
                }

                // area of incoming photons
                if ( single != entry && PositionCheck(Cols_[0], X1_LOW, X1_HIGH, Rows_[0], Y1_LOW, Y1_HIGH, ToTs_[0], Size_) )
                {
                    if ( bSmaller && (ToA - ToAs_[0]) < diffToA )
                    {
                        diffToA = ToA - ToAs_[0];
                        bFound = kTRUE;
                        col = Cols_[0];
                        row = Rows_[0];
                        index = single;
                    }
                    else if ( !bSmaller && (ToAs_[0] - ToA) < diffToA )
                    {
                        diffToA = ToAs_[0] - ToA;
                        bFound = kTRUE;
                        bSign  = kFALSE;
                        col = Cols_[0];
                        row = Rows_[0];
                        index = single;
                    }
                }
            }
            if (bFound)
            {
                if (bSign)
                    histSingle1->Fill((((Float_t) (diffToA) * 25.0 ) / 4096));
                else
                    histSingle1->Fill(-((Float_t) (diffToA) * 25.0 ) / 4096);

                if ((((Float_t) (diffToA) * 25.0 ) / 4096) < TOA_CHECK )
                {
                    if ((std::abs(((int) colMain) - ((int) col))+std::abs(((int) rowMain) - ((int) row))) < 3)
                    {
                        std::cout << "#################################################" << std::endl;
                        std::cout << "Something went wrong for entries: " << std::endl;
                        std::cout << "First: " << entry << ", second:" << index << std::endl;
                        std::cout << "Difference: " << (std::abs(((int) colMain) - ((int) col))+std::abs(((int) rowMain) - ((int) row))) << std::endl;
                        std::cout << "Timediff in ns: " << (((Float_t) (diffToA) * 25.0 ) / 4096) << std::endl;

                        rawTree_->GetEntry(entry);
                        ULong64_t min = ToAs_[0];
                        UInt_t  toot = ToTs_[0];
                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            if (ToAs_[cluster] < min)
                            {
                                min = ToAs_[cluster];
                                toot = ToTs_[cluster];
                            }
                        }

                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            mapDiff->Fill(((Float_t) ToTs_[cluster])/toot,((Float_t) ( ToAs_[cluster] - min)* 25.0 )/4096);
                        }
                        rawTree_->GetEntry(index);
                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            mapDiff->Fill(((Float_t) ToTs_[cluster])/toot,((Float_t) ( ToAs_[cluster] - min)* 25.0 )/4096);
                        }

                    }
                    histCol->Fill(pairs,std::abs(((int) colMain) - ((int) col)));
                    histRow->Fill(pairs,std::abs(((int) rowMain) - ((int) row)));

                    mapColRow->Fill(std::abs(((int) colMain) - ((int) col)), std::abs(((int) rowMain) - ((int) row)));

                    pairs++;
                }
            }
        }
        // choose fiber 2.
        // area of incoming photons
        else if ( PositionCheck(Cols_[0], X2_LOW, X2_HIGH, Rows_[0], Y2_LOW, Y2_HIGH, ToTs_[0], Size_))
        {
            ToA = ToAs_[0];
            diffToA = (ULong64_t) (163840) * MAX_DIFF;
            bSmaller = kTRUE;
            bFound   = kFALSE;
            bSign    = kTRUE;

            colMain = Cols_[0];
            rowMain = Rows_[0];

            // find pairs from fiber 1.
            for (Int_t entangled = std::max(entry - ENTRY_LOOP,0); entangled < std::min(entry + ENTRY_LOOP, Entries_); entangled++)
                //            for (Int_t entangled = entry - ENTRY_LOOP; entangled < entry + ENTRY_LOOP; entangled++)
            {
                rawTree_->GetEntry(entangled);
                if (entangled == entry)
                    bSmaller = kFALSE;

                // area of incoming photons
                if ( PositionCheck(Cols_[0], X1_LOW, X1_HIGH, Rows_[0], Y1_LOW, Y1_HIGH, ToTs_[0], Size_))
                {
                    if ( bSmaller && (ToA - ToAs_[0]) < diffToA )
                    {
                        diffToA = ToA - ToAs_[0];
                        bFound = kTRUE;
                        colNext = Cols_[0];
                        rowNext = Rows_[0];
                    }
                    else if ( !bSmaller && (ToAs_[0] - ToA) < diffToA )
                    {
                        diffToA = ToAs_[0] - ToA;
                        bFound = kTRUE;
                        bSign  = kFALSE;
                        colNext = Cols_[0];
                        rowNext = Rows_[0];
                    }
                }
            }
            if (bFound)
            {
                if (bSign)
                    histEntCheck->Fill( ((Float_t) (diffToA) * 25.0 ) / 4096);
                else
                    histEntCheck->Fill( ((Float_t) (diffToA) * 25.0 ) / 4096);
            }

            diffToA = (ULong64_t) (163840) * MAX_DIFF;
            bSmaller = kTRUE;
            bFound   = kFALSE;
            bSign    = kTRUE;
            index = 0;
            // sanity check (compare fiber 2. to itself)
            //            for (Int_t single = entry - ENTRY_LOOP; single < entry + ENTRY_LOOP; single++)
//            for (Int_t single = std::max(entry - ENTRY_LOOP,0); single < std::min(entry + ENTRY_LOOP, Entries_); single++)
            for (Int_t single = entry; single < std::min(entry + ENTRY_LOOP, Entries_); single++)
            {
                rawTree_->GetEntry(single);
                if (single == entry)
                {
                    bSmaller = kFALSE;
                    continue;
                }

                // area of incoming photons
                if ( single != entry && PositionCheck(Cols_[0], X2_LOW, X2_HIGH, Rows_[0], Y2_LOW, Y2_HIGH, ToTs_[0], Size_) )
                {
                    if ( bSmaller && (ToA - ToAs_[0]) < diffToA )
                    {
                        diffToA = ToA - ToAs_[0];
                        bFound = kTRUE;
                        col = Cols_[0];
                        row = Rows_[0];
                        index = single;
                    }
                    else if ( !bSmaller && (ToAs_[0] - ToA) < diffToA )
                    {
                        diffToA = ToAs_[0] - ToA;
                        bFound = kTRUE;
                        bSign  = kFALSE;
                        col = Cols_[0];
                        row = Rows_[0];
                        index = single;
                    }
                }
            }
            if (bFound)
            {
                if (bSign)
                    histSingle2->Fill(-((Float_t) (diffToA) * 25.0 ) / 4096);
                else
                    histSingle2->Fill( ((Float_t) (diffToA) * 25.0 ) / 4096);

                if ((((Float_t) (diffToA) * 25.0 ) / 4096) < TOA_CHECK )
                {
                    if ((std::abs(((int) colMain) - ((int) col))+std::abs(((int) rowMain) - ((int) row))) < 3)
                    {
                        std::cout << "#################################################" << std::endl;
                        std::cout << "Something went wrong for entries in 2nd: " << std::endl;
                        std::cout << "First: " << entry << ", second:" << index << std::endl;
                        std::cout << "Difference: " << (std::abs(((int) colMain) - ((int) col))+std::abs(((int) rowMain) - ((int) row))) << std::endl;
                        std::cout << "Timediff in ns: " << (((Float_t) (diffToA) * 25.0 ) / 4096) << std::endl;

                        rawTree_->GetEntry(entry);
                        ULong64_t min = ToAs_[0];
                        UInt_t  toot = ToTs_[0];
                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            if (ToAs_[cluster] < min)
                            {
                                min = ToAs_[cluster];
                                toot = ToTs_[cluster];
                            }
                        }

                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            mapDiff->Fill(((Float_t) ToTs_[cluster])/toot,((Float_t) ( ToAs_[cluster] - min)* 25.0 )/4096);
                        }
                        rawTree_->GetEntry(index);
                        for (UInt_t cluster = 0; cluster < Size_; cluster++)
                        {
                            mapDiff->Fill(((Float_t) ToTs_[cluster])/toot,((Float_t) ( ToAs_[cluster] - min)* 25.0 )/4096);
                        }

                    }
                    histCol->Fill(pairs,std::abs(((int) colMain) - ((int) col)));
                    histRow->Fill(pairs,std::abs(((int) rowMain) - ((int) row)));

                    mapColRow2->Fill(std::abs(((int) colMain) - ((int) col)), std::abs(((int) rowMain) - ((int) row)));

                    pairs++;
                }
            }
        }
    }


//    TCanvas* canvas = new TCanvas("Entanglement", pdfName, 1200, 800);
//    canvas->cd();

//    histEnt->Draw();
//    canvas->Print(pdfName);

    std::cout << "Saving files" << std::endl;

    histEnt->Write();
    histEntCheck->Write();
    histSingle1->Write();
    histSingle2->Write();

    histRate->Write();
    histRateCent->Write();

    histCol->Write();
    histRow->Write();

    mapColRow->Write();
    mapColRow2->Write();
    mapDiff->Write();

    histSpatialX->Write();
    histSpatialY->Write();
    histSpatialR->Write();

    mapSpatialDiff->Write();
    mapSpatialAdd ->Write();
    mapSpatial1   ->Write();
    mapSpatial2   ->Write();

    outputRoot_->Close();

//    fileRoot_->Close();
}

void Analysis::Mapper()
{
    TFile* files[MAX_FILES];

    for (Int_t file = 0; file < numberOfFiles_; file++)
    {
        files[file] = new TFile(inputFiles_[file], "READ");
        TString fileName = inputFiles_[file];

        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(0,fileName.First('_')+1);

        TString outputName = inputDir_ + "/processed/map_" + fileName + ".pdf";
        std::cout << "Mapped file " << outputName <<std::endl;

        std::cout << "Reading from file " << inputFiles_[file] <<std::endl;

        UInt_t id;

        TTree* entTree;

        TH2F* map[3][3];
        TDirectory* dir;

        int norm = 0;

        for (UInt_t i = 0; i < 3; i++)
        {
            for (UInt_t j = 0; j < 3; j++)
            {
                TString tmpName;
                tmpName.Form("%dx%d",i,j);
                std::cout << tmpName << std::endl;
//                        tmpName.Form("entry_%dx%d_%dx%d",i,j,x,y);
                map[i][j] = new TH2F(tmpName, tmpName, 3, -0.5, 2.5, 3, -0.5, 2.5);

                for (UInt_t x = 0; x < 3; x++)
                {
                    for (UInt_t y = 0; y < 3; y++)
                    {
                        TString tmpName2;
                        tmpName2.Form("entry_%dx%d_%dx%d",i,j,x,y);

                        std::cout << tmpName2 << std::endl;
                        dir = files[file]->GetDirectory(tmpName2);
                        dir->ls();


                        entTree = (TTree * ) dir->Get("entTree");
//                        entTree->SetBranchAddress("ID",&id);

                        int entries = entTree->Draw("ID","ID==0","goff");

                        if ( entries > norm)
                            norm = entries;

                        map[i][j]->Fill(x,y,entries);
                    }
                }
            }
        }

        TCanvas* canvas = new TCanvas("canvas", fileName, 1200, 1200);
        canvas->Divide(3,3);
        for (UInt_t i = 0; i < 3; i++)
        {
            for (UInt_t j = 0; j < 3; j++)
            {
                canvas->cd(1+ i + (j*3));
                gPad->SetLogz();
                map[i][j]->GetZaxis()->SetRangeUser(1,norm);
                map[i][j]->SetStats(kFALSE);
                map[i][j]->Draw("colz text");
            }
        }
        canvas->Print(outputName);
        canvas->Close();

    }
}

void Analysis::Correlation()
{
//    TString outputName = inputFiles_[0];
//    outputName.ReplaceAll("ent_THL230*","bell.pdf");

    TFile* files[MAX_FILES];

//    Int_t fileX = 0;
//    Int_t fileY = 0;

//    for (Int_t divisor = (Int_t) sqrt(numberOfFiles_); divisor > 0; divisor--)
//    {
//        if (numberOfFiles_ % divisor == 0)
//        {
//            fileX = divisor;
//            fileY = numberOfFiles_ / divisor;

//            std::cout << "Setting plot matrix as X: " << fileX << ", Y: " << fileY << std::endl;
//            break;
//        }
//    }

//    std::cout << "Canvas dimensions X: " << fileX << " Y: " << fileY << std::endl;


    for (Int_t file = 0; file < numberOfFiles_; file++)
    {

        files[file] = new TFile(inputFiles_[file], "READ");
        TString fileName = inputFiles_[file];

        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(fileName.Last('_'),200);
//        fileName.Remove(fileName.Last('_'),200);
        fileName.Remove(0,fileName.First('_')+1);

        TString outputName = inputDir_ + "/processed/" + fileName + ".pdf";
        std::cout << "Correlated file " << outputName <<std::endl;

        TCanvas* canvas = new TCanvas("canvas", fileName, 400, 400);
        std::ofstream csvFile;
        csvFile.open(inputDir_ + "/processed/"+fileName+".csv");
//        fileName.ReplaceAll("ENT_","");

        std::cout << "Reading from file " << inputFiles_[file] <<std::endl;

        UInt_t id, mainEntry, nextEntry;

        TTree* entTree = (TTree * ) files[file]->Get("entTree");
        entTree->SetBranchAddress("ID"  ,       &id);
        entTree->SetBranchAddress("mE"  ,&mainEntry);
        entTree->SetBranchAddress("nE"  ,&nextEntry);
        entTree->SetBranchAddress("Size",    &Size_);
        entTree->SetBranchAddress("Col" ,     Cols_);
        entTree->SetBranchAddress("Row" ,     Rows_);
        entTree->SetBranchAddress("ToT" ,     ToTs_);
        entTree->SetBranchAddress("ToA" ,     ToAs_);
        entTree->SetBranchAddress("Size2",  &Size2_);
        entTree->SetBranchAddress("Col2" ,   Cols2_);
        entTree->SetBranchAddress("Row2" ,   Rows2_);
        entTree->SetBranchAddress("ToT2" ,   ToTs2_);
        entTree->SetBranchAddress("ToA2" ,   ToAs2_);

        Long64_t entries = entTree->GetEntries();

        TH1F* ent = new TH1F("histEnt", "Entangled " + fileName, 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
        TH1F* single1 = new TH1F("histSingle1", "Fiber 1 " + fileName, 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
        TH1F* single2 = new TH1F("histSingle2", "Fiber 2 " + fileName, 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);

        Float_t dToA;
        Bool_t bSign;

        for (Long64_t entry = 0; entry < entries; entry++)
        {
            entTree->GetEntry(entry);

            if (entry % 1000 == 0)
            {
                std::cout << "Entangled entry: " << entry << std::endl;
            }

            if (ToAs2_[0] > ToAs_[0])
            {
                dToA = (Float_t) (ToAs2_[0] - ToAs_[0]);
                bSign = kTRUE;
            }
            else
            {
                dToA = (Float_t) (ToAs_[0] - ToAs2_[0]);
                bSign = kFALSE;
            }

            dToA *= 25.0/4096;

            if (id == 0)
            {
                if (bSign)
                    ent->Fill(dToA);
                else
                    ent->Fill(-dToA);
            }
            else if (id == 1)
            {
                if (bSign)
                    single1->Fill(dToA);
            }
            else if (id == 2)
            {
                if (!bSign)
                    single2->Fill(-dToA);
            }

        }
//        TH1F* ent     = (TH1F *) files[file]->Get("histEnt");
//        TH1F* single1 = (TH1F *) files[file]->Get("histSingle1");
//        TH1F* single2 = (TH1F *) files[file]->Get("histSingle2");

        canvas->cd(file+1);
        ent->SetTitle(fileName);
//        ent->GetXaxis()->SetRangeUser(-500,500);
        ent->SetLineColor(1);
        single1->SetLineColor(2);

        ent->Draw();
        single1->Draw("SAME");
        single2->Draw("SAME");

//        files[file]->Close();
        canvas->Print(outputName);
        canvas->Close();

        for (Int_t bin = 1; bin <= ent->GetNbinsX(); bin++)
        {
            if (bin < ent->GetNbinsX()/2)
                csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single2->GetBinContent(bin) << "\n";
            else if (bin > 1 + ent->GetNbinsX()/2)
                csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single1->GetBinContent(bin) << "\n";
            else
                csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << (single1->GetBinContent(bin)+single2->GetBinContent(bin)) << "\n";        }
        csvFile.close();

    }
//    canvas->Close();
}
