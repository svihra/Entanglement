#include "entangled.h"

#include <iostream>
#include <fstream>

#include <TTree.h>
#include <TThread.h>
#include <TH1.h>
#include <TCanvas.h>

Entangled::Entangled(TString fileName, TString tree, UInt_t maxEntries, Int_t startEntryPart, Int_t parts)
{
    if (fileName.EndsWith(".root") && !fileName.EndsWith("_processed.root"))
    {
        std::cout << "Starting init: " << std::endl;
        Init(fileName, tree, maxEntries, startEntryPart, parts);

        std::cout << "Starting entanglement processing: " << std::endl;
        Process();
        outputRoot_->Close();
    }
    else
    {
        std::cerr << "Wrong input file " << fileName << std::endl;
    }

}

Entangled::~Entangled()
{
    for (Int_t x1 = 0; x1 < X1_CUT; x1++)
    {
        for (Int_t y1 = 0; y1 < Y1_CUT; y1++)
        {
            for (Int_t x2 = 0; x2 < X2_CUT; x2++)
            {
                for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
                {
                    delete entTree_[x1][y1][x2][y2];
                }
            }
        }
    }
}

void Entangled::Init(TString file, TString tree, UInt_t maxEntries, Int_t startEntryPart, Int_t parts)
{
    std::cout << "Setting max entries as " << maxEntries << std::endl;
    maxEntries_ = maxEntries;

    std::cout << "Reading file" << std::endl;
    inputName_ = file;

    std::cout << "Reading tree" << std::endl;
    fileRoot_   = new TFile(inputName_, "READ");
    std::cout << "Getting tree" << std::endl;
    rawTree_ = (TTree *) fileRoot_->Get(tree);
    std::cout << " - setting branches" << std::endl;
    rawTree_->SetBranchAddress("Size", &Size_);
    rawTree_->SetBranchAddress("Col",  Cols_);
    rawTree_->SetBranchAddress("Row",  Rows_);
    rawTree_->SetBranchAddress("ToT",  ToTs_);
    rawTree_->SetBranchAddress("ToA",  ToAs_);

    Entries_ = rawTree_->GetEntries();

    std::cout << "Create writing file" << std::endl;
    outputName_ = inputName_;

    Part_ = startEntryPart;
    Parts_ = parts;
    if (Part_ == -1)
    {
        outputName_.ReplaceAll(".root","_"+tree+"_processed.root");
    }
    else
    {
        TString tmpOutput;
        tmpOutput.Form("_%d_%d_processed.root", Part_, Parts_-1);
        outputName_.ReplaceAll(".root","_"+tree+tmpOutput);
    }
    outputRoot_ = new TFile(outputName_,"RECREATE");

    id_         = 0;
    mainEntry_  = 0;
    nextEntry_  = 0;

    for (Int_t x1 = 0; x1 < X1_CUT; x1++)
    {
        for (Int_t y1 = 0; y1 < Y1_CUT; y1++)
        {
            for (Int_t i = 0; i < 2; i++)
            {
                area1_[x1][y1][i]   = X1_LOW + ((i+x1)*(X1_SIZE/X1_CUT));
                area1_[x1][y1][i+2] = Y1_LOW + ((i+y1)*(Y1_SIZE/Y1_CUT));
            }
        }
    }

    for (Int_t x2 = 0; x2 < X2_CUT; x2++)
    {
        for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
        {
            for (Int_t i = 0; i < 2; i++)
            {
                area2_[x2][y2][i]   = X2_LOW + ((i+x2)*(X2_SIZE/X2_CUT));
                area2_[x2][y2][i+2] = Y2_LOW + ((i+y2)*(Y2_SIZE/Y2_CUT));
            }
        }
    }

    for (Int_t x1 = 0; x1 < X1_CUT; x1++)
    {
        for (Int_t y1 = 0; y1 < Y1_CUT; y1++)
        {
            for (Int_t x2 = 0; x2 < X2_CUT; x2++)
            {
                for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
                {
                    std::cout <<"Areas: " << std::endl;
                    std::cout << "Area1: " << area1_[x1][y1][0] << "," << area1_[x1][y1][1] << "x" << area1_[x1][y1][2] << "," << area1_[x1][y1][3] << std::endl;
                    std::cout << "Area2: " << area2_[x2][y2][0] << "," << area2_[x2][y2][1] << "x" << area2_[x2][y2][2] << "," << area2_[x2][y2][3] << std::endl;

                    TString tmpName;
                    tmpName.Form("entry_%dx%d_%dx%d",x1,y1,x2,y2);
                    combinationDir_[x1][y1][x2][y2] = outputRoot_->mkdir(tmpName);
                    combinationDir_[x1][y1][x2][y2]->cd();

                    TTree* tmpTree = new TTree("entTree","entangled data for sections");
                    tmpTree->Branch("ID"   ,       &id_,"ID/i");
                    tmpTree->Branch("mE"   ,&mainEntry_,"mE/i");
                    tmpTree->Branch("nE"   ,&nextEntry_,"nE/i");
                    tmpTree->Branch("Size" ,     &Size_,"Size/i");
                    tmpTree->Branch("Col"  ,      Cols_,"Col[Size]/i");
                    tmpTree->Branch("Row"  ,      Rows_,"Row[Size]/i");
                    tmpTree->Branch("ToT"  ,      ToTs_,"ToT[Size]/i");
                    tmpTree->Branch("ToA"  ,      ToAs_,"ToA[Size]/l");
                    tmpTree->Branch("Size2",    &Size2_,"Size2/i");
                    tmpTree->Branch("Col2" ,     Cols2_,"Col2[Size2]/i");
                    tmpTree->Branch("Row2" ,     Rows2_,"Row2[Size2]/i");
                    tmpTree->Branch("ToT2" ,     ToTs2_,"ToT2[Size2]/i");
                    tmpTree->Branch("ToA2" ,     ToAs2_,"ToA2[Size2]/l");

                    entTree_[x1][y1][x2][y2] = tmpTree;
                }
            }
        }
    }
    outputRoot_->cd();
    entTreeAll_ = new TTree("entTree","entangled data for sections");
    entTreeAll_->Branch("ID"   ,       &id_,"ID/i");
    entTreeAll_->Branch("mE"   ,&mainEntry_,"mE/i");
    entTreeAll_->Branch("nE"   ,&nextEntry_,"nE/i");
    entTreeAll_->Branch("Size" ,     &Size_,"Size/i");
    entTreeAll_->Branch("Col"  ,      Cols_,"Col[Size]/i");
    entTreeAll_->Branch("Row"  ,      Rows_,"Row[Size]/i");
    entTreeAll_->Branch("ToT"  ,      ToTs_,"ToT[Size]/i");
    entTreeAll_->Branch("ToA"  ,      ToAs_,"ToA[Size]/l");
    entTreeAll_->Branch("Size2",    &Size2_,"Size2/i");
    entTreeAll_->Branch("Col2" ,     Cols2_,"Col2[Size2]/i");
    entTreeAll_->Branch("Row2" ,     Rows2_,"Row2[Size2]/i");
    entTreeAll_->Branch("ToT2" ,     ToTs2_,"ToT2[Size2]/i");
    entTreeAll_->Branch("ToA2" ,     ToAs2_,"ToA2[Size2]/l");
}

Bool_t Entangled::PositionCheck(UInt_t area[4])
{
    if ( !(ToTs_[0] <= CUT_TOT) && !(Size_ < CUT_SIZE) && Cols_[0] >= area[0] && Cols_[0] < area[1] && Rows_[0] >= area[2] && Rows_[0] < area[3] )
        return kTRUE;
    else
        return kFALSE;
}

UInt_t Entangled::FindPairs(UInt_t area[4], Int_t &entry, bool inverse)
{
    ULong64_t diffToA = (ULong64_t) (163.84 * MAX_DIFF);
    Bool_t bFound   = kFALSE;
    UInt_t nextEntry;

    // find pairs
    for (Int_t pair = std::max(entry - ENTRY_LOOP, 0); pair < std::min(entry + ENTRY_LOOP, Entries_); pair++)
    {
        rawTree_->GetEntry(pair);
        if (pair == entry)
        {
//            bSmaller = kFALSE;
            continue;
        }

        // area of incoming photons
        if ( !inverse && PositionCheck(area) )
        {
            if ( ToAs2_[0] > ToAs_[0] && (ToAs2_[0] - ToAs_[0]) < diffToA )
            {
                diffToA = ToAs2_[0] - ToAs_[0];
                nextEntry = pair;
                bFound = kTRUE;
            }
            else if ( ToAs2_[0] < ToAs_[0] && (ToAs_[0] - ToAs2_[0]) < diffToA )
            {
                diffToA = ToAs_[0] - ToAs2_[0];
                nextEntry = pair;
                bFound = kTRUE;
            }
        }
        else if ( inverse) // added to find pair from same area not within the area
        {
            if ( ToAs2_[0] > ToAs_[0] && (ToAs2_[0] - ToAs_[0]) < diffToA )
            {
                diffToA = ToAs2_[0] - ToAs_[0];
                nextEntry = pair;
                bFound = kTRUE;
            }
            else if ( ToAs2_[0] < ToAs_[0] && (ToAs_[0] - ToAs2_[0]) < diffToA )
            {
                diffToA = ToAs_[0] - ToAs2_[0];
                nextEntry = pair;
                bFound = kTRUE;
            }
        }
    }

    if (bFound)
    {
        return nextEntry;
    }
    return 0;
}

void Entangled::ScanEntry(Int_t &entry)
{
    rawTree_->GetEntry(entry);
    if (entry % 1000 == 0) std::cout << "Entry " << entry << " of " << Entries_ << " done!" << std::endl;

    // area of incoming photons
    if ( PositionCheck(area1All_) || PositionCheck(area2All_))
    {
        mainEntry_ = entry;
        Size2_ = Size_;
        for (UInt_t subentry = 0; subentry < Size_; subentry++)
        {
            Cols2_[subentry] = Cols_[subentry];
            Rows2_[subentry] = Rows_[subentry];
            ToAs2_[subentry] = ToAs_[subentry];
            ToTs2_[subentry] = ToTs_[subentry];
        }

        Int_t globalX1 = (Cols2_[0] - X1_LOW) / (X1_SIZE/X1_CUT);
        Int_t globalY1 = (Rows2_[0] - Y1_LOW) / (Y1_SIZE/Y1_CUT);
        Int_t globalX2 = (Cols2_[0] - X2_LOW) / (X2_SIZE/X2_CUT);
        Int_t globalY2 = (Rows2_[0] - Y2_LOW) / (Y2_SIZE/Y2_CUT);

        if (PositionCheck(area1All_))
        {
            if ((globalX1 >= 0 && globalX1 < X1_CUT) && (globalY1 >= 0 && globalY1 < Y1_CUT))
            {
                for (Int_t x2 = 0; x2 < X2_CUT; x2++)
                {
                    for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
                    {
                        nextEntry_ = FindPairs(area2_[x2][y2],entry);
                        if (nextEntry_ != 0)
                        {
                            id_ = 0;
                            rawTree_->GetEntry(nextEntry_);
                            entTree_[globalX1][globalY1][x2][y2]->Fill();
                        }
                    }
                }
            }

            nextEntry_ = FindPairs(area2All_,entry);
            if (nextEntry_ != 0)
            {
                id_ = 0;
                rawTree_->GetEntry(nextEntry_);
                entTreeAll_->Fill();
            }

            if ((globalX1 >= 0 && globalX1 < X1_CUT) && (globalY1 >= 0 && globalY1 < Y1_CUT))
            {
                for (Int_t x1 = 0; x1 < X1_CUT; x1++)
                {
                    for (Int_t y1 = 0; y1 < Y1_CUT; y1++)
                    {
                        nextEntry_ = FindPairs(area1_[x1][y1],entry);
                        if (nextEntry_ != 0)
                        {
                            id_ = 1;
                            rawTree_->GetEntry(nextEntry_);
                            entTree_[globalX1][globalY1][x1][y1]->Fill();
                        }
                    }
                }
            }

            nextEntry_ = FindPairs(area1All_,entry);
            if (nextEntry_ != 0)
            {
                id_ = 1;
                rawTree_->GetEntry(nextEntry_);
                entTreeAll_->Fill();
            }

        }
        else
        {
            if ((globalX2 >= 0 && globalX2 < X2_CUT) && (globalY2 >= 0 && globalY2 < Y2_CUT))
            {
                for (Int_t x2 = 0; x2 < X2_CUT; x2++)
                {
                    for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
                    {
                        nextEntry_ = FindPairs(area2_[x2][y2],entry);
                        if (nextEntry_ != 0)
                        {
                            id_ = 2;
                            rawTree_->GetEntry(nextEntry_);
                            entTree_[globalX2][globalY2][x2][y2]->Fill();
                        }
                    }
                }
            }

            nextEntry_ = FindPairs(area2All_,entry);
            if (nextEntry_ != 0)
            {
                id_ = 2;
                rawTree_->GetEntry(nextEntry_);
                entTreeAll_->Fill();
            }

        }
    }
}

void Entangled::PrintCsv()
{
    outputRoot_->cd();
    for (Int_t x1 = 0; x1 < X1_CUT; x1++)
    {
        for (Int_t y1 = 0; y1 < Y1_CUT; y1++)
        {
            for (Int_t x2 = 0; x2 < X2_CUT; x2++)
            {
                for (Int_t y2 = 0; y2 < Y2_CUT; y2++)
                {
                    TString csvName = outputName_;
                    TString pdfName = outputName_;

                    csvName.Remove(csvName.Last('.'),200);
                    csvName.Append(Form("_%dx%d_%dx%d.csv",x1,y1,x2,y2));
                    pdfName.Remove(pdfName.Last('.'),200);
                    pdfName.Append(Form("_%dx%d_%dx%d.pdf",x1,y1,x2,y2));

                    TCanvas* canvas = new TCanvas(Form("canvas_%dx%d_%dx%d",x1,y1,x2,y2), Form("canvas_%dx%d_%dx%d",x1,y1,x2,y2), 400, 400);

                    TH1F* ent = new TH1F(Form("histEnt_%dx%d_%dx%d",x1,y1,x2,y2), Form("Entangled %dx%d_%dx%d",x1,y1,x2,y2), 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
                    TH1F* single1 = new TH1F(Form("histSingle1_%dx%d_%dx%d",x1,y1,x2,y2), Form("Fiber 1 %dx%d_%dx%d",x1,y1,x2,y2), 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
                    TH1F* single2 = new TH1F(Form("histSingle2_%dx%d_%dx%d",x1,y1,x2,y2), Form("Fiber 2 %dx%d_%dx%d",x1,y1,x2,y2), 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);

                    Float_t dToA;
                    Bool_t bSign;

                    for (Int_t entry = 0; entry < entTree_[x1][y1][x2][y2]->GetEntries(); entry++)
                    {
                        entTree_[x1][y1][x2][y2]->GetEntry(entry);
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

                        if (id_ == 0)
                        {
                            if (bSign)
                                ent->Fill(dToA);
                            else
                                ent->Fill(-dToA);
                        }
                        else if (id_ == 1)
                        {
                            if (bSign)
                                single1->Fill(dToA);
                        }
                        else if (id_ == 2)
                        {
                            if (!bSign)
                                single2->Fill(-dToA);
                        }

                    }
                    ent->SetLineColor(1);
                    single1->SetLineColor(2);

                    ent->Draw();
                    single1->Draw("SAME");
                    single2->Draw("SAME");

                    canvas->Print(pdfName);
                    canvas->Close();

                    std::ofstream csvFile;
                    csvFile.open(csvName);

                    for (Int_t bin = 1; bin <= ent->GetNbinsX(); bin++)
                    {
                        if (bin < ent->GetNbinsX()/2)
                            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single2->GetBinContent(bin) << "\n";
                        else if (bin > 1 + ent->GetNbinsX()/2)
                            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single1->GetBinContent(bin) << "\n";
                        else
                            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << (single1->GetBinContent(bin)+single2->GetBinContent(bin)) << "\n";
                    }
                    csvFile.close();
                }
            }
        }
    }

    // general for whole areas
    TH1F* ent = new TH1F("histEnt", "Entangled", 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
    TH1F* single1 = new TH1F("histSingle1", "Fiber 1", 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);
    TH1F* single2 = new TH1F("histSingle2", "Fiber 2", 1 + ((2.0*MAX_DIFF)/1.5625), -MAX_DIFF-0.78125, MAX_DIFF+0.78125);

    for (Int_t entry = 0; entry < entTreeAll_->GetEntries(); entry++)
    {
        Float_t dToA;
        Bool_t bSign;

        entTreeAll_->GetEntry(entry);
        if (entry % 1000 == 0)
        {
            std::cout << "Entangled entry All: " << entry << std::endl;
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

        if (id_ == 0)
        {
            if (bSign)
                ent->Fill(dToA);
            else
                ent->Fill(-dToA);
        }
        else if (id_ == 1)
        {
            if (bSign)
                single1->Fill(dToA);
        }
        else if (id_ == 2)
        {
            if (!bSign)
                single2->Fill(-dToA);
        }
    }

    TString csvName = outputName_;
    csvName.Remove(csvName.Last('.'),200);
    csvName.Append(".csv");
    std::ofstream csvFile;
    csvFile.open(csvName);
    for (Int_t bin = 1; bin <= ent->GetNbinsX(); bin++)
    {
        if (bin < ent->GetNbinsX()/2)
            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single2->GetBinContent(bin) << "\n";
        else if (bin > 1 + ent->GetNbinsX()/2)
            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << single1->GetBinContent(bin) << "\n";
        else
            csvFile << ent->GetBinCenter(bin) << "," << ent->GetBinContent(bin) << "," << (single1->GetBinContent(bin)+single2->GetBinContent(bin)) << "\n";
    }
    csvFile.close();

    outputRoot_->Write();
}

void Entangled::Process()
{
    outputRoot_->cd();

    Int_t startEntry = 0;
    Int_t endEntry = Entries_;

    if (Part_ != -1)
    {
        startEntry = (Int_t) ((((Double_t) Part_)/Parts_)*Entries_);
        endEntry = (Int_t) ((((Double_t) (Part_+1))/Parts_)*Entries_);
    }

    std::cout << "Part: " << Part_ << " Starting at: " << startEntry << ", finishing at: " << endEntry << std::endl;

    for (Int_t entry = startEntry; entry < endEntry; entry++)
    {
        if (maxEntries_ != 0 && entry >= maxEntries_)
            break;
        ScanEntry(entry);
    }

    PrintCsv();
}
