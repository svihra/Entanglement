#include "video.h"
#include <iostream>

#define MAX_FRAMES   50    // n of frames
#define FRAME_STEP   10    // in ms

Video::Video()
{
    nFrames_ = MAX_FRAMES;
    timeStep_ = 1e6 * FRAME_STEP * 4096.0 / 25.0;

    LoadFiles();

    procTree_->GetEntry(0);
    timeStart_ = ToAs_[0];
    Id_ = 0;

//    RunLoop(procTree_, "movie/data_");
    RunLoop(entTree_ , "movie/entangled_", true);
    RunLoop(entExclTree_ , "movie/entangledExcl_", true);
}

void Video::LoadFiles(TString folder)
{
    folder_ = folder;
    TString fileName = "ENT_A080_B090_W0028_H11-180702-162628-1";
    procFile_ = new TFile(folder + fileName + ".root", "READ");

    procTree_ = (TTree *) procFile_->Get("proctree");
    std::cout << " - setting branches" << std::endl;
    procTree_->SetBranchAddress("Size", &Size_);
    procTree_->SetBranchAddress("Col",  Cols_);
    procTree_->SetBranchAddress("Row",  Rows_);
    procTree_->SetBranchAddress("ToT",  ToTs_);
    procTree_->SetBranchAddress("ToA",  ToAs_);

    entFile_ = new TFile(folder + fileName + "_proctree_processed.root", "READ");
    std::cout << " - entFile_ branches" << std::endl;

    entTree_ = (TTree * ) entFile_->Get("entTree");
    entTree_->SetBranchAddress("ID"   ,      &Id_);
    entTree_->SetBranchAddress("Size" ,    &Size_);
    entTree_->SetBranchAddress("Col"  ,     Cols_);
    entTree_->SetBranchAddress("Row"  ,     Rows_);
    entTree_->SetBranchAddress("ToT"  ,     ToTs_);
    entTree_->SetBranchAddress("ToA"  ,     ToAs_);
    entTree_->SetBranchAddress("Size2",   &Size2_);
    entTree_->SetBranchAddress("Col2" ,    Cols2_);
    entTree_->SetBranchAddress("Row2" ,    Rows2_);
    entTree_->SetBranchAddress("ToT2" ,    ToTs2_);
    entTree_->SetBranchAddress("ToA2" ,    ToAs2_);

    TDirectory* dir = entFile_->GetDirectory("entry_0x0_0x0");
    entExclTree_ = (TTree * ) dir->Get("entTree");
    entExclTree_->SetBranchAddress("ID"   ,      &Id_);
    entExclTree_->SetBranchAddress("Size" ,    &Size_);
    entExclTree_->SetBranchAddress("Col"  ,     Cols_);
    entExclTree_->SetBranchAddress("Row"  ,     Rows_);
    entExclTree_->SetBranchAddress("ToT"  ,     ToTs_);
    entExclTree_->SetBranchAddress("ToA"  ,     ToAs_);
    entExclTree_->SetBranchAddress("Size2",   &Size2_);
    entExclTree_->SetBranchAddress("Col2" ,    Cols2_);
    entExclTree_->SetBranchAddress("Row2" ,    Rows2_);
    entExclTree_->SetBranchAddress("ToT2" ,    ToTs2_);
    entExclTree_->SetBranchAddress("ToA2" ,    ToAs2_);

    std::cout << "done" << std::endl;
}

void Video::Plot(TH2F *map, TCanvas *can, TString name)
{
    can->cd();
    map->SetStats(kFALSE);
    map->Draw("colz");
    can->Print(name);
    can->Close();
}

void Video::RunLoop(TTree* tree, TString name, bool dual)
{
    UInt_t startPos = 0;
    TH2F* totInt     = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
    TH2F* toaInt     = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
    TH2F* entriesInt = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
    TCanvas* canTotInt     = new TCanvas("totInt" + name,"",512,512);
    TCanvas* canToaInt     = new TCanvas("toaInt" + name,"",512,512);
    TCanvas* canEntriesInt = new TCanvas("entInt" + name,"",512,512);

    for (UInt_t frame = 0; frame < nFrames_; frame++)
    {
        std::cout << "Frame: " << frame << "/" << nFrames_ << std::endl;
        TH2F* tot      = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
        TH2F* toa      = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
        TH2F* entries    = new TH2F("","",256,-0.5,255.5,256,-0.5,255.5);
        TCanvas* canTot     = new TCanvas("tot" + std::to_string(frame) + name, "",512,512);
        TCanvas* canToa     = new TCanvas("toa" + std::to_string(frame) + name, "",512,512);
        TCanvas* canEntries = new TCanvas("ent" + std::to_string(frame) + name, "",512,512);

        for (UInt_t entry = startPos; entry < tree->GetEntries(); entry++)
        {
            tree->GetEntry(entry);

            //
            // stop loop if further away than frame size
            if ( (ToAs_[0] - timeStart_) > ((ULong64_t) frame * timeStep_) )
            {
                startPos = entry;
                break;
            }
            else
            {
                if (Id_ == 0)
                {
                    for (UInt_t part = 0; part < Size_; part++)
                    {
                        tot       ->Fill(Rows_[part], Cols_[part], ToTs_[part]);
                        totInt    ->Fill(Rows_[part], Cols_[part], ToTs_[part]);
                        toa       ->Fill(Rows_[part], Cols_[part], (Float_t) (ToAs_[part] - timeStart_ - ((ULong64_t) frame * timeStep_))*25.0/(4096.0*1E6));
                        toaInt    ->Fill(Rows_[part], Cols_[part], (Float_t) (ToAs_[part] - timeStart_ - ((ULong64_t) frame * timeStep_))*25.0/(4096.0*16));
                        entries   ->Fill(Rows_[part], Cols_[part]             );
                        entriesInt->Fill(Rows_[part], Cols_[part]             );
                    }
                    if (dual)
                    {
                        for (UInt_t part = 0; part < Size2_; part++)
                        {
                            tot       ->Fill(Rows2_[part], Cols2_[part], ToTs2_[part]);
                            totInt    ->Fill(Rows2_[part], Cols2_[part], ToTs2_[part]);
                            toa       ->Fill(Rows2_[part], Cols2_[part], (Float_t) (ToAs2_[part] - timeStart_ - ((ULong64_t) frame * timeStep_))*25.0/(4096.0*1E6));
                            toaInt    ->Fill(Rows2_[part], Cols2_[part], (Float_t) (ToAs2_[part] - timeStart_ - ((ULong64_t) frame * timeStep_))*25.0/(4096.0*1E6));
                            entries   ->Fill(Rows2_[part], Cols2_[part]              );
                            entriesInt->Fill(Rows2_[part], Cols2_[part]              );
                        }
                    }
                }
            }
        }

        Plot(totInt    , canTotInt    , folder_ + name + "ToTInt_"     + std::to_string(frame) + ".png");
        Plot(toaInt    , canToaInt    , folder_ + name + "ToAInt_"     + std::to_string(frame) + ".png");
        Plot(entriesInt, canEntriesInt, folder_ + name + "EntriesInt_" + std::to_string(frame) + ".png");
        Plot(tot       , canTot       , folder_ + name + "ToT_"        + std::to_string(frame) + ".png");
        Plot(toa       , canToa       , folder_ + name + "ToA_"        + std::to_string(frame) + ".png");
        Plot(entries   , canEntries   , folder_ + name + "Entries_"    + std::to_string(frame) + ".png");

        delete canTot     ;
        delete canToa     ;
        delete canEntries ;
    }

    delete canTotInt     ;
    delete canToaInt     ;
    delete canEntriesInt ;
}
