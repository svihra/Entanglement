#include "dual.h"


// lab is first
// tower is second
Dual::Dual(TString file, TString file2, Float_t diff, Int_t time, Int_t time2, TString tree, UInt_t maxEntries)
{
    if (file.EndsWith(".root")  && !file.EndsWith("_processed.root")
     && file2.EndsWith(".root") && !file2.EndsWith("_processed.root"))
    {
        std::cout << "Starting init: " << std::endl;
        Init(file, file2, diff, time, time2, tree, maxEntries);

        std::cout << "Starting entanglement processing: " << std::endl;
        Process();
        outputRoot_->Write();
        outputRoot_->Close();
    }
    else
    {
        std::cerr << "Wrong input files: " << file << " / " << file2 << std::endl;
    }
}

void Dual::Init(TString file, TString file2, Float_t diff, Int_t time, Int_t time2, TString tree, UInt_t maxEntries)
{
    std::cout << "Setting time difference as " << diff << "us" << std::endl;
    ToAdiff_ = static_cast<ULong64_t>(163840 * diff); // converted to clocks

    std::cout << "Setting max entries as " << maxEntries << std::endl;
    maxEntries_ = maxEntries;

    std::cout << "Reading files: " << std::endl;
    std::cout << file << std::endl;
    std::cout << file2 << std::endl;
    inputName_ = file;
    inputName2_ = file2;

    std::cout << "Reading tree1" << std::endl;
    fileRoot_   = new TFile(inputName_, "READ");
    std::cout << "Getting tree1" << std::endl;
    tree_ = reinterpret_cast<TTree *>(fileRoot_->Get(tree));
    std::cout << " - setting branches" << std::endl;
    tree_->SetBranchAddress("Size"    , &Size_);
    tree_->SetBranchAddress("Col"     ,  Cols_);
    tree_->SetBranchAddress("Row"     ,  Rows_);
    tree_->SetBranchAddress("ToT"     ,  ToTs_);
    tree_->SetBranchAddress("ToA"     ,  ToAs_);
    tree_->SetBranchAddress("ToATrig" ,  ToATrigs_);    // l for long unsigned
    tree_->SetBranchAddress("TrigCntr", &TrigId_);
    tree_->SetBranchAddress("TrigTime", &TrigTime_);
    tree_->SetBranchAddress("TrigTimeNext", &TrigTimeNext_);

    Entries_ = tree_->GetEntries();

    std::cout << "Reading tree2" << std::endl;
    fileRoot2_   = new TFile(inputName2_, "READ");
    std::cout << "Getting tree2" << std::endl;
    tree2_ = reinterpret_cast<TTree *>(fileRoot2_->Get(tree));
    std::cout << " - setting branches" << std::endl;
    tree2_->SetBranchAddress("Size", &Size2_);
    tree2_->SetBranchAddress("Col",  Cols2_);
    tree2_->SetBranchAddress("Row",  Rows2_);
    tree2_->SetBranchAddress("ToT",  ToTs2_);
    tree2_->SetBranchAddress("ToA",  ToAs2_);
    tree2_->SetBranchAddress("ToATrig" ,  ToATrigs2_);    // l for long unsigned
    tree2_->SetBranchAddress("TrigCntr", &TrigId2_);
    tree2_->SetBranchAddress("TrigTime", &TrigTime2_);
    tree2_->SetBranchAddress("TrigTimeNext", &TrigTimeNext2_);

    Entries2_ = tree2_->GetEntries();

    tree_->GetEntry(0);
    tree2_->GetEntry(0);

    std::cout << "Time offset 1: " << time << "s, 2: " << time2 << "s" << std::endl;

    ToAzero_  = ToAs_[0]  + static_cast<ULong64_t>(163840000000 * static_cast<ULong64_t>(time));
    ToAzero2_ = ToAs2_[0] + static_cast<ULong64_t>(163840000000 * static_cast<ULong64_t>(time2));

    for (entry_ = 0; entry_ < Entries_; entry_++)
    {
        tree_->GetEntry(entry_);
        if (TrigTime_ > ToAzero_)
        {
            ToAzero_ = TrigTime_;
            TrigIdNext_ = TrigId_ + 1;
            std::cout << "Trig entry: " << TrigId_ << std::endl;
            break;
        }
    }

    for (entry2_ = 0; entry2_ < Entries2_; entry2_++)
    {
        tree2_->GetEntry(entry2_);
        if (TrigTime2_ > ToAzero2_)
        {
            ToAzero2_ = TrigTime2_;
            TrigIdNext2_ = TrigId2_ + 1;
            std::cout << "Trig entry2: " << TrigId2_ << std::endl;
            break;
        }
    }

    std::cout << "Zero time " << ToAzero_*25.0/4096000000 << "us" << std::endl;
    std::cout << "Zero time2 " << ToAzero2_*25.0/4096000000 << "us" << std::endl;

    std::cout << "Create writing file" << std::endl;
    outputName_ = inputName_;
    outputName_.ReplaceAll(".root","_" + std::to_string(time) + "_" + std::to_string(time2) + "_"+tree+"_processed.root");
    outputRoot_ = new TFile(outputName_,"RECREATE");

    outputRoot_->cd();
    entTree_ = new TTree("entTree", "entangled data, version 2");
    entTree_->Branch("ID"   ,       &id_,"ID/i");
    entTree_->Branch("mE"   ,&mainEntry_,"mE/i");
    entTree_->Branch("nE"   ,&nextEntry_,"nE/i");

    entTree_->Branch("Size" ,     &Size_,"Size/i");
    entTree_->Branch("Col"  ,      Cols_,"Col[Size]/i");
    entTree_->Branch("Row"  ,      Rows_,"Row[Size]/i");
    entTree_->Branch("ToT"  ,      ToTs_,"ToT[Size]/i");
    entTree_->Branch("ToA"  ,      ToAs_,"ToA[Size]/l");
    entTree_->Branch("ToAz" ,  &ToAzero_,"ToAz/l");
    entTree_->Branch("ToATrig" ,  ToATrigs_,"ToATrig[Size]/l" );
    entTree_->Branch("TrigCntr", &TrigId_  ,"&TrigCntr/i");
    entTree_->Branch("TrigTime", &TrigTime_,"&TrigTime/l");
    entTree_->Branch("TrigDiff", &TrigDiff_,"&TrigDiff/l");

    entTree_->Branch("Size2",    &Size2_,"Size2/i");
    entTree_->Branch("Col2" ,     Cols2_,"Col2[Size2]/i");
    entTree_->Branch("Row2" ,     Rows2_,"Row2[Size2]/i");
    entTree_->Branch("ToT2" ,     ToTs2_,"ToT2[Size2]/i");
    entTree_->Branch("ToA2" ,     ToAs2_,"ToA2[Size2]/l");
    entTree_->Branch("ToAz2", &ToAzero2_,"ToAz2/l");
    entTree_->Branch("ToATrig2" ,  ToATrigs2_,"ToATrig2[Size2]/l" );
    entTree_->Branch("TrigCntr2", &TrigId2_  ,"&TrigCntr2/i");
    entTree_->Branch("TrigTime2", &TrigTime2_,"&TrigTime2/l");
    entTree_->Branch("TrigDiff2", &TrigDiff2_,"&TrigDiff2/l");
}

void Dual::Process()
{
    outputRoot_->cd();

    Long64_t startEntry = entry_;
    Long64_t endEntry = Entries_;
    UInt_t    startTrig = TrigId_;

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Starting at: entry1 " << entry_ << ", entry2 " << entry2_ << std::endl;
    for (; entry_ < endEntry; entry_++)
    {
        if (maxEntries_ != 0 && entry_ >= startEntry + maxEntries_)
        {
            std::cout << "reached max entries" << std::endl;
            break;
        }
        if ((TrigId_ - startTrig) > 30)
        {
            std::cout << "reached 30s from start" << std::endl;
            break;
        }
        ScanEntry(entry_, entry2_);
    }
}

// --------------------------------------------------------------
// ------------------- ENTANGLEMENT SCAN ------------------------

Double_t Dual::FindDelta(ULong64_t &diff)
{
    // 16387e7 is 1 second converted to tpx3 clocks
    return static_cast<Double_t>(16384e7)/static_cast<Double_t>(diff);
}

Bool_t Dual::PositionCheck(UInt_t area[4])
{
    if ( !(ToTs_[0] <= CUT_TOT) && !(Size_ < CUT_SIZE) && Cols_[0] >= area[0] && Cols_[0] < area[1] && Rows_[0] >= area[2] && Rows_[0] < area[3] )
        return kTRUE;
    else
        return kFALSE;
}

Bool_t Dual::PositionCheck2(UInt_t area[4])
{
    if ( !(ToTs2_[0] <= CUT_TOT) && !(Size2_ < CUT_SIZE) && Cols2_[0] >= area[0] && Cols2_[0] < area[1] && Rows2_[0] >= area[2] && Rows2_[0] < area[3] )
        return kTRUE;
    else
        return kFALSE;
}

Long64_t Dual::FindPairs(UInt_t area[4], Long64_t &entry)
{
    ULong64_t diffToA = static_cast<ULong64_t>(163.84 * DUAL_MAX_DIFF);
    Bool_t bFound   = kFALSE;
    Long64_t nextEntry = 0;

    // find pairs
    for (Long64_t pair = std::max(entry - DUAL_ENTRY_LOOP, static_cast<Long64_t>(0)); pair < std::min(entry + DUAL_ENTRY_LOOP, Entries2_); pair++)
    {
        tree2_->GetEntry(pair);

        // area of incoming photons
        if ( ( TrigId2_ == (TrigIdNext2_ - 1)) && PositionCheck2(area))
        {
            if ( static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) + ToAdiff_ >= static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_))
            && ( static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) + ToAdiff_ - static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_))) < diffToA )
            {
                diffToA = static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) + ToAdiff_ - static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_));
                nextEntry = pair;
                bFound = kTRUE;
            }
            else if ( static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) + ToAdiff_ < static_cast<ULong64_t>(ToATrigs_[0] *FindDelta(TrigDiff_ ))
                 && ( static_cast<ULong64_t>(ToATrigs_[0] *FindDelta(TrigDiff_ )) - ToAdiff_ - static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_))) < diffToA )
            {
                diffToA = static_cast<ULong64_t>(ToATrigs_[0] *FindDelta(TrigDiff_ )) - ToAdiff_ - static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_));
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

bool Dual::FindEntry(Long64_t &entry2)
{
    tree2_->GetEntry(entry2);
    TrigDiff2_ = TrigTimeNext2_ - TrigTime2_;

    while(static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) < static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_))
       && TrigIdNext2_ != TrigId2_
       && entry2 < Entries2_ )
    {
        tree2_->GetEntry(entry2++);
    }

    if (TrigIdNext2_ == TrigId2_)
        return false;
    else
        return true;
}

void Dual::FindTrig(Long64_t &entry2)
{
    tree2_->GetEntry(entry2);

    while(TrigIdNext2_ != TrigId2_ && entry2 < Entries2_)
    {
        tree2_->GetEntry(entry2++);
    }

    TrigIdNext2_++;
}

void Dual::ScanEntry(Long64_t &entry, Long64_t &entry2)
{
    tree_->GetEntry(entry);
    TrigDiff_ = TrigTimeNext_ - TrigTime_;

    if (TrigIdNext_ == TrigId_)
    {
        TrigIdNext_++;
        FindTrig(entry2);
    }

    if (entry % 1000 == 0)
    {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Entry1 " << entry << " of " << Entries_ << " done!" << std::endl;
        std::cout << "Entry2 " << entry2 << " of " << Entries2_ << " done!" << std::endl;
    }

    // area of incoming photons + entries within checking range of triggers are ignored
    if ( PositionCheck(area1All_)
      && static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) > static_cast<ULong64_t>(163.84 * DUAL_MAX_DIFF)
      && static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) < static_cast<ULong64_t>(163.84 * (1e9 - DUAL_MAX_DIFF)))
    {
        mainEntry_ = entry;

        if (FindEntry(entry2))
        {
            // find entry from other file and pass
            nextEntry_ = FindPairs(area2All_, entry2);
            if (nextEntry_ != 0)
            {
                id_ = 0;
                tree2_->GetEntry(nextEntry_);
                entTree_->Fill();
            }
        }
    }
}
