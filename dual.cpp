#include "dual.h"


// lab is first
// tower is second
Dual::Dual(TString file, TString file2, UInt_t start, Int_t time, Int_t time2, TString tree, UInt_t maxEntries, TString name)
{
    dir_  = new TSystemFile(file, file);
    dir2_ = new TSystemFile(file2, file2);
    if ((dir_->IsDirectory() || (file.EndsWith(".root")  && !file.EndsWith("_processed.root")))
     && (dir2_->IsDirectory() || (file2.EndsWith(".root") && !file2.EndsWith("_processed.root"))))
    {
        std::cout << "Starting init: " << std::endl;
        if (Init(file, file2, start, time, time2, tree, maxEntries, name))
        {
            std::cout << "Starting entanglement processing: " << std::endl;
            Process();
            outputRoot_->Write();
            outputRoot_->Close();
        }
        else
        {
            std::cout << "Something wrong with reading files/dirs" << std::endl;
        }
    }
    else
    {
        std::cerr << "Wrong input files: " << file << " / " << file2 << std::endl;
    }
}

Bool_t Dual::Init(TString file, TString file2, UInt_t start, Int_t time, Int_t time2, TString tree, UInt_t maxEntries, TString name)
{
    std::cout << "Setting starting second as " << start << "s" << std::endl;
    TrigStart_ = start;

    ToAdiff_ = static_cast<ULong64_t>(163840 * 2.06); // converted to clocks - measured diff 2.05us

    std::cout << "Setting max entries as " << maxEntries << std::endl;
    maxEntries_ = maxEntries;

    std::cout << "Reading dirs/files: " << std::endl;
    std::cout << file << std::endl;
    std::cout << file2 << std::endl;

    inputName_ = file;
    inputName2_ = file2;

    treeChain_  = new TChain(tree);
    treeChain2_ = new TChain(tree);

    timeChain_  = new TChain("timetree");
    timeChain2_ = new TChain("timetree");

    if (!AddFiles(dir_, treeChain_, timeChain_))
        return kFALSE;

    if (!AddFiles(dir2_, treeChain2_, timeChain2_))
        return kFALSE;

    std::cout << "Reading chain1" << std::endl;

//    fileRoot_   = new TFile(inputName_, "READ");
//    std::cout << "Getting tree1" << std::endl;
//    tree_ = reinterpret_cast<TTree *>(fileRoot_->Get(tree));
    std::cout << " - setting branches" << std::endl;
    treeChain_->SetBranchAddress("Size"    , &Size_);
    treeChain_->SetBranchAddress("Col"     ,  Cols_);
    treeChain_->SetBranchAddress("Row"     ,  Rows_);
    treeChain_->SetBranchAddress("ToT"     ,  ToTs_);
    treeChain_->SetBranchAddress("ToA"     ,  ToAs_);
    treeChain_->SetBranchAddress("ToATrig" ,  ToATrigs_);    // l for long unsigned
    treeChain_->SetBranchAddress("TrigCntr", &TrigId_);
    treeChain_->SetBranchAddress("TrigTime", &TrigTime_);
    treeChain_->SetBranchAddress("TrigTimeNext", &TrigTimeNext_);

//    timeTree_ = reinterpret_cast<TTree *>(fileRoot_->Get("timetree"));
//    std::cout << " - setting branches" << std::endl;
    timeChain_->SetBranchAddress("TrigTime", &TrigWalk_);

    Entries_ = treeChain_->GetEntries();

    std::cout << "Reading chain2" << std::endl;
//    fileRoot2_   = new TFile(inputName2_, "READ");
//    std::cout << "Getting tree2" << std::endl;
//    tree2_ = reinterpret_cast<TTree *>(fileRoot2_->Get(tree));
    std::cout << " - setting branches" << std::endl;
    treeChain2_->SetBranchAddress("Size", &Size2_);
    treeChain2_->SetBranchAddress("Col",  Cols2_);
    treeChain2_->SetBranchAddress("Row",  Rows2_);
    treeChain2_->SetBranchAddress("ToT",  ToTs2_);
    treeChain2_->SetBranchAddress("ToA",  ToAs2_);
    treeChain2_->SetBranchAddress("ToATrig" ,  ToATrigs2_);    // l for long unsigned
    treeChain2_->SetBranchAddress("TrigCntr", &TrigId2_);
    treeChain2_->SetBranchAddress("TrigTime", &TrigTime2_);
    treeChain2_->SetBranchAddress("TrigTimeNext", &TrigTimeNext2_);

//    timeChain2_ = reinterpret_cast<TTree *>(fileRoot2_->Get("timetree"));
//    std::cout << " - setting branches" << std::endl;
    timeChain2_->SetBranchAddress("TrigTime", &TrigWalk2_);

    Entries2_ = treeChain2_->GetEntries();

    TrigIdTmp_  = 1;
    TrigIdTmp2_ = 1;
    std::cout << "Init trig " << TrigIdTmp_ << " " << TrigIdTmp2_ << std::endl;
    treeChain_->GetEntry(0);
    treeChain2_->GetEntry(0);
    std::cout << "Check trig " << TrigIdTmp_ << " " << TrigIdTmp2_ << std::endl;

    std::cout << "Time offset 1: " << time << "s, 2: " << time2 << "s" << std::endl;

    ToAzero_  = ToAs_[0]  + static_cast<ULong64_t>(163840000000 * static_cast<ULong64_t>(time));
    ToAzero2_ = ToAs2_[0] + static_cast<ULong64_t>(163840000000 * static_cast<ULong64_t>(time2));

    for (entry_ = 0; entry_ < Entries_; entry_++)
    {
        treeChain_->GetEntry(entry_);
        if (TrigTime_ > ToAzero_)
        {
            ToAzero_ = TrigTime_;
            std::cout << "Trig entry: " << TrigId_ << std::endl;
            break;
        }
    }

    TrigIdNext_ += TrigId_;
    TrigIdTmp_  += TrigId_;

    if (TrigStart_ >= TrigId_)
        TrigStart_ -= TrigId_;
    else
        TrigStart_ = 0;

    for (entry2_ = 0; entry2_ < Entries2_; entry2_++)
    {
        treeChain2_->GetEntry(entry2_);
        if (TrigTime2_ > ToAzero2_)
        {
            ToAzero2_ = TrigTime2_;
            std::cout << "Trig entry2: " << TrigId2_ << std::endl;
            break;
        }
    }

    TrigIdTmp2_ += TrigId2_;

    std::cout << "Shifting trig entry by " << TrigStart_ << std::endl;

    std::cout << "Zero time " << ToAzero_*25.0/4096000000000 << "s" << std::endl;
    std::cout << "Zero time2 " << ToAzero2_*25.0/4096000000000 << "s" << std::endl;

    std::cout << "Create writing file" << std::endl;
    outputName_ = inputName_;
    if (dir_->IsDirectory())
        outputName_.Append("/entangled" + name + "_" + std::to_string(time) + "_" + std::to_string(time2) + "_"+tree+"_processed.root");
    else
        outputName_.ReplaceAll(".root", name + "_" + std::to_string(time) + "_" + std::to_string(time2) + "_"+tree+"_processed.root");
    outputRoot_ = new TFile(outputName_, "RECREATE");

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

    return 1;
}

Bool_t Dual::AddFiles(TSystemFile* dir, TChain* chainDat, TChain* chainTime)
{
    TString name = dir->GetName();
    if (dir->IsDirectory())
    {
        TSystemDirectory folder(name,name);
        std::cout << "Directory: " << name << std::endl;
        TList *files = folder.GetListOfFiles();
        if (files)
        {
            files->Sort(kSortAscending);
            TSystemFile *file;
            TString fname;
            TIter next(files);
            while ((file=(TSystemFile*)next()))
            {
                fname = file->GetName();
                if (!file->IsDirectory() && fname.EndsWith(".root") && !fname.EndsWith("processed.root"))
                {
                    std::cout << fname << " at " << name << std::endl;
                    chainDat->Add(name + "/" + fname);
                    chainTime->Add(name + "/" + fname);
                }
            }
        }
        else
            return kFALSE;
    }
    else
    {
        std::cout << "File: " << name << std::endl;
        chainDat->Add(name);
        chainTime->Add(name);
    }
    return kTRUE;
}


void Dual::Process()
{
    outputRoot_->cd();

    Long64_t startEntry = entry_;
    Long64_t endEntry = Entries_;
    UInt_t   start = TrigId_;

    std::cout << "---------------------------------------------" << std::endl;
    std::cout << "Starting at: entry1 " << entry_ << ", entry2 " << entry2_ << std::endl;
    for (; entry_ < endEntry; entry_++)
    {
        treeChain_->GetEntry(entry_);
        treeChain2_->GetEntry(entry2_);

        if (maxEntries_ != 0 && entry_ >= startEntry + maxEntries_)
        {
            std::cout << "reached max entries" << std::endl;
            break;
        }

        if ((TrigId_ - start) < TrigStart_)
            continue;
        if ((TrigId_ - start) > 250 + TrigStart_)
        {
            std::cout << "reached 250s from " << TrigStart_ << std::endl;
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
        treeChain2_->GetEntry(pair);

        // area of incoming photons
        if ( ( TrigId2_ == (TrigIdTmp2_ - 1)) && PositionCheck2(area))
        {
            if ( static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) + ToAdiff_ >= static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_))
            && ( static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) + ToAdiff_ - static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_))) < diffToA )
            {
                diffToA = static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) + ToAdiff_ - static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_));
                nextEntry = pair;
                bFound = kTRUE;
            }
            else if ( static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)) + ToAdiff_ < static_cast<ULong64_t>(ToATrigs2_[0] *FindDelta(TrigDiff2_ ))
                 && ( static_cast<ULong64_t>(ToATrigs2_[0] *FindDelta(TrigDiff2_ )) - ToAdiff_ - static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_))) < diffToA )
            {
                diffToA = static_cast<ULong64_t>(ToATrigs2_[0] *FindDelta(TrigDiff2_ )) - ToAdiff_ - static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_));
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
    treeChain2_->GetEntry(entry2);

    TrigDiff2_ = TrigTimeNext2_ - TrigTime2_;
    UInt_t trig = TrigId2_;

    while((static_cast<ULong64_t>(ToATrigs2_[0]*FindDelta(TrigDiff2_)) < static_cast<ULong64_t>(ToATrigs_[0]*FindDelta(TrigDiff_)))
       && (TrigIdTmp2_ > TrigId2_)
       && (entry2 < Entries2_) )
    {
        treeChain2_->GetEntry(entry2++);
        if (TrigId2_ != trig)
        {
            std::cout << "Skipped trigger when searching starting entry2 " << TrigId2_ << std::endl;
            trig = TrigId2_;
        }
    }

    if (TrigIdTmp2_ == TrigId2_)
        return false;
    else
        return true;
}

void Dual::FindTrig(Long64_t &entry2)
{
    treeChain2_->GetEntry(entry2);

    UInt_t trig = TrigId2_;
    std::cout << "FT change " << TrigIdTmp2_ << " " << trig << " " << TrigId2_ << std::endl;

    while(TrigIdTmp2_ > TrigId2_ && entry2 < Entries2_)
    {
        treeChain2_->GetEntry(entry2++);

        if (TrigId2_ != trig)
        {
            std::cout << "FT Trig2 changed " << TrigIdTmp2_ << " " << trig << " " << TrigId2_ << std::endl;
            trig = TrigId2_;
        }
    }

    TrigIdTmp2_++;
}

void Dual::ScanEntry(Long64_t &entry, Long64_t &entry2)
{
    TrigDiff_ = TrigTimeNext_ - TrigTime_;

    if (TrigIdTmp_ <= TrigId_)
    {
        std::cout << "Entry1 " << entry << " of " << Entries_ << ", with trigger " << TrigId_ << " done!" << std::endl;
        std::cout << "Entry2 " << entry2 << " of " << Entries2_ << ", with trigger " << TrigId2_ << " done!" << std::endl;

        std::cout << "------------------------------" << std::endl;
        std::cout << "next/curr " << TrigIdTmp_ << " " << TrigId_ << std::endl;
        std::cout << "next2/curr2 " << TrigIdTmp2_ << " " << TrigId2_ << std::endl;
        TrigIdTmp_++;
        FindTrig(entry2);
        std::cout << "------------------------------" << std::endl;
    }

    if (entry % 1000 == 0)
    {
        std::cout << "------------------------------" << std::endl;
        std::cout << "Entry1 " << entry << " of " << Entries_ << ", with trigger " << TrigId_ << " done!" << std::endl;
        std::cout << "Entry2 " << entry2 << " of " << Entries2_ << ", with trigger " << TrigId2_ << " done!" << std::endl;
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
                treeChain2_->GetEntry(nextEntry_);
                entTree_->Fill();
            }
        }
    }
}
