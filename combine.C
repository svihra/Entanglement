#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TDirectory.h>

#include <fstream>

void save(TString fileName, TFile* fileRoot, Int_t x1, Int_t y1,Int_t x2, Int_t y2)
{
    TString csvName = fileName;
    csvName.Remove(csvName.Last('.'),200);
    csvName.ReplaceAll("combined_","");

    if (x1 != -1)
        csvName.Append(Form("_%dx%d_%dx%d.csv",x1,y1,x2,y2));
    else
        csvName.Append(".csv");

    TString entName, single1Name, single2Name;
    if (x1 != -1)
    {
        entName.Form("histEnt_%dx%d_%dx%d",x1,y1,x2,y2);
        single1Name.Form("histSingle1_%dx%d_%dx%d",x1,y1,x2,y2);
        single2Name.Form("histSingle2_%dx%d_%dx%d",x1,y1,x2,y2);
    }
    else
    {
        entName = "histEnt";
        single1Name = "histSingle1";
        single2Name = "histSingle2";
    }
    TH1F* ent = (TH1F*) fileRoot->Get(entName);
    TH1F* single1 = (TH1F*) fileRoot->Get(single1Name);
    TH1F* single2 = (TH1F*) fileRoot->Get(single2Name);

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

void combine(TString fileName)
{
    if (fileName.EndsWith("_combined_processed.root"))
    {
        std::cout << "Reading tree" << std::endl;
        TFile* fileRoot = new TFile(fileName, "READ");

        for (Int_t x1 = 0; x1 < 3; x1++)
        {
            for (Int_t y1 = 0; y1 < 3; y1++)
            {
                for (Int_t x2 = 0; x2 < 3; x2++)
                {
                    for (Int_t y2 = 0; y2 < 3; y2++)
                    {
                        save(fileName,fileRoot,x1,y1,x2,y2);
                    }
                }
            }
        }
        std::cout << "Saving single" << std::endl;
        save(fileName,fileRoot,-1,-1,-1,-1);
        std::cout << "Finished " << fileName << std::endl;
    }
}
