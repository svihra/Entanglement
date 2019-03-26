#include <iostream>
#include <TApplication.h>
#include <TStopwatch.h>

#include "analysis.h"
#include "dual.h"
#include "video.h"
#include "entangled.h"

using namespace std;

int main(int argc, char **argv)
{
    TStopwatch time;
    time.Start();

    int c;
    bool bash = false;
//    TString dir = ".";
    TString name = "";
    TString file = "";
    TString file2 = "";
    Int_t zero = 0;
    Int_t zero2 = 0;
    Float_t diff = 0;
    TString tree = "rawtree";
    UInt_t maxEntries = 0;
    Int_t startPart = -1;
    Int_t parts = 25;

    while( ( c = getopt (argc, argv, ":bn:d:f:g:t:e:s:p:z:x:") ) != -1 )
    {
        switch(c)
        {
            case 'b':
                bash = true;
                break;
            case 'n':
                if (optarg){ name = (TString) optarg; break;}
                else {std::cout << "Missing name" << std::endl;}
            case 'd':
                if (optarg){ diff = (Float_t) std::atof(optarg); break;}
                else {std::cout << "Missing diff value" << std::endl;}
            case 'f':
                if(optarg){ file = (TString) optarg; break;}
                else {std::cout << "Missing file name" << std::endl;}
            case 'g':
                if(optarg){ file2 = (TString) optarg; break;}
                else {std::cout << "Missing file name" << std::endl;}
            case 't':
                if(optarg){ tree = (TString) optarg; break;}
                else {std::cout << "Missing tree definition" << std::endl;}
            case 'e':
                if(optarg){ maxEntries = (UInt_t) std::atoll(optarg); break;}
                else {std::cout << "Define max entries " << std::endl;}
            case 's':
                if(optarg){ startPart = (UInt_t) std::atoi(optarg); break;}
                else {std::cout << "Define start part " << std::endl;}
            case 'p':
                if(optarg){ parts = (UInt_t) std::atoi(optarg); break;}
                else {std::cout << "Define number of parts " << std::endl;}
            case 'z':
                if(optarg){ zero  = (Int_t) std::atoi(optarg); break;}
                else {std::cout << "Define zero " << std::endl;}
            case 'x':
                if(optarg){ zero2 = (Int_t) std::atoi(optarg); break;}
                else {std::cout << "Define zero2 " << std::endl;}
            case '?':
                std::cout << "Received unknown argument " << c << std::endl;
            default :
                std::cout << "Exiting code. Usage: " << std::endl;
                std::cout << "-b executes in bash" << std::endl;
                std::cout << "-d timediff in us" << std::endl;
                std::cout << "-f,g [] name of file to process " << std::endl;
                std::cout << "-t [] rawtree/proctree" << std::endl;
                std::cout << "-e [] max number of entries" << std::endl;
                std::cout << "-s [] starting part" << std::endl;
                std::cout << "-p [] number of parts" << std::endl;
                std::cout << "-z,x [] zero in seconds" << std::endl;
                return -1;
            break;
        }
    }
    if (bash)
        new Entangled(file, tree, maxEntries, startPart, parts);
//        new Dual(file, file2, diff, zero, zero2, tree, maxEntries, name);
    else
    {
        TApplication theApp("App",&argc,argv);

        new Analysis();
//        new Video();
        theApp.Run();
    }
    time.Stop();
    time.Print();
    return 0;
}

