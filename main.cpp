#include <iostream>
#include "analysis.h"
#include "entangled.h"
#include <TApplication.h>
#include <TStopwatch.h>

using namespace std;

int main(int argc, char **argv)
{
    TStopwatch time;
    time.Start();

    int c;
    bool bash = false;
//    TString dir = ".";
    TString file = "";
    TString tree = "rawtree";
    UInt_t maxEntries = 0;

    while( ( c = getopt (argc, argv, ":bf:t:e:") ) != -1 )
    {
        switch(c)
        {
            case 'b':
                bash = true;
                break;
//            case 'd':
//                if (optarg){ dir = (TString) optarg; combine = true; break;}
//                else {std::cout << "Missing dir name" << std::endl;}
            case 'f':
                if(optarg){ file = (TString) optarg; break;}
                else {std::cout << "Missing file name" << std::endl;}
            case 't':
                if(optarg){ tree = (TString) optarg; break;}
                else {std::cout << "Missing tree definition" << std::endl;}
            case 'e':
                if(optarg){ maxEntries = (UInt_t) std::atoll(optarg); break;}
                else {std::cout << "Define max entries " << std::endl;}
            case '?':
                std::cout << "Received unknown argument " << c << std::endl;
            default :
                std::cout << "Exiting code. Usage: " << std::endl;
                std::cout << "-b executes in bash" << std::endl;
//                std::cout << "-d name of folder in which to combine" << std::endl;
                std::cout << "-f [] name of file to process " << std::endl;
                std::cout << "-t [] rawtree/proctree" << std::endl;
                std::cout << "-e [] max number of entries" << std::endl;
                return -1;
            break;
        }
    }
    if (bash)
        new Entangled(file, tree, maxEntries);
    else
    {
        TApplication theApp("App",&argc,argv);

        new Analysis();
        theApp.Run();
    }
    time.Stop();
    time.Print();
    return 0;
}

