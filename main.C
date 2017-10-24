// HEADERS
// C++
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdlib.h>

// ROOT
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

// USER DEFINITIONS
#include "TMDB.h"



int main(int argc, char **argv) {

  // split by ',' files from physics runs (L1TGCNtuple)
  std::string argStrPhysics = argv[1];
  std::vector<std::string> fileListPhysics;
  for (size_t i=0,n; i <= argStrPhysics.length(); i=n+1)
    {
      n = argStrPhysics.find_first_of(',',i);
      if (n == std::string::npos)
        n = argStrPhysics.length();
      std::string tmp = argStrPhysics.substr(i,n-i);
      fileListPhysics.push_back(tmp);
    }

  // open physics input files
  TChain* fChainPhysics = new TChain("physics");
  std::cout << "\nPhysics files:" << endl;
  for (unsigned int iFile=0; iFile<fileListPhysics.size(); ++iFile)
    {
      std::cout << "Open: " << fileListPhysics[iFile].c_str() << std::endl;
      fChainPhysics->Add(fileListPhysics[iFile].c_str());
    }


  gROOT->SetBatch(true);

  TApplication theApp("App", &argc, argv);
  gApplication->Init();

  argc=theApp.Argc();
  argv=theApp.Argv();

  std::string fout = argv[2];

// ########## Let's play! ###########

	TMDB * coin = new TMDB(fChainPhysics);
	coin->UserInit();							// Initializing attributes/variables

	// #### Physics functions
	coin->LoopPhysics	( fout );			// Read entries from physics runs (L1TGCNtuple)

	// #### Print functions
	// coin->print_latency(fout);			// print TMDB Latency

  return EXIT_SUCCESS;
}
