#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
//define lumi in pb-1
//#define LUMI 567.945
//#define LUMI 2100.0
//#define LUMI 3000.0
//#define LUMI 2154.5
//#define LUMI 2262.946
//#define LUMI 8000.0
//#define LUMI 7647.637518921
//#define LUMI 12905.3
#define LUMI 12877.0846508279992

//Fill chain from txt file
bool FillChain(TChain *chain, const TString &inputFileList, std::string tag)
{
  std::ifstream infile( inputFileList, std::ifstream::in );
  std::string buffer;

  if(!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1)
  {
    buffer.clear();
    infile >> buffer;

    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;
    //std::cout << (buffer.find(tag) != std::string::npos) << std::endl;
    if (buffer.find(tag) != std::string::npos) 
    {
      //std::cout << "found!" << '\n';
      chain->Add(buffer.c_str());
    }
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return true;
}

struct TTJetsSampleInfo
{
  std::string TTJetsTag;
  double weight;
  TChain *chain;
};

class TTJetsSampleWeight
{
 public:
  std::vector<TTJetsSampleInfo> TTJetsSampleInfos;
  void TTJetsSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, const TString &inputFileList );
};

void TTJetsSampleWeight::TTJetsSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, const TString &inputFileList)
{
  TTJetsSampleInfo oneInfo;

  oneInfo.TTJetsTag = tag;
  oneInfo.weight = xsec*lumi/nevents;
  //oneInfo.chain= new TChain("AUX");
  oneInfo.chain= new TChain("stopTreeMaker/AUX");
  if(!FillChain(oneInfo.chain, inputFileList, oneInfo.TTJetsTag))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }
  TTJetsSampleInfos.push_back(oneInfo);
  oneInfo = {};
}
