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
#define LUMI 3000
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
  
  void FillTTJetsSampleInfos( const TString &inputFileList );
 private:
  void TTJetsSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, const TString &inputFileList );
};

void TTJetsSampleWeight::FillTTJetsSampleInfos( const TString &inputFileList )
{
  TTJetsSampleInfos.clear();
 
  double W_Lept_BR = 0.1086*3;
  double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
  double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2
  //TTJets nominal
  TTJetsSampleInfo_push_back( "TTJets_", 831.76, 11339232, LUMI, inputFileList );
  //TTJets single lepton and di-lepton
  //TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromT_", 831.76*0.5*TTbar_SingleLept_BR, 60144642, LUMI, inputFileList );
  //TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromTbar_", 831.76*0.5*TTbar_SingleLept_BR, 59816364, LUMI, inputFileList );
  //TTJetsSampleInfo_push_back( "TTJets_DiLept_", 831.76*TTbar_DiLept_BR, 30498962, LUMI, inputFileList );
}

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
