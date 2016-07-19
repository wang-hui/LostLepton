#include "TTJetsReWeighting.h"

bool TTJetsSampleWeight::FillChain(TChain *chain, const TString &inputFileList, std::string tag)
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
