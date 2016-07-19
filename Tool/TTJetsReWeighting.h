#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "ConstantsSnippet.h"

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
  void TTJetsSampleInfo_push_back( std::string tag, double xsec, double nevents, double lumi, double kf, const TString &inputFileList );
 private:
  bool FillChain(TChain *chain, const TString &inputFileList, std::string tag);
};

