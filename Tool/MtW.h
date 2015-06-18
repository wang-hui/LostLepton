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

#include "LLBinFunction.h"

#define N_FILES 5

class BaseHistgram
{
 public:
  void BookHistgram(const char *);
  TFile *oFile;
  TH1D *h_b_mtw[N_FILES];
};

void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");

  h_b_mtw[0] = new TH1D("h_b_mtw_ttbar","",20,0,500);
  h_b_mtw[1] = new TH1D("h_b_mtw_t2tt_1","",20,0,500);
  h_b_mtw[2] = new TH1D("h_b_mtw_t2tt_2","",20,0,500);
  h_b_mtw[3] = new TH1D("h_b_mtw_t2tt_3","",20,0,500);
  h_b_mtw[4] = new TH1D("h_b_mtw_t2tt_4","",20,0,500);
}

//Fill chain from txt file
bool FillChain(TChain *chain, const TString &inputFileList)
{
  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open())
  {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return false;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1)
  {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
  }
  std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
  return true;
}

//##########functions to calculate Delta_R and Delta Phi###############
double DeltaPhi(double phi1, double phi2) 
{
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) 
{
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}
