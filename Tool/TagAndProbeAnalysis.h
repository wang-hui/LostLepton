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

#define PT_BINS 8
#define AC_BINS 8
#define NJETS_BINS 6

//############################begin to defin class AccRecoIsoEffs###################


class TagAndProbeEffs
{
 public:
  double mus_recoeff[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff[PT_BINS][AC_BINS] = {{0}};

  double mus_recoeff_err[PT_BINS][AC_BINS] = {{0}};
  double mus_isoeff_err[PT_BINS][AC_BINS] = {{0}};

  //here we define the muon/electron number we need to count in the loop
  double nmus_acc_bin[PT_BINS][AC_BINS] = {{0}};
  double nmus_reco[PT_BINS][AC_BINS] = {{0}};
  double nmus_iso[PT_BINS][AC_BINS] = {{0}};

  TFile *TagAndProbe2dPlots = new TFile("TagAndProbe2dPlots.root", "recreate");
  //double ptbins[9]={5.0,10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
  //double acbins[9]={0.0,5.0,10.0,20.0,40.0,60.0,80.0,100.0,120.0};
  //TH2D *mus_recoeffs2d  = new TH2D("mus_recoeffs","Muon RecoEffs",8,ptbins,8,acbins);
  //TH2D *mus_isoeffs2d  = new TH2D("mus_isoeffs","Muon IsoEffs",8,ptbins,8,acbins);
  //TH2D *els_recoeffs2d  = new TH2D("els_recoeffs","Electron RecoEffs",8,ptbins,8,acbins);
  //TH2D *els_isoeffs2d  = new TH2D("els_isoeffs","Electron IsoEffs",8,ptbins,8,acbins);

  TH2D *mus_nreco2d = new TH2D("mus_nreco_tagandprobe","TagAndProbe Muon NReco",8,0,8,8,0,8);
  TH2D *mus_niso2d = new TH2D("mus_niso_tagandprobe","TagAndProbe Muon NIso",8,0,8,8,0,8);
  TH2D *mus_isoeffs2d  = new TH2D("mus_isoeffs_tagandprobe","TagAndProbe Muon IsoEffs",8,0,8,8,0,8);

  void NumberstoEffs();
  void TagAndProbe2dPlotsGen();

 private:
  //here we define the variables we needed for normalization
  double XSec = 4746;
  double Lumi = 1000.0;
  double Nevents = 2829164;
  double scale = XSec*Lumi/Nevents;

  double get_stat_Error(double a,
                        double an
                       );
  
  double get_sys_Error(double r,
                       double p
                      );
};


double TagAndProbeEffs::get_stat_Error(double a,
                                       double an
                                      )
{
  double n;
  n = an - a;

  double err;
  err = 1000;

  double alpha;
  alpha = 1-0.6827;

  if( a>=0 && n>=0 )
  {
    err = std::sqrt(n/(a+n)/(a+n)*n/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,a+1,1)+a/(a+n)/(a+n)*a/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,n+1,1));
    return err;
  }
  else
  {
    return -1;
  }
}

double TagAndProbeEffs::get_sys_Error(double r,
                                      double p
                                     )
{
  double err;
  err = 1000;

  err = r - (1-(1-r)*(1+p));

  if(err < 0)
  {
    return -1;
  }
  return err;
}

//############finish the definition of class AccRecoEffs######################

class BaseHistgram
{
 public:
  void BookHistgram(const char *);

  TFile *oFile;
};

void BaseHistgram::BookHistgram(const char *outFileName)
{
  oFile = new TFile(outFileName, "recreate");
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


//############determine the pt bin number############
int Set_ptbin_number(
                     double gen_pt
                    )
{
  int ptbin_num;

  if(gen_pt < 10)
  {
    ptbin_num = 0;
  }
  else if(gen_pt >= 10 && gen_pt < 20)
  {
    ptbin_num = 1;
  }
  else if(gen_pt >= 20 && gen_pt < 30)
  {
    ptbin_num = 2;
  }
  else if(gen_pt >= 30 && gen_pt < 40)
  {
    ptbin_num = 3;
  }
  else if(gen_pt >= 40 && gen_pt < 50)
  {
    ptbin_num = 4;
  }
  else if(gen_pt >= 50 && gen_pt < 70)
  {
    ptbin_num = 5;
  }
  else if(gen_pt >= 70 && gen_pt < 100)
  {
    ptbin_num = 6;
  }
  else if(gen_pt >= 100 )
  {
    ptbin_num = 7;
  }
  
  return ptbin_num;
}

//############determine the activity bin number############

int Set_acbin_number(
                     double activity
                    )
{
  int acbin_num;

  if(activity < 5)
  {
    acbin_num = 0;
  }
  else if(activity >= 5 && activity < 10)
  {
    acbin_num = 1;
  }
  else if(activity >= 10 && activity < 20)
  {
    acbin_num = 2;
  }
  else if(activity >= 20 && activity < 40)
  {
    acbin_num = 3;
  }
  else if(activity >= 40 && activity < 60)
  {
    acbin_num = 4;
  }
  else if(activity >= 60 && activity < 80)
  {
    acbin_num = 5;
  }
  else if(activity >= 80 && activity < 100)
  {
    acbin_num = 6;
  }
  else if(activity >= 100 )
  {
    acbin_num = 7;
  }

  return acbin_num;
}

int Set_njetsbin_number(
                        int njets
                       )
{
  int njetsbin_num;

  if(njets == 4)
  {
    njetsbin_num = 0;
  }
  else if(njets == 5)
  {
    njetsbin_num = 1;
  }
  else if(njets == 6)
  {
    njetsbin_num = 2;
  }
  else if(njets == 7)
  {
    njetsbin_num = 3;
  }
  else if(njets ==8)
  {
    njetsbin_num = 4;
  }
  else if(njets >= 9)
  {
    njetsbin_num = 5;
  }

  return njetsbin_num;
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
