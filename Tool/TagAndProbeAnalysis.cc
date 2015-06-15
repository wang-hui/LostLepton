#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>

#include <vector>

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/NTupleReader.h"

#include "TStopwatch.h"
#include "TString.h"

#include "TagAndProbeAnalysis.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "TROOT.h"
//#include "TInterpreter.h"

using namespace std;

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr <<"Please give 2 arguments " << "runList " << " " << "outputFileName" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./TagAndProbeAnalysis runlist_TagAndProbe.txt test.root" << std::endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];

  //TChain *fChain = new TChain("stopTreeMaker/AUX");
  TChain *fChain = new TChain("muonEffs/fitter_tree");
  //define my TagAndProbeEffs class to stroe counts and efficiencies
  TagAndProbeEffs myTagAndProbeEffs;
  //define my histgram class
  //BaseHistgram myBaseHistgram;
  //myBaseHistgram.BookHistgram(outFileName);

  if(!FillChain(fChain, inputFileList))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }

  //clock to monitor the run time
  size_t t0 = clock();

  Float_t eta;
  fChain->SetBranchAddress("eta", &eta);
  Float_t pt;
  fChain->SetBranchAddress("pt", &pt);
  Float_t activity;
  fChain->SetBranchAddress("activity", &activity);
  Float_t miniiso;
  fChain->SetBranchAddress("miniIso", &miniiso);

  //NTupleReader tr(fChain);
  //tr.registerFunction(&Variables);

  std::cout<<"The loop begin: "<<std::endl;
  int nentries = (int)fChain->GetEntries();
  for (int i=0; i < nentries; i++) 
  {
    fChain->GetEvent(i);
    //std::cout << pt <<" yes" << std::endl;
    int ptbin_number = Set_ptbin_number(pt);
    int acbin_number = Set_acbin_number(activity);

    myTagAndProbeEffs.nmus_reco[ptbin_number][acbin_number]++;
    
    if (miniiso < 0.2)
    {
      myTagAndProbeEffs.nmus_iso[ptbin_number][acbin_number]++;
    }
  }
  /*
  while(tr.getNextEvent())
  {
    if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

    //float eta = tr.getVar<float>("eta");
    //int test = tr.getVar<int>("passingGlb");
    std::cout << eta <<"yes" << std::endl;

    
    double eta = tr.getVar<float>("eta");
    double pt = tr.getVar<float>("pt");
    double activity = tr.getVar<float>("activity");
    double miniiso = tr.getVar<float>("miniIso");
    double pairmass = tr.getVar<float>("mass");

    int ptbin_number = Set_ptbin_number(pt);
    int acbin_number = Set_acbin_number(activity);

    myTagAndProbeEffs.nmus_reco[ptbin_number][acbin_number]++;
    
    if (miniiso < 0.2)
    {
      myTagAndProbeEffs.nmus_iso[ptbin_number][acbin_number]++;
    }
  }
  */
  myTagAndProbeEffs.NumberstoEffs();
  myTagAndProbeEffs.TagAndProbe2dPlotsGen();
}

void TagAndProbeEffs::NumberstoEffs()
{
  int i_cal;
  int j_cal;
  
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      nmus_reco[i_cal][j_cal] = nmus_reco[i_cal][j_cal]*scale;
      nmus_iso[i_cal][j_cal] = nmus_iso[i_cal][j_cal]*scale;

      mus_isoeff[i_cal][j_cal]     = nmus_iso[i_cal][j_cal]/nmus_reco[i_cal][j_cal];
      mus_isoeff_err[i_cal][j_cal] = get_stat_Error(nmus_iso[i_cal][j_cal],nmus_reco[i_cal][j_cal]);
    }
  }
}

void TagAndProbeEffs::TagAndProbe2dPlotsGen()
{
  int i_cal;
  int j_cal;

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      mus_nreco2d->SetBinContent( i_cal+1 , j_cal+1, nmus_reco[i_cal][j_cal] );
      mus_niso2d->SetBinContent( i_cal+1 , j_cal+1, nmus_iso[i_cal][j_cal] );
      mus_isoeffs2d->SetBinContent( i_cal+1 , j_cal+1, mus_isoeff[i_cal][j_cal] );

      mus_isoeffs2d->SetBinError( i_cal+1 , j_cal+1, mus_isoeff_err[i_cal][j_cal] );
    }
  }


  mus_nreco2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_nreco2d->GetYaxis()->SetTitle("Activity");
  mus_niso2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_niso2d->GetYaxis()->SetTitle("Activity");
  mus_isoeffs2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_isoeffs2d->GetYaxis()->SetTitle("Activity");

  TagAndProbe2dPlots->Write();
}



