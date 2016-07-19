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
#include "SusyAnaTools/Tools/searchBins.h"

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
#include "TList.h"
#include "TPad.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TString.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "TROOT.h"
//#include "TInterpreter.h"

#include "MtW.h"

int main(int argc, char* argv[])
{
  if (argc < N_FILES+1)
  {
    std::cerr <<"Please give 6 arguments " << "runList1 " << " " << "runList2 " << " " << "runList3 " << " " << "runList4 " << " "<< "runList5 " << " " << "outputFileName" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./MtW runlist_ttbar.txt runlist_t2tt_1.txt runlist_t2tt_2.txt runlist_t2tt_3.txt runlist_t2tt_4.txt 123.root" << std::endl;
    return -1;
  }
  const char *outFileName   = argv[N_FILES+1];
  //define my histgram class
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  double Lumi = 10000;
  double weight[N_FILES];
  //ttbar weight
  weight[0] = Lumi*806.1/25446933;
  //t2tt 425 325
  weight[1] = Lumi*1.31169/1039030;
  //t2tt 500 325
  weight[2] = Lumi*0.51848/109591;
  //t2tt 650 325
  weight[3] = Lumi*0.107045/105672;
  //t2tt 850 100
  weight[4] = Lumi*0.0189612/102839;

  TChain *fChain[N_FILES];
  fChain[0] = new TChain("stopTreeMaker/AUX");
  fChain[1] = new TChain("stopTreeMaker/AUX");
  fChain[2] = new TChain("stopTreeMaker/AUX");
  fChain[3] = new TChain("stopTreeMaker/AUX");
  fChain[4] = new TChain("stopTreeMaker/AUX");


  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(spec);

  for(int i = 0 ; i < N_FILES ; i++ )
  {
    const char *inputFileList = argv[i+1];
 
    if(!FillChain(fChain[i], inputFileList))
    {
      std::cerr << "Cannot get the tree "<< i+1 << std::endl;
    }

    //clock to monitor the run time
    size_t t0 = clock();

    NTupleReader tr(fChain[i]);
    //initialize the type3Ptr defined in the customize.h
    AnaFunctions::prepareTopTagger();
    //The passBaselineFunc is registered here
    tr.registerFunction(&passBaselineFunc);
 
    std::cout << i+1 << " loop begin: " << std::endl;
    while(tr.getNextEvent())
    {
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      bool passBaselinelostlept = tr.getVar<bool>("passBaseline"+spec);
 
      if(
         passBaselinelostlept
        )
      {
        int nElectrons = tr.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = tr.getVar<int>("nMuons_CUT"+spec);

        double met = tr.getVar<double>("met");
        double metphi = tr.getVar<double>("metphi");

        int njets30 = tr.getVar<int>("cntNJetsPt30Eta24"+spec);
        int ntopjets = tr.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = tr.getVar<int>("cntCSVS"+spec);
        double MT2 = tr.getVar<double>("best_had_brJet_MT2"+spec);

        //get muon variables
        std::vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
        std::vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");

        double reco_mus_pt = -1, reco_mus_eta = 0, reco_mus_phi = 0;

        if (nElectrons == 0 && nMuons == 1)
        {

          for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
          {
            if( fabs(muonsLVec[im].Eta()) < (AnaConsts::muonsMiniIsoArr).maxAbsEta && muonsMiniIso[im] < (AnaConsts::muonsMiniIsoArr).maxIso )
            {
              reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
              reco_mus_eta = ( muonsLVec.at(im) ).Eta();
              reco_mus_phi = ( muonsLVec.at(im) ).Phi();
              double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
              double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * ( 1.0 - cos(deltaphi_mus) ) );
              (myBaseHistgram.h_b_mtw[i])->Fill(mtW_mus,weight[i]);
            }
          }
        }
      }
    }
  }

  (myBaseHistgram.oFile)->Write();
}
