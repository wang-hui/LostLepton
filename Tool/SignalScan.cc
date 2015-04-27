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

#include "TStopwatch.h"
#include "TString.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "LostLepton_MuCS_TTbar.h"
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
//#include "Activity.h"

using namespace std;

//const bool debug = true;

void passBaselineFunc(NTupleReader &tr)
{
  bool passBaseline_RelIso = true;
  bool passBaseline_MiniIso = true;
  //bool passBaseline_nolepveto = true;

  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);

  //Calculate number of leptons
  int nMuonsRelIso = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
  int nMuonsMiniIso = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);

  int nElectronsRelIso = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
  int nElectronsMiniIso = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);

  int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);

  //Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
  int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
  int cntNJetsPt30      = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);

  //Calculate deltaPhi
  std::vector<double> * dPhiVec = new std::vector<double>();
  (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);

  //Prepare jets and b-tag working points for top tagger
  //The jets stored in flat ntuples might have looser pt or eta requirement (here it's pt>10 GeV in flat ntuple),
  //while for the top tagger we need higher pt requirement as defined in AnaConsts::pt30Arr array.
  std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
  AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
  //if( debug ) std::cout<<"\njetsLVec_forTagger->size : "<<jetsLVec_forTagger->size()<<"  recoJetsBtag_forTagger->size : "<<recoJetsBtag_forTagger->size()<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass lepton veto?
  //bool passLeptVeto = true;
  if( nMuonsRelIso != AnaConsts::nMuonsSel ){ passBaseline_RelIso = false; }
  if( nElectronsRelIso != AnaConsts::nElectronsSel ){ passBaseline_RelIso = false; }

  if( nMuonsMiniIso != AnaConsts::nMuonsSel ){ passBaseline_MiniIso = false; }
  if( nElectronsMiniIso != AnaConsts::nElectronsSel ){ passBaseline_MiniIso = false; }
  //Isolated track veto is disabled for now
  //if( nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passLeptVeto = false; }
  //if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; }
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; }
  //if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passdPhis = false; }
  //if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passBJets = false; }
  //if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passMET = false; }
  //if( debug ) std::cout<<"met : "<<tr.getVar<double>("met")<<"  defaultMETcut : "<<AnaConsts::defaultMETcut<<"  passBaseline : "<<passBaseline<<std::endl;

  //std::cout << "AnaConsts::defaultMETcut = " << AnaConsts::defaultMETcut << std::endl;

  //Calculate top tagger related variables. 
  //Note that to save speed, only do the calculation after previous base line requirements.
  int bestTopJetIdx = -1;
  bool remainPassCSVS = false;
  int pickedRemainingCombfatJetIdx = -1;
  double bestTopJetMass = -1;
  int nTopCandSortedCnt = 0;
  double MT2 = -1;
  double mTcomb = -1;

  //if( passBaseline && cntNJetsPt30 >= AnaConsts::nJetsSel ){
  if( (passBaseline_RelIso || passBaseline_MiniIso) && cntNJetsPt30 >= AnaConsts::nJetsSel )
  {
    type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
    bestTopJetIdx = type3Ptr->bestTopJetIdx;
    remainPassCSVS = type3Ptr->remainPassCSVS;
    pickedRemainingCombfatJetIdx = type3Ptr->pickedRemainingCombfatJetIdx;
    if( bestTopJetIdx != -1 ) bestTopJetMass = type3Ptr->bestTopJetLVec.M();

    nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
    MT2 = type3Ptr->MT2;
    mTcomb = type3Ptr->mTbJet + 0.5*type3Ptr->mTbestTopJet;
  }

  //Pass top tagger requirement?
  bool passTagger = true;
  //bestTopJetIdx != -1 means at least 1 top candidate!
  if( bestTopJetIdx == -1 ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passTagger = false; }
  if( ! remainPassCSVS ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passTagger = false; }
  if( pickedRemainingCombfatJetIdx == -1 && jetsLVec_forTagger->size()>=6 ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passTagger = false; }
  if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ ) ){ passBaseline_RelIso = false; passBaseline_MiniIso = false; passTagger = false; }
  //if( debug ) std::cout<<"bestTopJetidx : "<<bestTopJetIdx<<"  remainPassCSVS : "<<remainPassCSVS<<"  pickedRemainingCombfatJetIdx : "<<pickedRemainingCombfatJetIdx<<"  bestTopJetMass : "<<bestTopJetMass<<"  passBaseline : "<<passBaseline<<std::endl;

  //Register all the calculated variables
  //tr.registerDerivedVar("nMuons_CUT2", nMuons);
  //tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
  //tr.registerDerivedVar("nIsoTrks_CUT", nIsoTrks);
  //tr.registerDerivedVar("cntNJetsPt50Eta24", cntNJetsPt50Eta24);
  tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);
  //tr.registerDerivedVec("dPhiVec", dPhiVec);
  //tr.registerDerivedVar("cntCSVS", cntCSVS);
  //tr.registerDerivedVec("jetsLVec_forTagger", jetsLVec_forTagger);
  //tr.registerDerivedVec("recoJetsBtag_forTagger", recoJetsBtag_forTagger);
  //tr.registerDerivedVar("cntNJetsPt30", cntNJetsPt30);
  //tr.registerDerivedVar("bestTopJetIdx", bestTopJetIdx);
  //tr.registerDerivedVar("remainPassCSVS", remainPassCSVS);
  //tr.registerDerivedVar("pickedRemainingCombfatJetIdx", pickedRemainingCombfatJetIdx);
  tr.registerDerivedVar("bestTopJetMass2", bestTopJetMass);

  //All the pass booleans are stored/registered into the NTupleReader and can be used later
  //tr.registerDerivedVar("passLeptVeto", passLeptVeto);
  //tr.registerDerivedVar("passnJets", passnJets);
  //tr.registerDerivedVar("passdPhis", passdPhis);
  //tr.registerDerivedVar("passBJets", passBJets);
  //tr.registerDerivedVar("passMET", passMET);
  //tr.registerDerivedVar("passTagger", passTagger);
  tr.registerDerivedVar("passBaseline_RelIso", passBaseline_RelIso);
  tr.registerDerivedVar("passBaseline_MiniIso", passBaseline_MiniIso);

  //if( debug ) std::cout<<"nTopCandSortedCnt : "<<nTopCandSortedCnt<<"  MT2 : "<<MT2<<"  mTcomb : "<<mTcomb<<"  passBaseline : "<<passBaseline<<std::endl;

  tr.registerDerivedVar("nTopCandSortedCnt", nTopCandSortedCnt);
  tr.registerDerivedVar("MT22", MT2);
  //tr.registerDerivedVar("mTcomb", mTcomb);

  //if( debug ) std::cout<<"passBaseline : "<<passBaseline<<"  passBaseline : "<<passBaseline<<std::endl;
}


int main(int argc, char* argv[])
{

  if (argc < 2)
  {
    std::cerr <<"Please give 2 arguments " << "runList " << " " << "outputFileName" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./LostLepton_MuCS_TTbar runlist_ttjets.txt isoplots.root" << std::endl;
    std::cerr <<" ./LostLepton runlist_ttjets.txt isoplots.root" << std::endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];

  TChain *fChain = new TChain("stopTreeMaker/AUX");
  //TChain *fChain = new TChain("AUX");

  if(!FillChain(fChain, inputFileList))
  {
    std::cerr << "Cannot get the tree " << std::endl;
  }

  //clock to monitor the run time
  size_t t0 = clock();

  NTupleReader tr(fChain);
  //initialize the type3Ptr defined in the customize.h
  AnaFunctions::prepareTopTagger();
  //The passBaselineFunc is registered here
  tr.registerFunction(&passBaselineFunc);

  int nevents_total = 0;
  int nevents_baseline_RelIso = 0;
  int nevents_baseline_MiniIso = 0;

  while(tr.getNextEvent())
  {
    if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    nevents_total++;  

    //bool passBaseline=tr.getVar<bool>("passBaseline");
    //if (passBaseline)
    //{
      //++nevents_baseline_ref;
    //}
    bool passBaseline_RelIso = tr.getVar<bool>("passBaseline_RelIso");
    bool passBaseline_MiniIso = tr.getVar<bool>("passBaseline_MiniIso");
    if (passBaseline_RelIso)
    {
      nevents_baseline_RelIso++;
    }
    if (passBaseline_MiniIso)
    {
      nevents_baseline_MiniIso++;
    }
  }

  std::cout << "Total:"<< nevents_total << std::endl;
  std::cout << "Baseline, RelIso:"<< nevents_baseline_RelIso << std::endl;
  std::cout << "Baseline, MiniIso:"<< nevents_baseline_MiniIso << std::endl;

  return 0;
}
