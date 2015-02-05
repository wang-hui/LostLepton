#include <iostream>
#include <vector>
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

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/customize.h"

#include "TStopwatch.h"
#include "TString.h"

#include "NTupleReader.h"
#include "LostLepton.h"
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

using namespace std;

const bool debug = true;

void passBaselineFunc(NTupleReader &tr){
   bool passBaseline = true;
   bool passBaseline_nolepveto = true;

// Form TLorentzVector of MET
   TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);

// Calculate number of leptons
   int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
   int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
   int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), AnaConsts::isoTrksArr);

// Calculate number of jets and b-tagged jets
   int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), AnaConsts::cutCSVS, AnaConsts::bTagArr);
   int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt50Eta24Arr);
   int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Eta24Arr);
   int cntNJetsPt30      = AnaFunctions::countJets(tr.getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);

// Calculate deltaPhi
   std::vector<double> * dPhiVec = new std::vector<double>();
   (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVar<double>("metphi"), 3, AnaConsts::dphiArr);

// Prepare jets and b-tag working points for top tagger
// The jets stored in flat ntuples might have looser pt or eta requirement (here it's pt>10 GeV in flat ntuple),
// while for the top tagger we need higher pt requirement as defined in AnaConsts::pt30Arr array.
   std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
   AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>("jetsLVec"), tr.getVec<double>("recoJetsBtag_0"), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
   //if( debug ) std::cout<<"\njetsLVec_forTagger->size : "<<jetsLVec_forTagger->size()<<"  recoJetsBtag_forTagger->size : "<<recoJetsBtag_forTagger->size()<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass lepton veto?
   bool passLeptVeto = true;
   if( nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passLeptVeto = false; }
   if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passLeptVeto = false; }
// Isolated track veto is disabled for now
//   if( nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passLeptVeto = false; }
   //if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass number of jets?
   bool passnJets = true;
   if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passBaseline_nolepveto = false; passnJets = false;}
   if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passBaseline_nolepveto = false; passnJets = false;}
   //if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass deltaPhi?
   bool passdPhis = true;
   if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passBaseline_nolepveto = false; passdPhis = false; }
   //if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass number of b-tagged jets?
   bool passBJets = true;
   if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBaseline_nolepveto = false; passBJets = false; }
   //if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

// Pass the baseline MET requirement?
   bool passMET = true;
   if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passBaseline_nolepveto = false; passMET = false; }
   //if( debug ) std::cout<<"met : "<<tr.getVar<double>("met")<<"  defaultMETcut : "<<AnaConsts::defaultMETcut<<"  passBaseline : "<<passBaseline<<std::endl;

   //std::cout << "AnaConsts::defaultMETcut = " << AnaConsts::defaultMETcut << std::endl;

// Calculate top tagger related variables. 
// Note that to save speed, only do the calculation after previous base line requirements.
   int bestTopJetIdx = -1;
   bool remainPassCSVS = false;
   int pickedRemainingCombfatJetIdx = -1;
   double bestTopJetMass = -1;
   int nTopCandSortedCnt = 0;
   double MT2 = -1;
   double mTcomb = -1;

   //if( passBaseline && cntNJetsPt30 >= AnaConsts::nJetsSel ){
   if( passBaseline_nolepveto && cntNJetsPt30 >= AnaConsts::nJetsSel ){
      type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
      bestTopJetIdx = type3Ptr->bestTopJetIdx;
      remainPassCSVS = type3Ptr->remainPassCSVS;
      pickedRemainingCombfatJetIdx = type3Ptr->pickedRemainingCombfatJetIdx;
      if( bestTopJetIdx != -1 ) bestTopJetMass = type3Ptr->bestTopJetLVec.M();

      nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
      MT2 = type3Ptr->MT2;
      mTcomb = type3Ptr->mTbJet + 0.5*type3Ptr->mTbestTopJet;
   }

// Pass top tagger requirement?
   bool passTagger = true;
// bestTopJetIdx != -1 means at least 1 top candidate!
   if( bestTopJetIdx == -1 ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
   if( ! remainPassCSVS ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
   if( pickedRemainingCombfatJetIdx == -1 && jetsLVec_forTagger->size()>=6 ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
   if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ ) ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
   //if( debug ) std::cout<<"bestTopJetidx : "<<bestTopJetIdx<<"  remainPassCSVS : "<<remainPassCSVS<<"  pickedRemainingCombfatJetIdx : "<<pickedRemainingCombfatJetIdx<<"  bestTopJetMass : "<<bestTopJetMass<<"  passBaseline : "<<passBaseline<<std::endl;

//// Register all the calculated variables
   tr.registerDerivedVar("nMuons_CUT2", nMuons);
   tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
//   tr.registerDerivedVar("nIsoTrks_CUT", nIsoTrks);
//
//   tr.registerDerivedVar("cntNJetsPt50Eta24", cntNJetsPt50Eta24);
   tr.registerDerivedVar("cntNJetsPt30Eta24", cntNJetsPt30Eta24);
//
//   tr.registerDerivedVec("dPhiVec", dPhiVec);
//
//   tr.registerDerivedVar("cntCSVS", cntCSVS);
//
//   tr.registerDerivedVec("jetsLVec_forTagger", jetsLVec_forTagger);
//   tr.registerDerivedVec("recoJetsBtag_forTagger", recoJetsBtag_forTagger);
//
//   tr.registerDerivedVar("cntNJetsPt30", cntNJetsPt30);
//
//   tr.registerDerivedVar("bestTopJetIdx", bestTopJetIdx);
//   tr.registerDerivedVar("remainPassCSVS", remainPassCSVS);
//   tr.registerDerivedVar("pickedRemainingCombfatJetIdx", pickedRemainingCombfatJetIdx);
   tr.registerDerivedVar("bestTopJetMass2", bestTopJetMass);
//
//// All the pass booleans are stored/registered into the NTupleReader and can be used later
//   tr.registerDerivedVar("passLeptVeto", passLeptVeto);
//   tr.registerDerivedVar("passnJets", passnJets);
//   tr.registerDerivedVar("passdPhis", passdPhis);
//   tr.registerDerivedVar("passBJets", passBJets);
//   tr.registerDerivedVar("passMET", passMET);
//   tr.registerDerivedVar("passTagger", passTagger);
   tr.registerDerivedVar("passBaseline", passBaseline);
   tr.registerDerivedVar("passBaseline_nolepveto", passBaseline_nolepveto);
//
//   //if( debug ) std::cout<<"nTopCandSortedCnt : "<<nTopCandSortedCnt<<"  MT2 : "<<MT2<<"  mTcomb : "<<mTcomb<<"  passBaseline : "<<passBaseline<<std::endl;
//
//   tr.registerDerivedVar("nTopCandSortedCnt", nTopCandSortedCnt);
   tr.registerDerivedVar("MT22", MT2);
//   tr.registerDerivedVar("mTcomb", mTcomb);
//
   //if( debug ) std::cout<<"passBaseline : "<<passBaseline<<"  passBaseline : "<<passBaseline<<std::endl;
}


int main()
{

  //char nBase[] = "/afs/cern.ch/user/h/hua/stop/AllHadronicSUSY/CMSSW_7_2_0/stop_ttbar_skimmed_tree.root";
  char nBase[] = "stop_ttbar_skimmed_tree.root";

  TChain *fChain = new TChain("AUX");

  size_t t0 = clock();

  char fname[512];
 
  sprintf(fname, nBase, 1);
  fChain->Add(fname);

  NTupleReader tr(fChain);

  // initialize the type3Ptr defined in the customize.h
  AnaFunctions::prepareTopTagger();

  // The passBaselineFunc is registered here
  tr.registerFunction(&passBaselineFunc);

  //define my AccRecoIsoEffs class to stroe counts and efficiencies
  AccRecoIsoEffs myAccRecoIsoEffs;

  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram("test.root");

  int nevents_muonCS= 0;
  int nevents_baseline= 0;
  int nevents_baseline_ref= 0;

  double mtwcorrfactor[8];
  mtwcorrfactor[0] = 1.05;
  mtwcorrfactor[1] = 1.06;
  mtwcorrfactor[2] = 1.04;
  mtwcorrfactor[3] = 1.11;
  mtwcorrfactor[4] = 1.15;
  mtwcorrfactor[5] = 1.17;
  mtwcorrfactor[6] = 1.20;
  mtwcorrfactor[7] = 1.62;

  double effiso[8];
  effiso[0] = 0.401575;
  effiso[1] = 0.643443;
  effiso[2] = 0.714286;
  effiso[3] = 0.764977;
  effiso[4] = 0.811828;
  effiso[5] = 0.873846;
  effiso[6] = 0.896552;
  effiso[7] = 0.939597;

  double effid[8];
  effid[0] = 0.900709;
  effid[1] = 0.949416;
  effid[2] = 0.929204;
  effid[3] = 0.973094;
  effid[4] = 0.96875;
  effid[5] = 0.964392;
  effid[6] = 0.9631;
  effid[7] = 0.892216;

  double effacc=0.896245;

  int countevt=0;

  while(tr.getNextEvent())
  {
    ++countevt;
    //if (countevt>100) break;
    
    if(tr.getEvtNum()%50000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    
    myAccRecoIsoEffs.nevents_tot++;

    bool passBaseline=tr.getVar<bool>("passBaseline");
    if (passBaseline)
    {
      ++nevents_baseline_ref;
    }

    bool passBaseline_nolepveto=tr.getVar<bool>("passBaseline_nolepveto");
    if (passBaseline_nolepveto)
    {
      myAccRecoIsoEffs.nevents_sel_base++;

      int nElectrons = tr.getVar<int>("nElectrons_CUT2");
      //int nElectrons = tr.getVar<int>("nElectrons_CUT");
      //int nElectrons2 = tr.getVar<int>("nElectrons_CUT2");
      //if (nElectrons != nElectrons2) std::cout << "nElectrons = " << nElectrons  << " , nElectrons2 = " << nElectrons2 << std::endl;

      int nMuons = tr.getVar<int>("nMuons_CUT2");
      vector<double> muonsRelIso = tr.getVec<double>("muonsRelIso");
      const double met=tr.getVar<double>("met");
      const double metphi=tr.getVar<double>("metphi");

      vector<int> W_emuVec = tr.getVec<int>("W_emuVec");
      vector<int> W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
      vector<int> emuVec_merge;
      emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
      emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
      emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
      int gen_emus_count = emuVec_merge.size();

      vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");

      int ngenmu=0;
      int ngenmuoutacc=0;
      int ngenmunotid=0;
      int ngenmunotiso=0;

      if(nElectrons == 0)
      {
        myAccRecoIsoEffs.nevents_sel_mus++;

        int reco_mus_count = muonsLVec.size();

        for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
        {
          //determine if this gen particle is Muon;
          vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          bool isMuon;
          isMuon = false;
          isMuon = ( ( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == 13 )||( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == -13 ) );

          if( isMuon )
          {
	    ++ngenmu;
            myAccRecoIsoEffs.nmus++;
            
            double gen_mus_eta, gen_mus_phi, gen_mus_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

            gen_mus_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_mus_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_mus_pt  = ( genDecayLVec.at ( genId ) ).Pt();

            if((std::abs(gen_mus_eta)) < 2.4 && gen_mus_pt > 5)
            {
              myAccRecoIsoEffs.nmus_acc++;

              int ptbin_number = Set_ptbin_number(gen_mus_pt);

              myAccRecoIsoEffs.nmus_acc_bin[ptbin_number]++;

              //loop over reco lepton information to find out smallest deltaR value
              vector<double> deltar_mus_pool;
              for(int reco_mus_i = 0 ; reco_mus_i < reco_mus_count ; reco_mus_i++)
              {
                double deltar_media;
                deltar_media = DeltaR(gen_mus_eta,
                                      gen_mus_phi,
                                      (muonsLVec.at(reco_mus_i)).Eta(),
                                      (muonsLVec.at(reco_mus_i)).Phi()
                                     );

                deltar_mus_pool.push_back(deltar_media);
              }

              double deltar;
              deltar = 1000;
              int mindeltar_index;

              if(deltar_mus_pool.size() > 0)
              {
                deltar = *(std::min_element(deltar_mus_pool.begin(), deltar_mus_pool.end()));
                mindeltar_index = std::min_element(deltar_mus_pool.begin(), deltar_mus_pool.end()) - deltar_mus_pool.begin();
              }
              //h_b_deltaR_mus->Fill(deltar);
              deltar_mus_pool.clear();

              bool ismatcheddeltaR;
              ismatcheddeltaR = (deltar < 0.02);

              if(ismatcheddeltaR
                 //&& 
                 //isgoodmuonid
                )
              {
                myAccRecoIsoEffs.nmus_reco[ptbin_number]++;

                bool mus_pass_iso;
                mus_pass_iso = false;               
                mus_pass_iso = ( muonsRelIso.at(mindeltar_index) < 0.2 );
                
                if(mus_pass_iso)
                {
                  myAccRecoIsoEffs.nmus_iso[ptbin_number]++;
                }//if isolated
		else
		{
		  ++ngenmunotiso;
		}
              }//if reconstructed
	      else
	      {
		++ngenmunotid;
	      }
            }//if accepted
	    else
	    {
	      ++ngenmuoutacc;
	    }
          }//if the gen particle is muon
        }//loop gen electrons/muons
      }//if no electrons

      if(nMuons == 0)
      {
        myAccRecoIsoEffs.nevents_sel_els++;

        vector<TLorentzVector> elesLVec = tr.getVec<TLorentzVector>("elesLVec");
        int reco_els_count = elesLVec.size();

        for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
        {
          //determine if this gen particle is electrons;
          vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          bool isElectron;
          isElectron = false;
          isElectron = ( ( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == 11 )||( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == -11 ) );
          
          if( isElectron )
          {
            myAccRecoIsoEffs.nels++;
          
            double gen_els_eta, gen_els_phi, gen_els_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          
            gen_els_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_els_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_els_pt  = ( genDecayLVec.at ( genId ) ).Pt();
    
            if((std::abs(gen_els_eta)) < 2.5 && gen_els_pt > 5)
            {
              myAccRecoIsoEffs.nels_acc++;

              int ptbin_number = Set_ptbin_number(gen_els_pt);

              myAccRecoIsoEffs.nels_acc_bin[ptbin_number]++;
       
              //loop over reco lepton information to determine the smallest deltar
              vector<double> deltar_els_pool;
              for(int reco_els_i = 0 ; reco_els_i < reco_els_count ; reco_els_i++)
              {
                double deltar_media;
                deltar_media = DeltaR(gen_els_eta,
                                      gen_els_phi,
                                      (elesLVec.at(reco_els_i)).Eta(),
                                      (elesLVec.at(reco_els_i)).Phi()
                                     );

                deltar_els_pool.push_back(deltar_media);
              }

              double deltar;
              deltar = 1000;
              int mindeltar_index;

              if(deltar_els_pool.size() > 0)
              {
                deltar = *(std::min_element(deltar_els_pool.begin(), deltar_els_pool.end()));
                mindeltar_index = std::min_element(deltar_els_pool.begin(), deltar_els_pool.end()) - deltar_els_pool.begin();
              }
              //h_b_deltaR_els->Fill(deltar);  
              deltar_els_pool.clear();

              bool ismatcheddeltaR;
              ismatcheddeltaR = (deltar < 0.02);
    
              if(ismatcheddeltaR
                 //&& 
                 //isgoodmuonid
                )
              {
                myAccRecoIsoEffs.nels_reco[ptbin_number]++;

                vector<double> elesRelIso = tr.getVec<double>("elesRelIso");
                bool els_pass_iso;
                els_pass_iso = false;
                els_pass_iso = ( elesRelIso.at(mindeltar_index) < 0.24 );

                if(els_pass_iso)
                {
                  myAccRecoIsoEffs.nels_iso[ptbin_number]++;
                }//if isolated
              }//if reconstructed
            }//if accepted
          }//if the gen particle is electron 
        }//loop gen electrons/muons
      }//if no muons

      const int njets30 = tr.getVar<int>("cntNJetsPt30Eta24");
      const double MT2 = tr.getVar<double>("MT22");
      const double bestTopJetMass = tr.getVar<double>("bestTopJetMass2");

      ///////////////////////////
      // expectation computation
      ///////////////////////////

      // exp 1 muon not iso
      if (nElectrons == 0 && nMuons==0 && ngenmu==1 && ngenmunotiso==1)
      {
	(myBaseHistgram.h_exp_mu_iso_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_iso_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_iso_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_iso_topmass)->Fill(bestTopJetMass);
      }

      // exp 1 muon not id
      if (nElectrons == 0 && nMuons==0 && ngenmu==1 && ngenmunotid==1)
      {
	(myBaseHistgram.h_exp_mu_id_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_id_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_id_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_id_topmass)->Fill(bestTopJetMass);
      }

      // exp 1 muon out acc
      if (nElectrons == 0 && nMuons==0 && ngenmu==1 && ngenmuoutacc==1)
      {
	(myBaseHistgram.h_exp_mu_acc_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_acc_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_acc_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_acc_topmass)->Fill(bestTopJetMass);
      }

      // exp 1 muon tot
      if (nElectrons == 0 && nMuons==0 && ngenmu==1 && (ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1))
      {
	(myBaseHistgram.h_exp_mu_all_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_all_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_all_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_all_topmass)->Fill(bestTopJetMass);
      }

     if (nElectrons == 0 && nMuons == 0)
      {
	++nevents_baseline;
      }

      // muon CS
     if (nElectrons == 0 && nMuons == 1)
      {
	// get muon variables
	double muon_pt5=0.0;
	double muon_phi5=0.0;
	for(unsigned int im=0; im<muonsLVec.size(); im++){
         double permuonpt = muonsLVec[im].Pt(), permuoneta = muonsLVec[im].Eta();
         if(fabs(permuoneta) < 2.4 && muonsRelIso[im] < 0.2 )
	 {
	   muon_pt5=permuonpt;
	   muon_phi5=muonsLVec[im].Phi();
	 }
	}

	double deltaphi5=muon_phi5-metphi;
	while (deltaphi5 > M_PI) deltaphi5 -= 2*M_PI;
	while (deltaphi5 <= -M_PI) deltaphi5 += 2*M_PI;
	const double mtW5=std::sqrt(2.0*muon_pt5*met*(1.0-cos(deltaphi5)));
 
	if (mtW5<100.0)
	{
	  ++nevents_muonCS;
    
	  //////////////////////////
	  // prediction computation
	  //////////////////////////

	  double EventWeight=1.0;
	  double EventWeight_iso=1.0;
	  double EventWeight_reco=1.0;
	  double EventWeight_acc=1.0;
    
	  // not iso
	  if (muon_pt5>=5.0 && muon_pt5<=10.0) EventWeight_iso=(1.0-effiso[0])/effiso[0];
	  if (muon_pt5>10.0 && muon_pt5<=20.0) EventWeight_iso=(1.0-effiso[1])/effiso[1];
	  if (muon_pt5>20.0 && muon_pt5<=30.0) EventWeight_iso=(1.0-effiso[2])/effiso[2];
	  if (muon_pt5>30.0 && muon_pt5<=40.0) EventWeight_iso=(1.0-effiso[3])/effiso[3];
	  if (muon_pt5>40.0 && muon_pt5<=50.0) EventWeight_iso=(1.0-effiso[4])/effiso[4];
	  if (muon_pt5>50.0 && muon_pt5<=70.0) EventWeight_iso=(1.0-effiso[5])/effiso[5];
	  if (muon_pt5>70.0 && muon_pt5<=100.0) EventWeight_iso=(1.0-effiso[6])/effiso[6];
	  if (muon_pt5>100.0) EventWeight_iso=(1.0-effiso[7])/effiso[7];

	  // mot reco
	  if (muon_pt5>=5.0 && muon_pt5<=10.0) EventWeight_reco=(1.0-effid[0])/effid[0]/effiso[0];
	  if (muon_pt5>10.0 && muon_pt5<=20.0) EventWeight_reco=(1.0-effid[1])/effid[1]/effiso[1];
	  if (muon_pt5>20.0 && muon_pt5<=30.0) EventWeight_reco=(1.0-effid[2])/effid[2]/effiso[2];
	  if (muon_pt5>30.0 && muon_pt5<=40.0) EventWeight_reco=(1.0-effid[3])/effid[3]/effiso[3];
	  if (muon_pt5>40.0 && muon_pt5<=50.0) EventWeight_reco=(1.0-effid[4])/effid[4]/effiso[4];
	  if (muon_pt5>50.0 && muon_pt5<=70.0) EventWeight_reco=(1.0-effid[5])/effid[5]/effiso[5];
	  if (muon_pt5>70.0 && muon_pt5<=100.0) EventWeight_reco=(1.0-effid[6])/effid[6]/effiso[6];
	  if (muon_pt5>100.0) EventWeight_reco=(1.0-effid[7])/effid[7]/effiso[7];

	  // out of acceptance
	  if (muon_pt5>=5.0 && muon_pt5<=10.0) EventWeight_acc=1.0/effiso[0]/effid[0]*(1.0-effacc)/effacc;
	  if (muon_pt5>10.0 && muon_pt5<=20.0) EventWeight_acc=1.0/effiso[1]/effid[1]*(1.0-effacc)/effacc;
	  if (muon_pt5>20.0 && muon_pt5<=30.0) EventWeight_acc=1.0/effiso[2]/effid[2]*(1.0-effacc)/effacc;
	  if (muon_pt5>30.0 && muon_pt5<=40.0) EventWeight_acc=1.0/effiso[3]/effid[3]*(1.0-effacc)/effacc;
	  if (muon_pt5>40.0 && muon_pt5<=50.0) EventWeight_acc=1.0/effiso[4]/effid[4]*(1.0-effacc)/effacc;
	  if (muon_pt5>50.0 && muon_pt5<=70.0) EventWeight_acc=1.0/effiso[5]/effid[5]*(1.0-effacc)/effacc;
	  if (muon_pt5>70.0 && muon_pt5<=100.0) EventWeight_acc=1.0/effiso[6]/effid[6]*(1.0-effacc)/effacc;
	  if (muon_pt5>100.0) EventWeight_acc=1.0/effiso[7]/effid[7]*(1.0-effacc)/effacc;

	  // mtwcorrfactor
	  if (muon_pt5>=5.0 && muon_pt5<=10.0) EventWeight*=mtwcorrfactor[0];
	  if (muon_pt5>10.0 && muon_pt5<=20.0) EventWeight*=mtwcorrfactor[1];
	  if (muon_pt5>20.0 && muon_pt5<=30.0) EventWeight*=mtwcorrfactor[2];
	  if (muon_pt5>30.0 && muon_pt5<=40.0) EventWeight*=mtwcorrfactor[3];
	  if (muon_pt5>40.0 && muon_pt5<=50.0) EventWeight*=mtwcorrfactor[4];
	  if (muon_pt5>50.0 && muon_pt5<=70.0) EventWeight*=mtwcorrfactor[5];
	  if (muon_pt5>70.0 && muon_pt5<=100.0) EventWeight*=mtwcorrfactor[6];
	  if (muon_pt5>100.0) EventWeight*=mtwcorrfactor[7];

	  // dimuon correction factor
	  // need to be added!!!

	  // Fill muon iso closure plots
	  (myBaseHistgram.h_pred_mu_iso_met)->Fill(met, EventWeight_iso*EventWeight);
	  (myBaseHistgram.h_pred_mu_iso_njets)->Fill(njets30, EventWeight_iso*EventWeight);
	  (myBaseHistgram.h_pred_mu_iso_mt2)->Fill(MT2, EventWeight_iso*EventWeight);
	  (myBaseHistgram.h_pred_mu_iso_topmass)->Fill(bestTopJetMass, EventWeight_iso*EventWeight);
	  // Fill muon id closure plots
	  (myBaseHistgram.h_pred_mu_id_met)->Fill(met, EventWeight_reco*EventWeight);
	  (myBaseHistgram.h_pred_mu_id_njets)->Fill(njets30, EventWeight_reco*EventWeight);
	  (myBaseHistgram.h_pred_mu_id_mt2)->Fill(MT2, EventWeight_reco*EventWeight);
	  (myBaseHistgram.h_pred_mu_id_topmass)->Fill(bestTopJetMass, EventWeight_reco*EventWeight);
	  // Fill muon acc closure plots
	  (myBaseHistgram.h_pred_mu_acc_met)->Fill(met, EventWeight_acc*EventWeight);
	  (myBaseHistgram.h_pred_mu_acc_njets)->Fill(njets30, EventWeight_acc*EventWeight);
	  (myBaseHistgram.h_pred_mu_acc_mt2)->Fill(MT2, EventWeight_acc*EventWeight);
	  (myBaseHistgram.h_pred_mu_acc_topmass)->Fill(bestTopJetMass, EventWeight_acc*EventWeight);
	  // Fill all muon closure plots
	  (myBaseHistgram.h_pred_mu_all_met)->Fill(met, (EventWeight_iso+EventWeight_reco+EventWeight_acc)*EventWeight);
	  (myBaseHistgram.h_pred_mu_all_njets)->Fill(njets30, (EventWeight_iso+EventWeight_reco+EventWeight_acc)*EventWeight);
	  (myBaseHistgram.h_pred_mu_all_mt2)->Fill(MT2, (EventWeight_iso+EventWeight_reco+EventWeight_acc)*EventWeight);
	  (myBaseHistgram.h_pred_mu_all_topmass)->Fill(bestTopJetMass, (EventWeight_iso+EventWeight_reco+EventWeight_acc)*EventWeight);

  	} // mtW5<100.0 (muon CS)
      } // nElectrons == 0 && nMuons == 1
    }
    
    //<< "\t" << tr.getVar<double>("joe") << "\t" << tr.getVar<int>("five") << "\t" << tr.getVec<double>("muonsMtw").size() << "\t" << tr.getVec<double>("threeNum")[2] << std::endl;
  }

  myAccRecoIsoEffs.printOverview();

  myAccRecoIsoEffs.NumberstoEffs();

  myAccRecoIsoEffs.printAccRecoIsoEffs();

  (myBaseHistgram.oFile)->Write();

  const double ttbarCrossSection=806.1;
  const double lumi=1000.0;
  const double ntoteventsttbar=25446993.0;

  std::cout << "nevents_muonCS = " << nevents_muonCS << std::endl;
  std::cout << "nevents_muonCS_norm (1fb-1) = " << nevents_muonCS*ttbarCrossSection*lumi/ntoteventsttbar << std::endl;

  std::cout << "nevents_baseline = " << nevents_baseline << std::endl;
  std::cout << "nevents_baseline_ref = " << nevents_baseline_ref << std::endl;
  std::cout << "nevents_baseline_norm (1fb-1) = " << nevents_baseline*ttbarCrossSection*lumi/ntoteventsttbar << std::endl;

  return 0;
}


void AccRecoIsoEffs::printOverview()
{
  std::cout << "Lost Lepton Sample(ttbar sample) Overview:" << std::endl;
  std::cout << "nevents_tot = " << nevents_tot << std::endl;
  std::cout << "nevents_sel_base: = " << nevents_sel_base << std::endl;

  //const double ttbarCrossSection=818.8;
  //const double lumi=1000.0;
  //const double ntoteventsttbar=17184825.0;
  //const double ntoteventsttbar=3029767.0;
  
  std::cout << "nevents_sel_mus = " << nevents_sel_mus << std::endl;
  std::cout << "nevents_sel_els = " << nevents_sel_els << std::endl;

  return ;
}

void AccRecoIsoEffs::NumberstoEffs()
{
  mus_acc=nmus_acc/nmus;
  els_acc=nels_acc/nels;

  mus_acc_err = std::sqrt( get_stat_Error(nmus_acc,nmus)*get_stat_Error(nmus_acc,nmus) + get_sys_Error(mus_acc,0.09)*get_sys_Error(mus_acc,0.09) );
  els_acc_err = std::sqrt( get_stat_Error(nels_acc,nels)*get_stat_Error(nels_acc,nels) + get_sys_Error(els_acc,0.09)*get_sys_Error(els_acc,0.09) ); 

  //to be determined??
  int i_cal;
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    mus_recoeff[i_cal]=nmus_reco[i_cal]/nmus_acc_bin[i_cal];
    mus_recoeff_err[i_cal]=get_stat_Error(nmus_reco[i_cal],nmus_acc_bin[i_cal]);
    els_recoeff[i_cal]=nels_reco[i_cal]/nels_acc_bin[i_cal];
    els_recoeff_err[i_cal]=get_stat_Error(nels_reco[i_cal],nels_acc_bin[i_cal]);


    mus_isoeff[i_cal]=nmus_iso[i_cal]/nmus_reco[i_cal];
    mus_isoeff_err[i_cal]=get_stat_Error(nmus_iso[i_cal],nmus_reco[i_cal]);
    els_isoeff[i_cal]=nels_iso[i_cal]/nels_reco[i_cal];
    els_isoeff_err[i_cal]=get_stat_Error(nels_iso[i_cal],nels_reco[i_cal]);
  }
  
  return ;
}

void AccRecoIsoEffs::printAccRecoIsoEffs()
{
  int i_cal = 0;
  
  std::cout<<std::endl<<"Muon information: "<<std::endl;

  std::cout<<"number of muons from top: "<<nmus<<std::endl;

  std::cout<<"number of muons from top, accepted: "<<nmus_acc<<std::endl;

  std::cout<<"number of muons from top, accepted, bins: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nmus_acc_bin[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"number of muons from top, reconstructed: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nmus_reco[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"number of muons from top, isolated: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nmus_iso[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"muons from top, acceptance: "<<mus_acc<<"("<<mus_acc_err<<")"<<std::endl;

  std::cout<<"muons from top, reconstruction efficiency: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<mus_recoeff[i_cal]<<"("<<mus_recoeff_err[i_cal]<<")"<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"muons from top, isolation efficiency: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<mus_isoeff[i_cal]<<"("<<mus_isoeff_err[i_cal]<<")"<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }
 
  std::cout<<std::endl<<"Electron information: "<<std::endl;

  std::cout<<"number of electrons from top: "<<nels<<std::endl;

  std::cout<<"number of electrons from top, accepted: "<<nels_acc<<std::endl;

  std::cout<<"number of electrons from top, accepted, bins: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nels_acc_bin[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"number of electrons from top, reconstructed: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nels_reco[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"number of electrons from top, isolated: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<nels_iso[i_cal]<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"electrons from top, acceptance: "<<els_acc<<"("<<els_acc_err<<")"<<std::endl;

  std::cout<<"electrons from top, reconstruction efficiency: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<els_recoeff[i_cal]<<"("<<els_recoeff_err[i_cal]<<")"<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  std::cout<<"electrons from top, isolation efficiency: ";
  for(i_cal=0;i_cal<PT_BINS;i_cal++)
  {
    std::cout<<els_isoeff[i_cal]<<"("<<els_isoeff_err[i_cal]<<")"<<" ";
    if(i_cal==PT_BINS-1)
    {
      std::cout<<std::endl;
    }
  }

  return ;
}

