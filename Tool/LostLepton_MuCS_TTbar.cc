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

using namespace std;

//const bool debug = true;

void passBaselineFunc(NTupleReader &tr)
{
  bool passBaseline = true;
  bool passBaseline_nolepveto = true;

  //Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>("met"), 0, tr.getVar<double>("metphi"), 0);

  //Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsRelIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesRelIso"), tr.getVec<double>("elesMtw"), AnaConsts::elesArr);
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
  bool passLeptVeto = true;
  if( nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passLeptVeto = false; }
  if( nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passLeptVeto = false; }
  //Isolated track veto is disabled for now
  //if( nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passLeptVeto = false; }
  //if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passBaseline_nolepveto = false; passnJets = false;}
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passBaseline_nolepveto = false; passnJets = false;}
  //if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passBaseline_nolepveto = false; passdPhis = false; }
  //if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBaseline_nolepveto = false; passBJets = false; }
  //if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

  //Pass the baseline MET requirement?
  bool passMET = true;
  if( tr.getVar<double>("met") < AnaConsts::defaultMETcut ){ passBaseline = false; passBaseline_nolepveto = false; passMET = false; }
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
  if( passBaseline_nolepveto && cntNJetsPt30 >= AnaConsts::nJetsSel )
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
  if( bestTopJetIdx == -1 ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
  if( ! remainPassCSVS ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
  if( pickedRemainingCombfatJetIdx == -1 && jetsLVec_forTagger->size()>=6 ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
  if( ! (bestTopJetMass > AnaConsts::lowTopCut_ && bestTopJetMass < AnaConsts::highTopCut_ ) ){ passBaseline = false; passBaseline_nolepveto = false; passTagger = false; }
  //if( debug ) std::cout<<"bestTopJetidx : "<<bestTopJetIdx<<"  remainPassCSVS : "<<remainPassCSVS<<"  pickedRemainingCombfatJetIdx : "<<pickedRemainingCombfatJetIdx<<"  bestTopJetMass : "<<bestTopJetMass<<"  passBaseline : "<<passBaseline<<std::endl;

  //Register all the calculated variables
  tr.registerDerivedVar("nMuons_CUT2", nMuons);
  tr.registerDerivedVar("nElectrons_CUT2", nElectrons);
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
  tr.registerDerivedVar("passBaseline", passBaseline);
  tr.registerDerivedVar("passBaseline_nolepveto", passBaseline_nolepveto);

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

  TChain *fChain = new TChain("AUX");

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
  //define my AccRecoIsoEffs class to stroe counts and efficiencies
  AccRecoIsoEffs myAccRecoIsoEffs;
  //define my histgram class
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  int nevents_muonCS= 0;
  int nevents_baseline= 0;
  //int nevents_baseline_ref= 0;

  double mtwcorrfactor[8];
  mtwcorrfactor[0] = 1.01;
  mtwcorrfactor[1] = 1.04;
  mtwcorrfactor[2] = 1.07;
  mtwcorrfactor[3] = 1.20;
  mtwcorrfactor[4] = 1.15;
  mtwcorrfactor[5] = 1.20;
  mtwcorrfactor[6] = 1.28;
  mtwcorrfactor[7] = 1.83;

  //first loop, to generate Acc, reco and Iso effs and also fill expected histgram
  std::cout<<"First loop begin: "<<std::endl;
  while(tr.getNextEvent())
  {
    if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    
    myAccRecoIsoEffs.nevents_tot++;

    //bool passBaseline=tr.getVar<bool>("passBaseline");
    //if (passBaseline)
    //{
      //++nevents_baseline_ref;
    //}

    bool passBaseline_nolepveto=tr.getVar<bool>("passBaseline_nolepveto");
    if (passBaseline_nolepveto)
    {
      myAccRecoIsoEffs.nevents_sel_base++;

      //nMuons in flatree means no iso cut muons; CUT2 we add iso
      int nElectrons = tr.getVar<int>("nElectrons_CUT2");
      int nMuons = tr.getVar<int>("nMuons_CUT2");

      double met=tr.getVar<double>("met");
      double metphi=tr.getVar<double>("metphi");

      vector<int> W_emuVec = tr.getVec<int>("W_emuVec");
      vector<int> W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
      vector<int> emuVec_merge;
      emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
      emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
      emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
      int gen_emus_count = emuVec_merge.size();

      int ngenmunotiso = 0, ngenmunotid = 0, ngenmuoutacc = 0;
      int ngenelnotiso = 0, ngenelnotid = 0, ngeneloutacc = 0;

      int ngenmu = 0;
      int ngenel = 0;

      if(nElectrons == 0)
      {
        myAccRecoIsoEffs.nevents_sel_mus++;
       
        //get reco level information of muons
        vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
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
            ngenmu++;
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

                vector<double> muonsRelIso = tr.getVec<double>("muonsRelIso");

                bool mus_pass_iso;
                mus_pass_iso = false;               
                mus_pass_iso = ( muonsRelIso.at(mindeltar_index) < 0.2 );
                
                if(mus_pass_iso)
                {
                  myAccRecoIsoEffs.nmus_iso[ptbin_number]++;
                }//if isolated
                else
                {
                  ngenmunotiso++;
                }
              }//if reconstructed
              else
              {
                ngenmunotid++;
              }
            }//if accepted
            else
            {
              ngenmuoutacc++;
            }
          }//if the gen particle is muon
        }//loop gen electrons/muons
      }//if no electrons

      if(nMuons == 0)
      {
        myAccRecoIsoEffs.nevents_sel_els++;

        //get reco level information of electrons
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
            ngenel++;
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
                else
                {
                  ngenelnotiso++;
                }
              }//if reconstructed
              else
              {
                ngenelnotid++;
              }
            }//if accepted
            else
            {
              ngeneloutacc++;
            }
          }//if the gen particle is electron 
        }//loop gen electrons/muons
      }//if no muons
      
      //book the variable we need to check
      int njets30 = tr.getVar<int>("cntNJetsPt30Eta24");
      int ntopjets = tr.getVar<int>("nTopCandSortedCnt");
      double MT2 = tr.getVar<double>("MT22");
      double bestTopJetMass = tr.getVar<double>("bestTopJetMass2");
      double ht = tr.getVar<double>("ht");
      double mht = tr.getVar<double>("mht");

      ///////////////////////////
      // expectation computation
      ///////////////////////////

      // exp 1 muon not iso
      if ( nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmunotiso==1) || ( ngenmu==2 && ngenmunotiso==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_iso_mus++; 
    
	(myBaseHistgram.h_exp_mu_iso_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_iso_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_iso_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_iso_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_mu_iso_ht)->Fill(ht);
        (myBaseHistgram.h_exp_mu_iso_mht)->Fill(mht);
        (myBaseHistgram.h_exp_mu_iso_ntopjets)->Fill(ntopjets);
      }

      // exp 1 muon not id
      if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmunotid==1) || ( ngenmu==2 && ngenmunotid==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_id_mus++;

	(myBaseHistgram.h_exp_mu_id_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_id_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_id_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_id_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_mu_id_ht)->Fill(ht);
        (myBaseHistgram.h_exp_mu_id_mht)->Fill(mht);
        (myBaseHistgram.h_exp_mu_id_ntopjets)->Fill(ntopjets);
      }

      // exp 1 muon out acc
      if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmuoutacc==1) || ( ngenmu==2 && ngenmuoutacc==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_acc_mus++;

	(myBaseHistgram.h_exp_mu_acc_met)->Fill(met);
	(myBaseHistgram.h_exp_mu_acc_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_mu_acc_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_mu_acc_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_mu_acc_ht)->Fill(ht);
        (myBaseHistgram.h_exp_mu_acc_mht)->Fill(mht);
        (myBaseHistgram.h_exp_mu_acc_ntopjets)->Fill(ntopjets);
      }

      // exp 1 muon tot
      if (nElectrons == 0 && nMuons==0 && ngenmu==1 && (ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1))
      //if (nElectrons == 0 && nMuons==0 && ngenmu==1)
      {
        myAccRecoIsoEffs.nevents_exp_all_mus++;
        myAccRecoIsoEffs.nevents_single_mus++;
	
        (myBaseHistgram.h_exp_musingle_all_met)->Fill(met);
	(myBaseHistgram.h_exp_musingle_all_njets)->Fill(njets30);
	(myBaseHistgram.h_exp_musingle_all_mt2)->Fill(MT2);
	(myBaseHistgram.h_exp_musingle_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_musingle_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_musingle_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_musingle_all_ntopjets)->Fill(ntopjets);

        (myBaseHistgram.h_exp_mu_all_met)->Fill(met);
        (myBaseHistgram.h_exp_mu_all_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_mu_all_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_mu_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_mu_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_mu_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets);
      }

      if ( nElectrons == 0 && nMuons==0 && ngenmu==2 && ( ngenmuoutacc==2 || ngenmunotid==2 || ngenmunotiso==2 || ( ngenmuoutacc==1 && ngenmunotid==1 ) || (ngenmuoutacc==1 && ngenmunotiso==1 ) || ( ngenmunotiso==1 && ngenmunotid==1 ) ) )
      //if ( nElectrons == 0 && nMuons==0 && ngenmu==2 )
      {
        myAccRecoIsoEffs.nevents_di_mus++;
        myAccRecoIsoEffs.nevents_exp_all_mus++;

        (myBaseHistgram.h_exp_mu_all_met)->Fill(met);
        (myBaseHistgram.h_exp_mu_all_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_mu_all_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_mu_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_mu_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_mu_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets);
      }

      // exp 1 electron not iso
      if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotiso==1) || ( ngenel==2 && ngenelnotiso==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_iso_els++;

        (myBaseHistgram.h_exp_el_iso_met)->Fill(met);
        (myBaseHistgram.h_exp_el_iso_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_el_iso_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_el_iso_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_el_iso_ht)->Fill(ht);
        (myBaseHistgram.h_exp_el_iso_mht)->Fill(mht);
        (myBaseHistgram.h_exp_el_iso_ntopjets)->Fill(ntopjets);
      }

      // exp 1 electron not id
      if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotid==1) || ( ngenel==2 && ngenelnotid==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_id_els++;

        (myBaseHistgram.h_exp_el_id_met)->Fill(met);
        (myBaseHistgram.h_exp_el_id_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_el_id_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_el_id_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_el_id_ht)->Fill(ht);
        (myBaseHistgram.h_exp_el_id_mht)->Fill(mht);
        (myBaseHistgram.h_exp_el_id_ntopjets)->Fill(ntopjets);
      }

      // exp 1 electron not acc
      if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngeneloutacc==1) || ( ngenel==2 && ngeneloutacc==2) ) )
      {
        myAccRecoIsoEffs.nevents_exp_acc_els++;

        (myBaseHistgram.h_exp_el_acc_met)->Fill(met);
        (myBaseHistgram.h_exp_el_acc_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_el_acc_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_el_acc_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_el_acc_ht)->Fill(ht);
        (myBaseHistgram.h_exp_el_acc_mht)->Fill(mht);
        (myBaseHistgram.h_exp_el_acc_ntopjets)->Fill(ntopjets);
      }

      // exp 1 electron tot
      if (nElectrons == 0 && nMuons==0 && ngenel==1 && (ngeneloutacc==1 || ngenelnotid==1 || ngenelnotiso==1))
      //if (nElectrons == 0 && nMuons==0 && ngenel==1)
      {
        myAccRecoIsoEffs.nevents_exp_all_els++;
        myAccRecoIsoEffs.nevents_single_els++;

        (myBaseHistgram.h_exp_elsingle_all_met)->Fill(met);
        (myBaseHistgram.h_exp_elsingle_all_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_elsingle_all_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_elsingle_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_elsingle_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_elsingle_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_elsingle_all_ntopjets)->Fill(ntopjets);  

        (myBaseHistgram.h_exp_el_all_met)->Fill(met);
        (myBaseHistgram.h_exp_el_all_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_el_all_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_el_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_el_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_el_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets);
      }

      if ( nElectrons == 0 && nMuons==0 && ngenel==2 && ( ngeneloutacc==2 || ngenelnotid==2 || ngenelnotiso==2 || ( ngeneloutacc==1 && ngenelnotid==1 ) || (ngeneloutacc==1 && ngenelnotiso==1 ) || ( ngenelnotiso==1 && ngenelnotid==1 ) ) )
      //if (nElectrons == 0 && nMuons==0 && ngenel==2)
      {
        myAccRecoIsoEffs.nevents_exp_all_els++;
        myAccRecoIsoEffs.nevents_di_els++;

        (myBaseHistgram.h_exp_el_all_met)->Fill(met);
        (myBaseHistgram.h_exp_el_all_njets)->Fill(njets30);
        (myBaseHistgram.h_exp_el_all_mt2)->Fill(MT2);
        (myBaseHistgram.h_exp_el_all_topmass)->Fill(bestTopJetMass);
        (myBaseHistgram.h_exp_el_all_ht)->Fill(ht);
        (myBaseHistgram.h_exp_el_all_mht)->Fill(mht);
        (myBaseHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets);
      }
    }//baseline, nolepveto
  }//end of first loop

  //All numbers counted, now calculated effs and print out 
  myAccRecoIsoEffs.printOverview();
  myAccRecoIsoEffs.NumberstoEffs();
  myAccRecoIsoEffs.EffstoWeights();
  myAccRecoIsoEffs.GetDiLeptonFactor();
  myAccRecoIsoEffs.printAccRecoIsoEffs();

  NTupleReader trCS(fChain);
  //initialize the type3Ptr defined in the customize.h
  AnaFunctions::prepareTopTagger();
  //The passBaselineFunc is registered here
  trCS.registerFunction(&passBaselineFunc);

  //second loop, to select CS sample and make prediction
  std::cout<<"Second loop begin: "<<std::endl;
  while(trCS.getNextEvent())
  {
    if(trCS.getEvtNum()%20000 == 0) std::cout << trCS.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

    bool passBaseline_nolepveto=trCS.getVar<bool>("passBaseline_nolepveto");

    if(passBaseline_nolepveto)
    {
      int nElectrons = trCS.getVar<int>("nElectrons_CUT2");
      int nMuons = trCS.getVar<int>("nMuons_CUT2");

      double met = trCS.getVar<double>("met");
      double metphi = trCS.getVar<double>("metphi");

      int njets30 = trCS.getVar<int>("cntNJetsPt30Eta24");
      int ntopjets = tr.getVar<int>("nTopCandSortedCnt");
      double MT2 = trCS.getVar<double>("MT22");
      double bestTopJetMass = trCS.getVar<double>("bestTopJetMass2");
      double ht = trCS.getVar<double>("ht");
      double mht = trCS.getVar<double>("mht");

      //muon CS
      //nMuons in flatree means no iso cut muons; CUT2 we add iso
      if (nElectrons == 0 && nMuons == 1)
      {
        //get muon variables
	vector<TLorentzVector> muonsLVec = trCS.getVec<TLorentzVector>("muonsLVec");
        vector<double> muonsRelIso = trCS.getVec<double>("muonsRelIso");

        double reco_mus_pt = 0, reco_mus_eta = 0, reco_mus_phi = 0;

        for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
        {
          if( fabs(muonsLVec[im].Eta()) < 2.4 && muonsRelIso[im] < 0.2 )
	  {
            reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
            reco_mus_eta = ( muonsLVec.at(im) ).Eta();
            reco_mus_phi = ( muonsLVec.at(im) ).Phi();
	  }
	}

        double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
        double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * (1.0 - cos(deltaphi_mus) ) );

        if ( mtW_mus < 100.0 )
        {
	  //////////////////////////
	  // prediction computation
	  //////////////////////////
	  double EventWeight_mus = 1.0;
          int ptbin_number = Set_ptbin_number(reco_mus_pt);

	  //mtwcorrfactor
	  EventWeight_mus = EventWeight_mus * mtwcorrfactor[ptbin_number];
	  //dimuon correction factor
          EventWeight_mus = EventWeight_mus * myAccRecoIsoEffs.corrfactor_di_mus;

          //muon prediction from muon CS
	  //Fill muon iso closure plots
	  (myBaseHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus);
	  //Fill muon id closure plots
	  (myBaseHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus);
	  //Fill muon acc closure plots
	  (myBaseHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus);
	  //Fill all muon closure plots
	  double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number] + myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number] + myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number];

	  (myBaseHistgram.h_pred_mu_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_topmass)->Fill(bestTopJetMass, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

          //total events flow for muons, prediction
          myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus;
        
          myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[ptbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[ptbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[ptbin_number]*EventWeight_mus;

          //begin to predict lost electrons from muon CS
          double EventWeight_els = 1.0;
          //mtwcorrfactor
          EventWeight_els = EventWeight_els * mtwcorrfactor[ptbin_number];
          //dielectron correction factor
          EventWeight_els = EventWeight_els * myAccRecoIsoEffs.corrfactor_di_els;
  
          //electron prediction from muon CS
          //Fill electron iso closure plots
          (myBaseHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els);
          //Fill electron id closure plots
          (myBaseHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els);
          //Fill electron acc closure plots
          (myBaseHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els);
          //Fill all electron closure plots
          double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number] + myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number] + myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number];

          (myBaseHistgram.h_pred_el_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_topmass)->Fill(bestTopJetMass, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);
          //total events flow for electrons, prediction
          myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els;

          myAccRecoIsoEffs.nevents_pred_all_els_err += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[ptbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[ptbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[ptbin_number]*EventWeight_els;
        }//mtW5_els<100 (muon CS)
      }//nElectrons == 0 && nMuons == 1 (muon CS)
    }//baseline_nolepveto
  }

  myAccRecoIsoEffs.NormalizeFlowNumber();
  myAccRecoIsoEffs.printNormalizeFlowNumber();

  //write into histgram
  (myBaseHistgram.oFile)->Write();
  
  //const double ttbarCrossSection=806.1;
  //const double lumi=1000.0;
  //const double ntoteventsttbar=25446993.0;
  //std::cout << "nevents_muonCS = " << nevents_muonCS << std::endl;
  //std::cout << "nevents_muonCS_norm (1fb-1) = " << nevents_muonCS*ttbarCrossSection*lumi/ntoteventsttbar << std::endl;
  //std::cout << "nevents_baseline = " << nevents_baseline << std::endl;
  //std::cout << "nevents_baseline_ref = " << nevents_baseline_ref << std::endl;
  //std::cout << "nevents_baseline_norm (1fb-1) = " << nevents_baseline*ttbarCrossSection*lumi/ntoteventsttbar << std::endl;
  return 0;
}


void AccRecoIsoEffs::printOverview()
{
  std::cout << "Lost Lepton Sample(ttbar sample) Overview:" << std::endl;
  std::cout << "nevents_tot = " << nevents_tot << std::endl;
  std::cout << "nevents_sel_base: = " << nevents_sel_base << std::endl;

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

  int i_cal;
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
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

void AccRecoIsoEffs::EffstoWeights()
{
  int i_cal;
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    mus_EventWeight_iso[i_cal]  = (1.0 - mus_isoeff[i_cal])/mus_isoeff[i_cal];
    mus_EventWeight_reco[i_cal] = (1.0 - mus_recoeff[i_cal])/mus_recoeff[i_cal]/mus_isoeff[i_cal]; 
    mus_EventWeight_acc[i_cal]  = 1.0/mus_isoeff[i_cal]/mus_recoeff[i_cal]*(1.0 - mus_acc)/mus_acc;
    
    els_EventWeight_iso[i_cal]  = 1.0/mus_isoeff[i_cal]/mus_recoeff[i_cal]*(1.0 - els_acc)/mus_acc;;
    els_EventWeight_reco[i_cal] = 1.0/mus_isoeff[i_cal]*(1 - els_recoeff[i_cal])/mus_recoeff[i_cal]*els_acc/mus_acc; 
    els_EventWeight_acc[i_cal]  = (1.0 - els_isoeff[i_cal])/mus_isoeff[i_cal]*els_recoeff[i_cal]/mus_recoeff[i_cal]*els_acc/mus_acc;
  }

  return ;
}

void AccRecoIsoEffs::GetDiLeptonFactor()
{
  corrfactor_di_mus = ( nevents_single_mus + nevents_di_mus ) / ( nevents_single_mus + 2 * nevents_di_mus );
  corrfactor_di_els = ( nevents_single_els + nevents_di_els ) / ( nevents_single_els + 2 * nevents_di_els );

  double alpha;
  alpha = 1-0.6827;

  double mus1dev = 0, mus2dev = 0;
  mus1dev = nevents_di_mus/(nevents_single_mus + 2 * nevents_di_mus)/(nevents_single_mus + 2 * nevents_di_mus);
  mus2dev = (0-nevents_single_mus)/(nevents_single_mus + 2 * nevents_di_mus)/(nevents_single_mus + 2 * nevents_di_mus);

  corrfactor_di_mus_err = std::sqrt(mus1dev*mus1dev*ROOT::Math::gamma_quantile_c( alpha/2 , nevents_single_mus + 1, 1 ) + mus2dev*mus2dev*ROOT::Math::gamma_quantile_c( alpha/2 , nevents_di_mus + 1 , 1 ));

  double els1dev = 0, els2dev = 0;
  els1dev = nevents_di_els/(nevents_single_els + 2 * nevents_di_els)/(nevents_single_els + 2 * nevents_di_els);
  els2dev = (0-nevents_single_els)/(nevents_single_els + 2 * nevents_di_els)/(nevents_single_els + 2 * nevents_di_els);

  corrfactor_di_els_err = std::sqrt(els1dev*els1dev*ROOT::Math::gamma_quantile_c( alpha/2 , nevents_single_els + 1, 1 ) + els2dev*els2dev*ROOT::Math::gamma_quantile_c( alpha
/2 , nevents_di_els + 1 , 1 ));
}

void AccRecoIsoEffs::NormalizeFlowNumber()
{
  double alpha;
  alpha = 1-0.6827;

  nevents_exp_all_mus_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_all_mus + 1, 1 ) );
  nevents_exp_acc_mus_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_acc_mus + 1, 1 ) );
  nevents_exp_id_mus_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_id_mus + 1, 1 ) );
  nevents_exp_iso_mus_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_iso_mus + 1, 1 ) );

  nevents_exp_all_els_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_all_els + 1, 1 ) );
  nevents_exp_acc_els_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_acc_els + 1, 1 ) );
  nevents_exp_id_els_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_id_els + 1, 1 ) );
  nevents_exp_iso_els_err = std::sqrt( ROOT::Math::gamma_quantile_c( alpha/2 , nevents_exp_iso_els + 1, 1 ) );

  nevents_pred_all_mus_err = std::sqrt(nevents_pred_all_mus_err);
  nevents_pred_acc_mus_err = std::sqrt(nevents_pred_acc_mus_err);
  nevents_pred_id_mus_err = std::sqrt(nevents_pred_id_mus_err);
  nevents_pred_iso_mus_err = std::sqrt(nevents_pred_iso_mus_err);

  nevents_pred_all_els_err = std::sqrt(nevents_pred_all_els_err);
  nevents_pred_acc_els_err = std::sqrt(nevents_pred_acc_els_err);
  nevents_pred_id_els_err = std::sqrt(nevents_pred_id_els_err);
  nevents_pred_iso_els_err = std::sqrt(nevents_pred_iso_els_err);


  double Nevents = 25446993;
  double XSec = 806.1;
  double Lumi = 1000;
  double scale = XSec*Lumi/Nevents;

  nevents_exp_all_mus *= scale;
  nevents_exp_acc_mus *= scale;
  nevents_exp_id_mus *= scale;
  nevents_exp_iso_mus *= scale;

  nevents_exp_all_els *= scale;
  nevents_exp_acc_els *= scale;
  nevents_exp_id_els *= scale;
  nevents_exp_iso_els *= scale;

  nevents_pred_all_mus *= scale;
  nevents_pred_acc_mus *= scale;
  nevents_pred_id_mus *= scale;
  nevents_pred_iso_mus *= scale;

  nevents_pred_all_els *= scale;
  nevents_pred_acc_els *= scale;
  nevents_pred_id_els *= scale;
  nevents_pred_iso_els *= scale;

  nevents_exp_all_mus_err *= scale;
  nevents_exp_acc_mus_err *= scale;
  nevents_exp_id_mus_err *= scale;
  nevents_exp_iso_mus_err *= scale;

  nevents_exp_all_els_err *= scale;
  nevents_exp_acc_els_err *= scale;
  nevents_exp_id_els_err *= scale;
  nevents_exp_iso_els_err *= scale;

  nevents_pred_all_mus_err *= scale;
  nevents_pred_acc_mus_err *= scale;
  nevents_pred_id_mus_err *= scale;
  nevents_pred_iso_mus_err *= scale;

  nevents_pred_all_els_err *= scale;
  nevents_pred_acc_els_err *= scale;
  nevents_pred_id_els_err *= scale;
  nevents_pred_iso_els_err *= scale;
}

void AccRecoIsoEffs::printNormalizeFlowNumber()
{
  std::cout<<"Normalized Flow Number information: "<<std::endl;
  
  std::cout<<"mus,exp,all: "<<nevents_exp_all_mus<<"("<<nevents_exp_all_mus_err<<")"<<std::endl;
  std::cout<<"mus,exp,acc: "<<nevents_exp_acc_mus<<"("<<nevents_exp_acc_mus_err<<")"<<std::endl;
  std::cout<<"mus,exp,id: "<<nevents_exp_id_mus<<"("<<nevents_exp_id_mus_err<<")"<<std::endl;
  std::cout<<"mus,exp,iso: "<<nevents_exp_iso_mus<<"("<<nevents_exp_iso_mus_err<<")"<<std::endl;

  std::cout<<"els,exp,all: "<<nevents_exp_all_els<<"("<<nevents_exp_all_els_err<<")"<<std::endl;
  std::cout<<"els,exp,acc: "<<nevents_exp_acc_els<<"("<<nevents_exp_acc_els_err<<")"<<std::endl;
  std::cout<<"els,exp,id: "<<nevents_exp_id_els<<"("<<nevents_exp_id_els_err<<")"<<std::endl;
  std::cout<<"els,exp,iso: "<<nevents_exp_iso_els<<"("<<nevents_exp_iso_els_err<<")"<<std::endl;

  std::cout<<"mus,pred,all: "<<nevents_pred_all_mus<<"("<<nevents_pred_all_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,acc: "<<nevents_pred_acc_mus<<"("<<nevents_pred_acc_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,id: "<<nevents_pred_id_mus<<"("<<nevents_pred_id_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,iso: "<<nevents_pred_iso_mus<<"("<<nevents_pred_iso_mus_err<<")"<<std::endl;

  std::cout<<"els,pred,all: "<<nevents_pred_all_els<<"("<<nevents_pred_all_els_err<<")"<<std::endl;
  std::cout<<"els,pred,acc: "<<nevents_pred_acc_els<<"("<<nevents_pred_acc_els_err<<")"<<std::endl;
  std::cout<<"els,pred,id: "<<nevents_pred_id_els<<"("<<nevents_pred_id_els_err<<")"<<std::endl;
  std::cout<<"els,pred,iso: "<<nevents_pred_iso_els<<"("<<nevents_pred_iso_els_err<<")"<<std::endl;

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

  std::cout<<"correction factor from di muons: "<< corrfactor_di_mus << "("<< corrfactor_di_mus_err <<")"<<std::endl;
 

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

  std::cout<<"correction factor from di electrons: "<< corrfactor_di_els <<"("<< corrfactor_di_els_err <<")"<<std::endl;

  return ;
}

