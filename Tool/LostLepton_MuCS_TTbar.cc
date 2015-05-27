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

#include "Baseline.h"
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
#include "Activity.h"

using namespace std;

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

  //TChain *fChain = new TChain("stopTreeMaker/AUX");
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
  //define activity variables
  Activity myActivity;
  //define my AccRecoIsoEffs class to stroe counts and efficiencies
  AccRecoIsoEffs myAccRecoIsoEffs;
  //define my histgram class
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram(outFileName);

  int nevents_muonCS= 0;
  int nevents_baseline= 0;
  //int nevents_baseline_ref= 0;

  //double mtwcorrfactor[8];
  //mtwcorrfactor[0] = 1.01;
  //mtwcorrfactor[1] = 1.04;
  //mtwcorrfactor[2] = 1.07;
  //mtwcorrfactor[3] = 1.20;
  //mtwcorrfactor[4] = 1.15;
  //mtwcorrfactor[5] = 1.20;
  //mtwcorrfactor[6] = 1.28;
  //mtwcorrfactor[7] = 1.83;

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

      double met = tr.getVar<double>("met");
      double metphi = tr.getVar<double>("metphi");
      int njets30 = tr.getVar<int>("cntNJetsPt30Eta24");
      const double ht = tr.getVar<double>("ht");

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

        vector<TLorentzVector> jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
        vector<double> recoJetschargedHadronEnergyFraction = tr.getVec<double>("recoJetschargedHadronEnergyFraction");
        vector<double> recoJetschargedEmEnergyFraction = tr.getVec<double>("recoJetschargedEmEnergyFraction");

        //for( unsigned int i = 0 ; i < jetsLVec.size() ; i++ )
        //{
          //(myBaseHistgram.h_b_jet_pt)->Fill( ( jetsLVec.at(i) ).Pt() );
        //}

        for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
        {
          //determine if this gen particle is Muon;
          vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          bool isMuon;
          isMuon = false;
          isMuon = ( ( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == 13 ) || ( genDecayPdgIdVec.at ( emuVec_merge.at ( gen_emus_i ) ) == -13 ) );

          if( isMuon )
          {
            ngenmu++;

            int njetsbin_number = Set_njetsbin_number(njets30);
            myAccRecoIsoEffs.nmus[njetsbin_number]++;
            
            double gen_mus_eta, gen_mus_phi, gen_mus_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

            gen_mus_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_mus_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_mus_pt  = ( genDecayLVec.at ( genId ) ).Pt();

            myActivity.getVariables(
                                    gen_mus_eta,
                                    gen_mus_phi,
                                    jetsLVec,
                                    recoJetschargedHadronEnergyFraction,
                                    recoJetschargedEmEnergyFraction
                                   );
            double activity = myActivity.getMuActivity();
            myActivity.reset();
            (myBaseHistgram.h_b_activity_mus)->Fill(activity);
            (myBaseHistgram.h_b_njets_mus)->Fill(njets30);
            //std::cout << activity  << std::endl;
            (myBaseHistgram.h_b_njets30_pt_mus)->Fill( njets30 , gen_mus_pt );
            (myBaseHistgram.h_b_njets30_eta_mus)->Fill( njets30 , gen_mus_eta );

            if( njets30 == 4 )
            {
              (myBaseHistgram.h_b_njets30_4_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_4_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_4_ht_mus)->Fill(ht);
            }
            else if( njets30 == 5 )
            {
              (myBaseHistgram.h_b_njets30_5_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_5_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_5_ht_mus)->Fill(ht);
            }
            else if( njets30 == 6 )
            {
              (myBaseHistgram.h_b_njets30_6_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_6_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_6_ht_mus)->Fill(ht);
            }
            else if( njets30 == 7 )
            {
              (myBaseHistgram.h_b_njets30_7_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_7_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_7_ht_mus)->Fill(ht);
            }
            else if( njets30 == 8 )
            {
              (myBaseHistgram.h_b_njets30_8_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_8_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_8_ht_mus)->Fill(ht);
            }
            else if( njets30 == 9 )
            {
              (myBaseHistgram.h_b_njets30_9_pt_mus)->Fill( gen_mus_pt );
              (myBaseHistgram.h_b_njets30_9_eta_mus)->Fill( gen_mus_eta );
	      (myBaseHistgram.h_b_njets30_9_ht_mus)->Fill(ht);
            }

            if((std::abs(gen_mus_eta)) < AnaConsts::muonsArr[1] && gen_mus_pt > AnaConsts::muonsArr[2])
            {
              myAccRecoIsoEffs.nmus_acc[njetsbin_number]++;

              int ptbin_number = Set_ptbin_number(gen_mus_pt);
              int acbin_number = Set_acbin_number(activity);

              myAccRecoIsoEffs.nmus_acc_bin[ptbin_number][acbin_number]++;

              for( int i = 0 ; i < genDecayPdgIdVec.size() ; i++ )
              {
                int genindice = genDecayPdgIdVec.at(i);
                if( ( genindice == 13 || genindice  == -13 ) && i != genId )
                {
                  double deltar_study;
                  deltar_study = DeltaR(
                                        gen_mus_eta,
                                        gen_mus_phi,
                                        ( genDecayLVec.at ( i ) ).Eta(),
                                        ( genDecayLVec.at ( i ) ).Phi()
                                       );
                  (myBaseHistgram.h_b_deltaR_genup_mus)->Fill(deltar_study);
                }
                else
                  continue;
              }

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
              (myBaseHistgram.h_b_deltaR_mus)->Fill( deltar );
              (myBaseHistgram.h_b_deltaR_pt_mus)->Fill( deltar , gen_mus_pt );

              deltar_mus_pool.clear();

              bool ismatcheddeltaR;
              ismatcheddeltaR = (deltar < 0.2);

              if(ismatcheddeltaR
                 //&& 
                 //isgoodmuonid
                )
              {
                myAccRecoIsoEffs.nmus_reco[ptbin_number][acbin_number]++;
                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                double reco_mus_pt = (muonsLVec.at(mindeltar_index)).Pt();
                int ptbin_number_allreco = Set_ptbin_number(reco_mus_pt);

                double reco_mus_eta = (muonsLVec.at(mindeltar_index)).Eta();
                double reco_mus_phi = (muonsLVec.at(mindeltar_index)).Phi();
                myActivity.getVariables(
                                        reco_mus_eta,
                                        reco_mus_phi,
                                        jetsLVec,
                                        recoJetschargedHadronEnergyFraction,
                                        recoJetschargedEmEnergyFraction
                                       );
                double activity_allreco = myActivity.getMuActivity();
                myActivity.reset();
                int acbin_number_allreco = Set_acbin_number(activity_allreco);

                myAccRecoIsoEffs.nmus_reco_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //vector<double> muonsRelIso = tr.getVec<double>("muonsRelIso");
                vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");

                bool mus_pass_iso;
                mus_pass_iso = false;               
                mus_pass_iso = ( muonsMiniIso.at(mindeltar_index) < AnaConsts::muonsArr[4] );
                
                if(mus_pass_iso)
                {
                  myAccRecoIsoEffs.nmus_iso[ptbin_number][acbin_number]++;
                  myAccRecoIsoEffs.nmus_iso_allreco[ptbin_number_allreco][acbin_number_allreco]++;
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
        
        vector<TLorentzVector> jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
        vector<double> recoJetschargedHadronEnergyFraction = tr.getVec<double>("recoJetschargedHadronEnergyFraction");
        vector<double> recoJetschargedEmEnergyFraction = tr.getVec<double>("recoJetschargedEmEnergyFraction");

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

            int njetsbin_number = Set_njetsbin_number(njets30);
            myAccRecoIsoEffs.nels[njetsbin_number]++;
          
            double gen_els_eta, gen_els_phi, gen_els_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          
            gen_els_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_els_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_els_pt  = ( genDecayLVec.at ( genId ) ).Pt();
    
            myActivity.getVariables(
                                    gen_els_eta,
                                    gen_els_phi,
                                    jetsLVec,
                                    recoJetschargedHadronEnergyFraction,
                                    recoJetschargedEmEnergyFraction
                                   );
            double activity = myActivity.getElActivity();
            myActivity.reset();
            (myBaseHistgram.h_b_activity_els)->Fill(activity);
            (myBaseHistgram.h_b_njets_els)->Fill(njets30);

            (myBaseHistgram.h_b_njets30_pt_els)->Fill( njets30 , gen_els_pt );
            (myBaseHistgram.h_b_njets30_eta_els)->Fill( njets30 , gen_els_eta );

            if( njets30 == 4 )
            {
              (myBaseHistgram.h_b_njets30_4_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_4_eta_els)->Fill( gen_els_eta );
            }
            else if( njets30 == 5 )
            {
              (myBaseHistgram.h_b_njets30_5_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_5_eta_els)->Fill( gen_els_eta );
            }
            else if( njets30 == 6 )
            {
              (myBaseHistgram.h_b_njets30_6_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_6_eta_els)->Fill( gen_els_eta );
            }
            else if( njets30 == 7 )
            {
              (myBaseHistgram.h_b_njets30_7_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_7_eta_els)->Fill( gen_els_eta );
            }
            else if( njets30 == 8 )
            {
              (myBaseHistgram.h_b_njets30_8_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_8_eta_els)->Fill( gen_els_eta );
            }
            else if( njets30 == 9 )
            {
              (myBaseHistgram.h_b_njets30_9_pt_els)->Fill( gen_els_pt );
              (myBaseHistgram.h_b_njets30_9_eta_els)->Fill( gen_els_eta );
            }

            if( (std::abs(gen_els_eta)) < AnaConsts::elesArr[1] && gen_els_pt > AnaConsts::elesArr[2] )
            {
              myAccRecoIsoEffs.nels_acc[njetsbin_number]++;

              int ptbin_number = Set_ptbin_number(gen_els_pt);
              int acbin_number = Set_acbin_number(activity);

              myAccRecoIsoEffs.nels_acc_bin[ptbin_number][acbin_number]++;
       
              for( int i = 0 ; i < genDecayPdgIdVec.size() ; i++ )
              {
                int genindice = genDecayPdgIdVec.at(i);
                if( ( genindice == 11 || genindice  == -11 ) && i != genId )
                {
                  double deltar_study;
                  deltar_study = DeltaR(
                                        gen_els_eta,
                                        gen_els_phi,
                                        ( genDecayLVec.at ( i ) ).Eta(),
                                        ( genDecayLVec.at ( i ) ).Phi()
                                       );
                  (myBaseHistgram.h_b_deltaR_genup_els)->Fill(deltar_study);
                }
                else
                  continue;
              }

              //loop over reco lepton information to determine the smallest deltar
              vector<double> deltar_els_pool;
              for(int reco_els_i = 0 ; reco_els_i < reco_els_count ; reco_els_i++)
              {
                double deltar_media;
                deltar_media = DeltaR(
                                      gen_els_eta,
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
              (myBaseHistgram.h_b_deltaR_els)->Fill(deltar);
              (myBaseHistgram.h_b_deltaR_pt_els)->Fill( deltar , gen_els_pt );

              deltar_els_pool.clear();

              bool ismatcheddeltaR;
              ismatcheddeltaR = (deltar < 0.2);
    
              if(ismatcheddeltaR
                 //&& 
                 //isgoodmuonid
                )
              {
                myAccRecoIsoEffs.nels_reco[ptbin_number][acbin_number]++;

                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                double reco_els_pt = (elesLVec.at(mindeltar_index)).Pt();
                int ptbin_number_allreco = Set_ptbin_number(reco_els_pt);

                double reco_els_eta = (elesLVec.at(mindeltar_index)).Eta();
                double reco_els_phi = (elesLVec.at(mindeltar_index)).Phi();
                myActivity.getVariables(
                                        reco_els_eta,
                                        reco_els_phi,
                                        jetsLVec,
                                        recoJetschargedHadronEnergyFraction,
                                        recoJetschargedEmEnergyFraction
                                       );
                double activity_allreco = myActivity.getElActivity();
                myActivity.reset();
                int acbin_number_allreco = Set_acbin_number(activity_allreco);

                myAccRecoIsoEffs.nels_reco_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //vector<double> elesRelIso = tr.getVec<double>("elesRelIso");
                vector<double> elesMiniIso = tr.getVec<double>("elesMiniIso");

                bool els_pass_iso;
                els_pass_iso = false;
                els_pass_iso = ( elesMiniIso.at(mindeltar_index) < AnaConsts::elesArr[4] );

                if(els_pass_iso)
                {
                  myAccRecoIsoEffs.nels_iso[ptbin_number][acbin_number]++;
                  myAccRecoIsoEffs.nels_iso_allreco[ptbin_number_allreco][acbin_number_allreco]++;
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
      
      //mtW correction factor calculation
      if (nElectrons == 0 && nMuons == 1)
      {
        vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
        vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");

        double reco_mus_pt = 0, reco_mus_eta = 0, reco_mus_phi = 0;

        for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
        {
          if( fabs(muonsLVec[im].Eta()) < 2.4 && muonsMiniIso[im] < 0.2 )
          {
            reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
            reco_mus_eta = ( muonsLVec.at(im) ).Eta();
            reco_mus_phi = ( muonsLVec.at(im) ).Phi();
          }
        }

        double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
        double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * ( 1.0 - cos(deltaphi_mus) ) );
        (myBaseHistgram.h_mtw_mus)->Fill(mtW_mus);

        int ptbin_number_allreco = Set_ptbin_number(reco_mus_pt);

        myAccRecoIsoEffs.mtwall[ptbin_number_allreco]++;        
        if( mtW_mus < 100 )
        {
          myAccRecoIsoEffs.mtw100[ptbin_number_allreco]++;
        }
      }

      //book the variable we need to check
      int ntopjets = tr.getVar<int>("nTopCandSortedCnt");
      double MT2 = tr.getVar<double>("MT22");
      double bestTopJetMass = tr.getVar<double>("bestTopJetMass2");
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
  myAccRecoIsoEffs.NumberstoEffs();
  myAccRecoIsoEffs.EffsPlotsGen();
  myAccRecoIsoEffs.EffstoWeights();
  myAccRecoIsoEffs.GetDiLeptonFactor();
  myAccRecoIsoEffs.printAccRecoIsoEffs();
  myAccRecoIsoEffs.printEffsHeader();

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
      int ntopjets = trCS.getVar<int>("nTopCandSortedCnt");
      double MT2 = trCS.getVar<double>("MT22");
      double bestTopJetMass = trCS.getVar<double>("bestTopJetMass2");
      double ht = trCS.getVar<double>("ht");
      double mht = trCS.getVar<double>("mht");

      //muon CS
      //nMuons in flatree means no iso cut muons; CUT2 we add iso
      if (nElectrons == 0 && nMuons == 1)
      {
        //counting the events for muon control sample
        myAccRecoIsoEffs.nevents_cs_mus++;
        //get muon variables
	vector<TLorentzVector> muonsLVec = trCS.getVec<TLorentzVector>("muonsLVec");
        //vector<double> muonsRelIso = trCS.getVec<double>("muonsRelIso");
        vector<double> muonsMiniIso = trCS.getVec<double>("muonsMiniIso");
        //get jet variables for activity calculation
        vector<TLorentzVector> jetsLVec = trCS.getVec<TLorentzVector>("jetsLVec");
        vector<double> recoJetschargedHadronEnergyFraction = trCS.getVec<double>("recoJetschargedHadronEnergyFraction");
        vector<double> recoJetschargedEmEnergyFraction = trCS.getVec<double>("recoJetschargedEmEnergyFraction");

        double reco_mus_pt = 0, reco_mus_eta = 0, reco_mus_phi = 0;

        for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
        {
          if( fabs(muonsLVec[im].Eta()) < 2.4 && muonsMiniIso[im] < 0.2 )
	  {
            reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
            reco_mus_eta = ( muonsLVec.at(im) ).Eta();
            reco_mus_phi = ( muonsLVec.at(im) ).Phi();
	  }
	}

        double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
        double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * ( 1.0 - cos(deltaphi_mus) ) );

        myActivity.getVariables(
                                reco_mus_eta,
                                reco_mus_phi,
                                jetsLVec,
                                recoJetschargedHadronEnergyFraction,
                                recoJetschargedEmEnergyFraction
                               );
        double activity = myActivity.getMuActivity();
        myActivity.reset();

        if ( mtW_mus < 100.0 )
        {
	  //////////////////////////
	  // prediction computation
	  //////////////////////////
	  double EventWeight_mus = 1.0;
          int njetsbin_number = Set_njetsbin_number(njets30);
          int ptbin_number = Set_ptbin_number(reco_mus_pt);
          int acbin_number = Set_acbin_number(activity);

	  //mtwcorrfactor
	  EventWeight_mus = EventWeight_mus * myAccRecoIsoEffs.mtwcorrfactor[ptbin_number];
	  //dimuon correction factor
          EventWeight_mus = EventWeight_mus * myAccRecoIsoEffs.corrfactor_di_mus;

          //muon prediction from muon CS
	  //Fill muon iso closure plots
	  (myBaseHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_iso_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  //Fill muon id closure plots
	  (myBaseHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_id_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  //Fill muon acc closure plots
	  (myBaseHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_acc_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	  //Fill all muon closure plots
	  double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number];

	  (myBaseHistgram.h_pred_mu_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
	  (myBaseHistgram.h_pred_mu_all_topmass)->Fill(bestTopJetMass, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
          (myBaseHistgram.h_pred_mu_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

          //total events flow for muons, prediction
          myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
        
          myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
          myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;

          //begin to predict lost electrons from muon CS
          double EventWeight_els = 1.0;
          //mtwcorrfactor
          EventWeight_els = EventWeight_els * myAccRecoIsoEffs.mtwcorrfactor[ptbin_number];
          //dielectron correction factor
          EventWeight_els = EventWeight_els * myAccRecoIsoEffs.corrfactor_di_els;
  
          //electron prediction from muon CS
          //Fill electron iso closure plots
          (myBaseHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          //Fill electron id closure plots
          (myBaseHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          //Fill electron acc closure plots
          (myBaseHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_topmass)->Fill(bestTopJetMass, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          (myBaseHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
          //Fill all electron closure plots
          double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number];

          (myBaseHistgram.h_pred_el_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_topmass)->Fill(bestTopJetMass, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
          (myBaseHistgram.h_pred_el_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);
          //total events flow for electrons, prediction
          myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;

          myAccRecoIsoEffs.nevents_pred_all_els_err += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
          myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
        }//mtW5_els<100 (muon CS)
      }//nElectrons == 0 && nMuons == 1 (muon CS)
    }//baseline_nolepveto
  }

  myAccRecoIsoEffs.printOverview();
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

  std::cout << "nevents_cs_mus = "<< nevents_cs_mus << std::endl;

  return ;
}

void AccRecoIsoEffs::NumberstoEffs()
{
  //mus_acc_err = std::sqrt( get_stat_Error(nmus_acc,nmus)*get_stat_Error(nmus_acc,nmus) + get_sys_Error(mus_acc,0.09)*get_sys_Error(mus_acc,0.09) );
  //els_acc_err = std::sqrt( get_stat_Error(nels_acc,nels)*get_stat_Error(nels_acc,nels) + get_sys_Error(els_acc,0.09)*get_sys_Error(els_acc,0.09) ); 

  int i_cal;
  int j_cal;

  for(i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++)
  {
    mus_acc[i_cal] = nmus_acc[i_cal]/nmus[i_cal];
    mus_acc_err[i_cal] = std::sqrt( get_stat_Error(nmus_acc[i_cal],nmus[i_cal])*get_stat_Error(nmus_acc[i_cal],nmus[i_cal]) + get_sys_Error(mus_acc[i_cal],0.09)*get_sys_Error(mus_acc[i_cal],0.09) );

    els_acc[i_cal] = nels_acc[i_cal]/nels[i_cal];
    els_acc_err[i_cal] = std::sqrt( get_stat_Error(nels_acc[i_cal],nels[i_cal])*get_stat_Error(nels_acc[i_cal],nels[i_cal]) + get_sys_Error(els_acc[i_cal],0.09)*get_sys_Error(els_acc[i_cal],0.09) );
  }

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      mus_recoeff[i_cal][j_cal]     = nmus_reco[i_cal][j_cal]/nmus_acc_bin[i_cal][j_cal];
      mus_recoeff_err[i_cal][j_cal] = get_stat_Error(nmus_reco[i_cal][j_cal],nmus_acc_bin[i_cal][j_cal]);
      els_recoeff[i_cal][j_cal]     = nels_reco[i_cal][j_cal]/nels_acc_bin[i_cal][j_cal];
      els_recoeff_err[i_cal][j_cal] = get_stat_Error(nels_reco[i_cal][j_cal],nels_acc_bin[i_cal][j_cal]);

      mus_isoeff[i_cal][j_cal]     = nmus_iso[i_cal][j_cal]/nmus_reco[i_cal][j_cal];
      mus_isoeff_err[i_cal][j_cal] = get_stat_Error(nmus_iso[i_cal][j_cal],nmus_reco[i_cal][j_cal]);
      els_isoeff[i_cal][j_cal]     = nels_iso[i_cal][j_cal]/nels_reco[i_cal][j_cal];
      els_isoeff_err[i_cal][j_cal] = get_stat_Error(nels_iso[i_cal][j_cal],nels_reco[i_cal][j_cal]);

      mus_isoeff_allreco[i_cal][j_cal]     = nmus_iso_allreco[i_cal][j_cal]/nmus_reco_allreco[i_cal][j_cal];
      mus_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nmus_iso_allreco[i_cal][j_cal],nmus_reco_allreco[i_cal][j_cal]);
      els_isoeff_allreco[i_cal][j_cal]     = nels_iso_allreco[i_cal][j_cal]/nels_reco_allreco[i_cal][j_cal];
      els_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nels_iso_allreco[i_cal][j_cal],nels_reco_allreco[i_cal][j_cal]);
    }
  }
  
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    mtwcorrfactor[i_cal] = mtwall[i_cal]/mtw100[i_cal];
    mtwcorrfactor_err[i_cal] = get_stat_Error(mtw100[i_cal],mtwall[i_cal]);
  }

  return ;
}

void AccRecoIsoEffs::EffsPlotsGen()
{
  int i_cal;
  int j_cal;

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      mus_recoeffs2d->SetBinContent( i_cal+1 , j_cal+1, mus_recoeff[i_cal][j_cal] );
      mus_isoeffs2d->SetBinContent( i_cal+1 , j_cal+1, mus_isoeff_allreco[i_cal][j_cal] );
      els_recoeffs2d->SetBinContent( i_cal+1 , j_cal+1, els_recoeff[i_cal][j_cal] );
      els_isoeffs2d->SetBinContent( i_cal+1 , j_cal+1, els_isoeff_allreco[i_cal][j_cal] );

      mus_recoeffs2d->SetBinError( i_cal+1 , j_cal+1, mus_recoeff_err[i_cal][j_cal] );
      mus_isoeffs2d->SetBinError( i_cal+1 , j_cal+1, mus_isoeff_err_allreco[i_cal][j_cal] );
      els_recoeffs2d->SetBinError( i_cal+1 , j_cal+1, els_recoeff_err[i_cal][j_cal] );
      els_isoeffs2d->SetBinError( i_cal+1 , j_cal+1, els_isoeff_err_allreco[i_cal][j_cal] );
    }       
  }

  mus_recoeffs2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_recoeffs2d->GetYaxis()->SetTitle("Activity");
  mus_isoeffs2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_isoeffs2d->GetYaxis()->SetTitle("Activity");
  els_recoeffs2d->GetXaxis()->SetTitle("Electron Pt [GeV]");
  els_recoeffs2d->GetYaxis()->SetTitle("Activity");
  els_isoeffs2d->GetXaxis()->SetTitle("Electron Pt [GeV]");
  els_isoeffs2d->GetYaxis()->SetTitle("Activity");

  Effs2dPlots->Write();
}

void AccRecoIsoEffs::EffstoWeights()
{
  int i_cal;
  int j_cal;
  int k_cal;

  for(i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < PT_BINS ; j_cal++)
    {
      for(k_cal = 0 ; k_cal < AC_BINS ; k_cal++)
      {
        //mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - mus_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal];
        //mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1.0 - mus_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] ); 
        //mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - mus_acc[i_cal])/mus_acc[i_cal] );
    
        //els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - els_acc[i_cal])/mus_acc[i_cal] );
        //els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1 - els_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] )* (els_acc[i_cal]/mus_acc[i_cal]); 
        //els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - els_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal] ) * (els_recoeff[j_cal][k_cal]/mus_recoeff[j_cal][k_cal])* (els_acc[i_cal]/mus_acc[i_cal]);

        mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - mus_isoeff_allreco[j_cal][k_cal])/mus_isoeff_allreco[j_cal][k_cal];
        mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * ( (1.0 - mus_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] );
        mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - mus_acc[i_cal])/mus_acc[i_cal] );

        els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - els_acc[i_cal])/mus_acc[i_cal] );
        els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * ( (1 - els_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] )* (els_acc[i_cal]/mus_acc[i_cal]);
        els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - els_isoeff_allreco[j_cal][k_cal])/mus_isoeff_allreco[j_cal][k_cal] ) * (els_recoeff[j_cal][k_cal]/mus_recoeff[j_cal][k_cal])* (els_acc[i_cal]/mus_acc[i_cal]);
      }
    }
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

  std::cout<<"mus,pred,all: "<<nevents_pred_all_mus<<"("<<nevents_pred_all_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,acc: "<<nevents_pred_acc_mus<<"("<<nevents_pred_acc_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,id: "<<nevents_pred_id_mus<<"("<<nevents_pred_id_mus_err<<")"<<std::endl;
  std::cout<<"mus,pred,iso: "<<nevents_pred_iso_mus<<"("<<nevents_pred_iso_mus_err<<")"<<std::endl;

  std::cout<<"mus,percentage,all: "<< nevents_pred_all_mus/nevents_exp_all_mus - 1 << std::endl;
  std::cout<<"mus,percentage,acc: "<< nevents_pred_acc_mus/nevents_exp_acc_mus - 1 << std::endl;
  std::cout<<"mus,percentage,id: " << nevents_pred_id_mus/nevents_exp_id_mus - 1   << std::endl;
  std::cout<<"mus,percentage,iso: "<< nevents_pred_iso_mus/nevents_exp_iso_mus - 1 << std::endl;

  std::cout<<"els,exp,all: "<<nevents_exp_all_els<<"("<<nevents_exp_all_els_err<<")"<<std::endl;
  std::cout<<"els,exp,acc: "<<nevents_exp_acc_els<<"("<<nevents_exp_acc_els_err<<")"<<std::endl;
  std::cout<<"els,exp,id: "<<nevents_exp_id_els<<"("<<nevents_exp_id_els_err<<")"<<std::endl;
  std::cout<<"els,exp,iso: "<<nevents_exp_iso_els<<"("<<nevents_exp_iso_els_err<<")"<<std::endl;

  std::cout<<"els,pred,all: "<<nevents_pred_all_els<<"("<<nevents_pred_all_els_err<<")"<<std::endl;
  std::cout<<"els,pred,acc: "<<nevents_pred_acc_els<<"("<<nevents_pred_acc_els_err<<")"<<std::endl;
  std::cout<<"els,pred,id: "<<nevents_pred_id_els<<"("<<nevents_pred_id_els_err<<")"<<std::endl;
  std::cout<<"els,pred,iso: "<<nevents_pred_iso_els<<"("<<nevents_pred_iso_els_err<<")"<<std::endl;

  std::cout<<"els,percentage,all: "<< nevents_pred_all_els/nevents_exp_all_els - 1 << std::endl;
  std::cout<<"els,percentage,acc: "<< nevents_pred_acc_els/nevents_exp_acc_els - 1 << std::endl;
  std::cout<<"els,percentage,id: " << nevents_pred_id_els/nevents_exp_id_els - 1   << std::endl;
  std::cout<<"els,percentage,iso: "<< nevents_pred_iso_els/nevents_exp_iso_els - 1 << std::endl;


}

void AccRecoIsoEffs::printAccRecoIsoEffs()
{
  int i_cal = 0;
  int j_cal = 0;
  std::cout.precision(3);

  std::cout << "mtW correction factor: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    std::cout << mtwcorrfactor[i_cal] << "(" << mtwcorrfactor_err[i_cal] << ")"<< " ";;
    if( i_cal == PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }


  std::cout << std::endl << "Muon information: " << std::endl;

  std::cout << "number of muons from top: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << nmus[i_cal] << " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, accepted: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << nmus_acc[i_cal] << " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, accepted, bins: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nmus_acc_bin[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal == PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, reconstructed: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nmus_reco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, reconstructed (allreco): " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nmus_reco_allreco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, isolated: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nmus_iso[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of muons from top, isolated (allreco): " << std::endl;

  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nmus_iso_allreco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }


  std::cout << "muons from top, acceptance: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << mus_acc[i_cal] << "(" << mus_acc_err[i_cal] << ")"<< " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  } 

  std::cout << "muons from top, reconstruction efficiency: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << mus_recoeff[i_cal][j_cal] << "(" << mus_recoeff_err[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout<<"muons from top, isolation efficiency: " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << mus_isoeff[i_cal][j_cal] << "(" << mus_isoeff_err[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout<<"muons from top, isolation efficiency (allreco): " << std::endl;
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << mus_isoeff_allreco[i_cal][j_cal] << "(" << mus_isoeff_err_allreco[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "correction factor from di muons: "<< corrfactor_di_mus << "(" << corrfactor_di_mus_err << ")" << std::endl;

  std::cout << std::endl << "Electron information: " << std::endl;

  std::cout << "number of electrons from top: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << nels[i_cal] << " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of electrons from top, accepted: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << nels_acc[i_cal] << " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  }


  std::cout << "number of electrons from top, accepted, bins: " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nels_acc_bin[i_cal][j_cal] <<" ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout<<std::endl;
    }
  }

  std::cout << "number of electrons from top, reconstructed: " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nels_reco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of electrons from top, reconstructed (allreco): " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nels_reco_allreco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }


  std::cout << "number of electrons from top, isolated: " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nels_iso[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "number of electrons from top, isolated (allreco): " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << nels_iso_allreco[i_cal][j_cal] << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "electrons from top, acceptance: " << std::endl;
  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    std::cout << els_acc[i_cal] << "(" << els_acc_err[i_cal] << ")"<< " ";
    if( i_cal == NJETS_BINS-1 )
    {
      std::cout << std::endl;
    }
  }


  std::cout << "electrons from top, reconstruction efficiency: " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << els_recoeff[i_cal][j_cal] << "(" << els_recoeff_err[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal == PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "electrons from top, isolation efficiency: " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << els_isoeff[i_cal][j_cal] << "(" << els_isoeff_err[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "electrons from top, isolation efficiency (allreco): " << std::endl;
  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      std::cout << els_isoeff_allreco[i_cal][j_cal] << "(" << els_isoeff_err_allreco[i_cal][j_cal] << ")" << " ";
      if( j_cal == AC_BINS-1 )
      {
        std::cout << std::endl;
      }
    }
    if( i_cal==PT_BINS-1 )
    {
      std::cout << std::endl;
    }
  }

  std::cout << "correction factor from di electrons: " << corrfactor_di_els << "(" << corrfactor_di_els_err << ")" <<std::endl;

  return ;
}


void AccRecoIsoEffs::printEffsHeader()
{
  ofstream EffsHeader;
  EffsHeader.open ("EffsHeader_MuCS.h");

  int i_cal = 0;
  int j_cal = 0;

  EffsHeader << "  const double ttbar_mtwcorrfactor[" << PT_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    if( i_cal == 0 ) { EffsHeader << "{"; }
    EffsHeader << mtwcorrfactor[i_cal];
    if( i_cal != PT_BINS-1 ) { EffsHeader << ","; }
    if( i_cal == PT_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }

  EffsHeader << "  const double ttbar_mus_acc[" << NJETS_BINS << "] = "; 
  for( i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    if( i_cal == 0 ) { EffsHeader << "{"; }
    EffsHeader << mus_acc[i_cal];
    if( i_cal != NJETS_BINS-1 ) { EffsHeader << ","; }
    if( i_cal == NJETS_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }

  EffsHeader << "  const double ttbar_mus_recoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < AC_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << mus_recoeff[i_cal][j_cal];            
      if( j_cal != AC_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }

  EffsHeader << "  const double ttbar_mus_isoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < AC_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << mus_isoeff_allreco[i_cal][j_cal];
      if( j_cal != AC_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }

  EffsHeader << "  const double ttbar_els_acc[" << NJETS_BINS << "] = ";                            
  for( i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    if( i_cal == 0 ) { EffsHeader << "{"; }
    EffsHeader << els_acc[i_cal];
    if( i_cal != NJETS_BINS-1 ) { EffsHeader << ","; }
    if( i_cal == NJETS_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }


  EffsHeader << "  const double ttbar_els_recoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < AC_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << els_recoeff[i_cal][j_cal];
      if( j_cal != AC_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }

  EffsHeader << "  const double ttbar_els_isoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < AC_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << els_isoeff_allreco[i_cal][j_cal];
      if( j_cal != AC_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == AC_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }


  EffsHeader << "  const double ttbar_corrfactor_di_mus = " << corrfactor_di_mus << ";" << std::endl;
  EffsHeader << "  const double ttbar_corrfactor_di_els = " << corrfactor_di_els << ";" << std::endl;

  EffsHeader.close();
}
