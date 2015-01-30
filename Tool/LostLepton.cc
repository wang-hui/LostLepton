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

#include "NTupleReader.h"
#include "LostLepton.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TLorentzVector.h"

//here we define the number of Pt bins for reco and iso efficiencies
#define PT_BINS 8

using namespace std;

int main()
{
  char nBase[] = "/afs/cern.ch/user/h/hua/stop/AllHadronicSUSY/CMSSW_7_2_0/stop_ttbar_skimmed_tree.root";

  TChain *fChain = new TChain("AUX");

  size_t t0 = clock();

  char fname[512];
 
  sprintf(fname, nBase, 1);
  fChain->Add(fname);

  NTupleReader tr(fChain);

  //here we define the overall information for ttbar sample
  int nevents_tot = 0;
  int nevents_sel_base = 0;
  int nevents_sel_mus = 0;
  int nevents_sel_els = 0;

  //here we defin the muon/electron number we need to count in the loop
  double nmus,nmus_acc,nels,nels_acc;
  nmus = 0;
  nmus_acc = 0;
  nels = 0;
  nels_acc = 0;

  double nmus_iso[PT_BINS] = {0}, nmus_reco[PT_BINS] = {0}, nels_iso[PT_BINS] = {0}, nels_reco[PT_BINS] = {0};
  double nmus_acc_bin[PT_BINS] = {0}, nels_acc_bin[PT_BINS] = {0};

  //here we define the baseline cut variables
  vector<TLorentzVector> jetsLVec;
  double met = 0;  
  double MT2 = 0;
  unsigned int remainPassCSVS = false;
  double linearCombmTbJetPlusmTbestTopJet = 0;
  double bestTopJetMass = 0;

  double metphi;
  double jet1_met_phi_diff, jet2_met_phi_diff, jet3_met_phi_diff;

  while(tr.getNextEvent())
  {    
    if(tr.getEvtNum()%10000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    
    nevents_tot++;

    //get baseline cut variables from root file
    jetsLVec = tr.getVec<TLorentzVector>("jetsLVec");
    met = tr.getVar<double>("met");
    MT2 = tr.getVar<double>("MT2");
    remainPassCSVS = tr.getVar<unsigned int>("remainPassCSVS");
    linearCombmTbJetPlusmTbestTopJet = tr.getVar<double>("linearCombmTbJetPlusmTbestTopJet");
    bestTopJetMass = tr.getVar<double>("bestTopJetMass");

    //to avoid the segmentation error, we need to we have at least 4 jets in the final state
    if( jetsLVec.size() < 5)
      continue;

    //std::cout << jetsLVec.at(3).Pt() << "\t" << jetsLVec.at(4).Pt() << std::endl;

    metphi = tr.getVar<double>("metphi");

    jet1_met_phi_diff = std::abs(DeltaPhi((jetsLVec.at(0)).Phi(), metphi));
    jet2_met_phi_diff = std::abs(DeltaPhi((jetsLVec.at(1)).Phi(), metphi));
    jet3_met_phi_diff = std::abs(DeltaPhi((jetsLVec.at(2)).Phi(), metphi));

    //base line cut here
    if(
       (jetsLVec.at(3)).Pt() > 50 &&                                                 
       (jetsLVec.at(4)).Pt() > 30 &&                                                 
       met>=200.0 &&                                                                  
       MT2>=300.0 &&                                                                  
       remainPassCSVS &&                                                              
       linearCombmTbJetPlusmTbestTopJet>=500.0 &&                                     
       bestTopJetMass>=80.0 &&                                                        
       bestTopJetMass<=270.0 &&
       jet1_met_phi_diff>=0.5 &&
       jet2_met_phi_diff>=0.5 &&
       jet3_met_phi_diff>=0.3
      )
    {
      nevents_sel_base++;

      int nElectrons = tr.getVar<int>("nElectrons");
      int nMuons = tr.getVar<int>("nMuons");

      vector<int> W_emuVec = tr.getVec<int>("W_emuVec");
      vector<int> W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
      vector<int> emuVec_merge;
      emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
      emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
      emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
      int gen_emus_count = emuVec_merge.size();

      if(nElectrons == 0)
      {
        nevents_sel_mus++;

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
            nmus++;
            
            double gen_mus_eta, gen_mus_phi, gen_mus_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");

            gen_mus_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_mus_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_mus_pt  = ( genDecayLVec.at ( genId ) ).Pt();

            if((std::abs(gen_mus_eta)) < 2.4 && gen_mus_pt > 5)
            {
              nmus_acc++;
              int ptbin_number;

              if(gen_mus_pt < 10)
              {
                ptbin_number = 0;
              }
              else if(gen_mus_pt >= 10 && gen_mus_pt < 20)
              {
                ptbin_number = 1;
              }
              else if(gen_mus_pt >= 20 && gen_mus_pt < 30)
              {
                ptbin_number = 2;
              }
              else if(gen_mus_pt >= 30 && gen_mus_pt < 40)
              {
                ptbin_number = 3;
              }
              else if(gen_mus_pt >= 40 && gen_mus_pt < 50)
              {
                ptbin_number = 4;
              }
              else if(gen_mus_pt >= 50 && gen_mus_pt < 70)
              {
                ptbin_number = 5;
              }
              else if(gen_mus_pt >= 70 && gen_mus_pt < 100)
              {
                ptbin_number = 6;
              }
              else if(gen_mus_pt >= 100 )
              {
                ptbin_number = 7;
              }

              nmus_acc_bin[ptbin_number]++;

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
                nmus_reco[ptbin_number]++;

                vector<double> muonsRelIso = tr.getVec<double>("muonsRelIso");
                bool mus_pass_iso;
                mus_pass_iso = false;               
                mus_pass_iso = ( muonsRelIso.at(mindeltar_index) < 0.2 );
                
                if(mus_pass_iso)
                {
                  nmus_iso[ptbin_number]++;
                }//if isolated
              }//if reconstructed
            }//if accepted
          }//if the gen particle is muon
        }//loop gen electrons/muons
      }//if no electrons

      if(nMuons == 0)
      {
        nevents_sel_els++;

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
            nels++;
          
            double gen_els_eta, gen_els_phi, gen_els_pt;
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          
            gen_els_eta = ( genDecayLVec.at ( genId ) ).Eta();
            gen_els_phi = ( genDecayLVec.at ( genId ) ).Phi();
            gen_els_pt  = ( genDecayLVec.at ( genId ) ).Pt();
    
            if((std::abs(gen_els_eta)) < 2.5 && gen_els_pt > 5)
            {
              nels_acc++;
              int ptbin_number;

              if(gen_els_pt < 10)
              {
                ptbin_number = 0;
              }
              else if(gen_els_pt >= 10 && gen_els_pt < 20)
              {
                ptbin_number = 1;
              }
              else if(gen_els_pt >= 20 && gen_els_pt < 30)
              {
                ptbin_number = 2;
              }
              else if(gen_els_pt >= 30 && gen_els_pt < 40)
              {
                ptbin_number = 3;
              }
              else if(gen_els_pt >= 40 && gen_els_pt < 50)
              {
                ptbin_number = 4;
              }
              else if(gen_els_pt >= 50 && gen_els_pt < 70)
              {
                ptbin_number = 5;
              }
              else if(gen_els_pt >= 70 && gen_els_pt < 100)
              {
                ptbin_number = 6;
              }
              else if(gen_els_pt >= 100 )
              {
                ptbin_number = 7;
              }

              nels_acc_bin[ptbin_number]++;
       
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
                nels_reco[ptbin_number]++;

                vector<double> elesRelIso = tr.getVec<double>("elesRelIso");
                bool els_pass_iso;
                els_pass_iso = false;
                els_pass_iso = ( elesRelIso.at(mindeltar_index) < 0.24 );

                if(els_pass_iso)
                {
                  nels_iso[ptbin_number]++;
                }//if isolated
              }//if reconstructed
            }//if accepted
          }//if the gen particle is electron 
        }//loop gen electrons/muons
      }//if no muons
    }
    
    //clear the vector at the end of loop
    jetsLVec.clear();

    //<< "\t" << tr.getVar<double>("joe") << "\t" << tr.getVar<int>("five") << "\t" << tr.getVec<double>("muonsMtw").size() << "\t" << tr.getVec<double>("threeNum")[2] << std::endl;
  }

  std::cout << "Lost Lepton Sample(ttbar sample) Overview:" << std::endl;
  std::cout << "nevents_tot = " << nevents_tot << std::endl;
  std::cout << "nevents_sel_base: = " << nevents_sel_base << std::endl;

  //const double ttbarCrossSection=818.8;
  //const double lumi=1000.0;
  //const double ntoteventsttbar=17184825.0;
  //const double ntoteventsttbar=3029767.0;
  
  std::cout << "nevents_sel_mus = " << nevents_sel_mus << std::endl;
  std::cout << "nevents_sel_els = " << nevents_sel_els << std::endl;

  //define acceptance, reco eff and iso eff to be calculated
  double mus_acc, els_acc;
  double mus_acc_err, els_acc_err;
  double mus_recoeff[PT_BINS]={0}, els_recoeff[PT_BINS]={0};
  double mus_isoeff[PT_BINS]={0}, els_isoeff[PT_BINS]={0};
  double mus_recoeff_err[PT_BINS]={0}, els_recoeff_err[PT_BINS]={0};
  double mus_isoeff_err[PT_BINS]={0}, els_isoeff_err[PT_BINS]={0};
  
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

  return 0;
}

