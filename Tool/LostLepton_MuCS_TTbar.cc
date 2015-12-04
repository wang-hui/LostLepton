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

#include "Activity.h"
#include "LostLepton_MuCS_TTbar.h"
#include "TTJetsReWeighting.h"
#include "v151201_EffsHeader_MuCS.h"

const double isotrackvetoeff = 0.563499421;
const bool applyisotrkveto = false;
//const double isotrackvetoeff = 1;


void LoopLLCal( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight )
{
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram("test.root");

  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  std::cout << "Efficiencies Calculation: " << std::endl;

  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  {  
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    AnaFunctions::prepareTopTagger();
    //The passBaseline is registered here
    tr.registerFunction(&mypassBaselineFunc);    
 
    double thisweight = (*iter_TTJetsSampleInfos).weight;
    std::cout << "Weight " << thisweight << std::endl;
    int neventc=0;
    while(tr.getNextEvent())
    //while(tr.getNextEvent() && neventc<10000)
    {
      ++neventc;
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    
      myAccRecoIsoEffs.nevents_tot+=thisweight;

      //baseline cut without lepton veto
      bool passBaselinelostlept = tr.getVar<bool>("passBaseline"+spec);

      if (passBaselinelostlept)
      {
        myAccRecoIsoEffs.nevents_sel_base+=thisweight;

        int nElectrons = tr.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = tr.getVar<int>("nMuons_CUT"+spec);


        //searchbin variables
        int ntopjets = tr.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = tr.getVar<int>("cntCSVS"+spec);
        double MT2 = tr.getVar<double>("best_had_brJet_MT2"+spec);
        double met = tr.getVar<double>("met");
        double metphi = tr.getVar<double>("metphi");
	const double ht = tr.getVar<double>("ht");

        int ngenmu = 0;
        int ngenel = 0;

        int njets30 = tr.getVar<int>("cntNJetsPt30Eta24"+spec);

        const int nIsoTrks = tr.getVar<int>("nIsoTrks_CUT"+spec);
        //std::cout << "nIsoTrks = " << nIsoTrks << std::endl;

        //merge electron, muon generation level information
        std::vector<int> W_emuVec = tr.getVec<int>("W_emuVec");
        std::vector<int> W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
        std::vector<int> emuVec_merge;
        emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
        emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
        emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
        int gen_emus_count = emuVec_merge.size();

        std::vector<double> W_emu_pfActivityVec = tr.getVec<double>("W_emu_pfActivityVec");
        std::vector<double> W_tau_emu_pfActivityVec = tr.getVec<double>("W_tau_emu_pfActivityVec");
        std::vector<double> emu_pfActivityVec_merge;
        emu_pfActivityVec_merge.reserve( W_emu_pfActivityVec.size() + W_tau_emu_pfActivityVec.size() );
        emu_pfActivityVec_merge.insert( emu_pfActivityVec_merge.end(), W_emu_pfActivityVec.begin(), W_emu_pfActivityVec.end() );
        emu_pfActivityVec_merge.insert( emu_pfActivityVec_merge.end(), W_tau_emu_pfActivityVec.begin(), W_tau_emu_pfActivityVec.end() );

	double thisweight = (*iter_TTJetsSampleInfos).weight;

        if(nElectrons == 0)
        {
          myAccRecoIsoEffs.nevents_sel_mus+=thisweight;

          //get gen level information of leptons
          std::vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          std::vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          //get reco level information of muons
          std::vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
          std::vector<double> muonspfActivity = tr.getVec<double>("muonspfActivity");
          std::vector<int> muonsFlagMedium = tr.getVec<int>("muonsFlagMedium");
          std::vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");

          for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
          {
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            double genmuonpfActivity = emu_pfActivityVec_merge.at(gen_emus_i);

            LostLeptonObj myLostMuonObj;
            
            myLostMuonObj.SetMyLL(
                                   //LL variables that will be always useful 
                                   genDecayPdgIdVec.at ( genId ),
                             	   ( genDecayLVec.at ( genId ) ),
                                   genmuonpfActivity,
                                   muonsFlagMedium,
                                   muonsLVec,
                                   muonspfActivity,
                                   muonsMiniIso
                                 );

            if( myLostMuonObj.isMu )
            {
	      ngenmu++;
              const int njetsbin_number = Set_njetsbin_number(njets30);
	      const int HTbin_number = Set_HTbin_number(ht);
	      //std::cout << "HTbin_number = " << HTbin_number << std::endl;            

              myAccRecoIsoEffs.nmus[njetsbin_number][HTbin_number]+=thisweight;
              //myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
              myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]++;
            
              if( myLostMuonObj.passAcc )
              {
                myAccRecoIsoEffs.nmus_acc[njetsbin_number][HTbin_number]+=thisweight;
                //myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
                myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]++;

                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);

                myAccRecoIsoEffs.nmus_acc_bin[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nmus_acc_bin_MC[ptbin_number][acbin_number]++;
              }

              //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
              if( myLostMuonObj.passId )
              {
                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);

                myAccRecoIsoEffs.nmus_reco[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nmus_reco_MC[ptbin_number][acbin_number]++;

                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostMuonObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostMuonObj.reco_activity);

                myAccRecoIsoEffs.nmus_reco_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight;
                myAccRecoIsoEffs.nmus_reco_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //std::cout << myLostMuonObj.gen_activity <<"," << myLostMuonObj.reco_activity << std::endl;
                (myBaseHistgram.h_id_genactivity_mus)->Fill(myLostMuonObj.gen_activity,thisweight);
                (myBaseHistgram.h_id_recoactivity_mus)->Fill(myLostMuonObj.reco_activity,thisweight);
              }

              if( myLostMuonObj.passIso )
              {
                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);

                myAccRecoIsoEffs.nmus_iso[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nmus_iso_MC[ptbin_number][acbin_number]++;
              
                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostMuonObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostMuonObj.reco_activity);
                
                myAccRecoIsoEffs.nmus_iso_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight;
                myAccRecoIsoEffs.nmus_iso_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //std::cout << muonsMiniIso.size() << "," << myLostMuonObj.reco_index << "ISO:" << muonsMiniIso.at(myLostMuonObj.reco_index) << std::endl;
              }
              //check warning function when we calculate the efficienies!
              if ( nMuons == 0 && myLostMuonObj.passIso )
              { 
		std::cout << "Warning: mu pass iso but nMuons=0!" << std::endl;
		std::cout << "reco_mus_pt = " << myLostMuonObj.reco_pt << std::endl;
		std::cout << "reco_mus_eta = " << myLostMuonObj.reco_eta << std::endl;
              }
            }//if the gen particle is muon
          }//loop gen electrons/muons

	  // dilepton computation
	  if (nMuons==0)
	  {
	    if (ngenmu==1) myAccRecoIsoEffs.nevents_single_mus+=thisweight;
	    if (ngenmu==2) myAccRecoIsoEffs.nevents_di_mus+=thisweight;
	  }


        }//if no electrons

        if(nMuons == 0)
        {
          myAccRecoIsoEffs.nevents_sel_els+=thisweight;

          //get gen level information of leptons
          std::vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          std::vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          //get reco level information of muons
          std::vector<TLorentzVector> elesLVec = tr.getVec<TLorentzVector>("elesLVec");
          std::vector<double> elespfActivity = tr.getVec<double>("elespfActivity");
          std::vector<int> elesFlagVeto = tr.getVec<int>("elesFlagVeto");
          std::vector<double> elesMiniIso = tr.getVec<double>("elesMiniIso");

          for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
          {
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            double genelectronpfActivity = emu_pfActivityVec_merge.at(gen_emus_i);

            LostLeptonObj myLostElectronObj;

            myLostElectronObj.SetMyLL(
                                       //LL variables that will be always useful 
                                       genDecayPdgIdVec.at ( genId ),
                                       ( genDecayLVec.at ( genId ) ),
                                       genelectronpfActivity,
                                       elesFlagVeto,
                                       elesLVec,
                                       elespfActivity,
                                       elesMiniIso
                                     );

            if( myLostElectronObj.isEl )
            {
	      ngenel++;
              int njetsbin_number = Set_njetsbin_number(njets30);
            
              myAccRecoIsoEffs.nels[njetsbin_number]+=thisweight;
              myAccRecoIsoEffs.nels_MC[njetsbin_number]++;
            
              if( myLostElectronObj.passAcc )
              {
                myAccRecoIsoEffs.nels_acc[njetsbin_number]+=thisweight;
                myAccRecoIsoEffs.nels_acc_MC[njetsbin_number]++;

                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);

                myAccRecoIsoEffs.nels_acc_bin[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nels_acc_bin_MC[ptbin_number][acbin_number]++;
              }

              //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
              if( myLostElectronObj.passId )
              {
                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);

                myAccRecoIsoEffs.nels_reco[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nels_reco_MC[ptbin_number][acbin_number]++;

                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostElectronObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostElectronObj.reco_activity);

                myAccRecoIsoEffs.nels_reco_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight;
                myAccRecoIsoEffs.nels_reco_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;

                (myBaseHistgram.h_id_genactivity_els)->Fill(myLostElectronObj.gen_activity,thisweight);
                (myBaseHistgram.h_id_recoactivity_els)->Fill(myLostElectronObj.reco_activity,thisweight);
              }

              if( myLostElectronObj.passIso )
              {
                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);

                myAccRecoIsoEffs.nels_iso[ptbin_number][acbin_number]+=thisweight;
                myAccRecoIsoEffs.nels_iso_MC[ptbin_number][acbin_number]++;
              
                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostElectronObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostElectronObj.reco_activity);
                
                myAccRecoIsoEffs.nels_iso_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight;
                myAccRecoIsoEffs.nels_iso_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
              }
              //check warning function when we calculate the efficienies!
              if ( nElectrons == 0 && myLostElectronObj.passIso )
              { 
		std::cout << "Warning: el pass iso but nElectrons=0!" << std::endl;
		std::cout << "reco_els_pt = " << myLostElectronObj.reco_pt << std::endl;
		std::cout << "reco_els_eta = " << myLostElectronObj.reco_eta << std::endl;
              }
            }//if the gen particle is muon
          }//loop gen electrons/muons

	  // dilepton computation
	  if (nElectrons == 0)
	  {
	    if (ngenel==1) myAccRecoIsoEffs.nevents_single_els+=thisweight;
	    if (ngenel==2) myAccRecoIsoEffs.nevents_di_els+=thisweight;
	  }


        }//if no muons
      
        //loop over muon CS, mtW correction factor calculation and other calculations
        if (nElectrons == 0 && nMuons == 1)
        {
          //mtw correction factor calculation
          std::vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
          std::vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");
	  std::vector<int> muonsFlagMedium = tr.getVec<int>("muonsFlagMedium");

          double reco_mus_pt = 0, reco_mus_eta = 0, reco_mus_phi = 0;
 	  int nisomuons=0;

          for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
          {
            if(muonsFlagMedium[im] && muonsLVec[im].Pt()>(AnaConsts::muonsMiniIsoArr).minPt && fabs(muonsLVec[im].Eta()) < (AnaConsts::muonsMiniIsoArr).maxAbsEta && muonsMiniIso[im] < (AnaConsts::muonsMiniIsoArr).maxIso )
            {
              reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
              reco_mus_eta = ( muonsLVec.at(im) ).Eta();
              reco_mus_phi = ( muonsLVec.at(im) ).Phi();
	      ++nisomuons;
            }
          }
	  if (nisomuons!=1) std::cout << "Error: nisomuons != 1: nisomuons = " << nisomuons << std::endl;

          double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
          double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * ( 1.0 - cos(deltaphi_mus) ) );
          (myBaseHistgram.h_mtw_mus)->Fill(mtW_mus,thisweight);

          int ptbin_number_allreco = Set_ptbin_number(reco_mus_pt);

          myAccRecoIsoEffs.mtwall[ptbin_number_allreco]+=thisweight;
          myAccRecoIsoEffs.mtwall_MC[ptbin_number_allreco]+=thisweight*thisweight;
          if( mtW_mus < 100.0 )
          {
            myAccRecoIsoEffs.mtw100[ptbin_number_allreco]+=thisweight;
            myAccRecoIsoEffs.mtw100_MC[ptbin_number_allreco]+=thisweight*thisweight;
          }

          //muon CS statistics
          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbin_id]+=thisweight;
          }
        }
      }//baseline, nolepveto
    }//TTjets samples class
  }//end of first loop

  myAccRecoIsoEffs.NumberstoEffs();
  myAccRecoIsoEffs.EffsPlotsGen();
  //myAccRecoIsoEffs.EffstoWeights();
  myAccRecoIsoEffs.GetDiLeptonFactor();
  myAccRecoIsoEffs.printAccRecoIsoEffs();
  myAccRecoIsoEffs.printEffsHeader();

  (myBaseHistgram.oFile)->Write();
  (myBaseHistgram.oFile)->Close();

  return ;
}

void LoopLLExp( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight )
{
  ClosureHistgram myClosureHistgram;
  myClosureHistgram.BookHistgram("ExpLL.root");
  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  std::cout << "Expectation: " << std::endl;
  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  {  
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    AnaFunctions::prepareTopTagger();
    //The passBaseline is registered here
    tr.registerFunction(&mypassBaselineFunc);    
 
    double thisweight = (*iter_TTJetsSampleInfos).weight;
    std::cout << "Weight " << thisweight << std::endl;
    int neventc=0;

    //while(tr.getNextEvent() && neventc<10000)
    while(tr.getNextEvent())
    {
      ++neventc;
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;
    
      myAccRecoIsoEffs.nevents_tot+=thisweight;

      //baseline cut without lepton veto
      bool passBaselinelostlept = tr.getVar<bool>("passBaseline"+spec);

      if ( 
          passBaselinelostlept 
         )
      {
        myAccRecoIsoEffs.nevents_sel_base+=thisweight;

        int nElectrons = tr.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = tr.getVar<int>("nMuons_CUT"+spec);

        double met = tr.getVar<double>("met");
        double metphi = tr.getVar<double>("metphi");
        int njets30 = tr.getVar<int>("cntNJetsPt30Eta24"+spec);
        const double ht = tr.getVar<double>("ht");
        int ntopjets = tr.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = tr.getVar<int>("cntCSVS"+spec);
        double MT2 = tr.getVar<double>("best_had_brJet_MT2"+spec);
        double mht = tr.getVar<double>("mht");

        const int nIsoTrks = tr.getVar<int>("nIsoTrks_CUT"+spec);
        //std::cout << "nIsoTrks = " << nIsoTrks << std::endl;

        //merge electron, muon generation level information
        std::vector<int> W_emuVec = tr.getVec<int>("W_emuVec");
        std::vector<int> W_tau_emuVec = tr.getVec<int>("W_tau_emuVec");
        std::vector<int> emuVec_merge;
        emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
        emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
        emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
        int gen_emus_count = emuVec_merge.size();

        std::vector<double> W_emu_pfActivityVec = tr.getVec<double>("W_emu_pfActivityVec");
        std::vector<double> W_tau_emu_pfActivityVec = tr.getVec<double>("W_tau_emu_pfActivityVec");
        std::vector<double> emu_pfActivityVec_merge;
        emu_pfActivityVec_merge.reserve( W_emu_pfActivityVec.size() + W_tau_emu_pfActivityVec.size() );
        emu_pfActivityVec_merge.insert( emu_pfActivityVec_merge.end(), W_emu_pfActivityVec.begin(), W_emu_pfActivityVec.end() );
        emu_pfActivityVec_merge.insert( emu_pfActivityVec_merge.end(), W_tau_emu_pfActivityVec.begin(), W_tau_emu_pfActivityVec.end() );

        int ngenmunotiso = 0, ngenmunotid = 0, ngenmuoutacc = 0;
        int ngenelnotiso = 0, ngenelnotid = 0, ngeneloutacc = 0;

        int ngenmu = 0;
        int ngenel = 0;

        if(nElectrons == 0)
        {
          myAccRecoIsoEffs.nevents_sel_mus+=thisweight;

          //get gen level information of leptons
          std::vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          std::vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          //get reco level information of muons
          std::vector<TLorentzVector> muonsLVec = tr.getVec<TLorentzVector>("muonsLVec");
          std::vector<double> muonspfActivity = tr.getVec<double>("muonspfActivity");
          std::vector<int> muonsFlagMedium = tr.getVec<int>("muonsFlagMedium");
          std::vector<double> muonsMiniIso = tr.getVec<double>("muonsMiniIso");

          for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
          {
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            double genmuonpfActivity = emu_pfActivityVec_merge.at(gen_emus_i);

            LostLeptonObj myLostMuonObj;
            
            myLostMuonObj.SetMyLL(
                                   //LL variables that will be always useful 
                                   genDecayPdgIdVec.at ( genId ),
                             	   ( genDecayLVec.at ( genId ) ),
                                   genmuonpfActivity,
                                   muonsFlagMedium,
                                   muonsLVec,
                                   muonspfActivity,
                                   muonsMiniIso
                                 );

            if( myLostMuonObj.isMu ) 
            { 
              ngenmu++;
              if( myLostMuonObj.passAcc ) 
              {
                if( myLostMuonObj.passId ) 
                {
                  if( myLostMuonObj.passIso ){ }
                  else ngenmunotiso++;
                }
                else ngenmunotid++;
              }
              else ngenmuoutacc++;
            }
          }//loop gen electrons/muons
        }//if no electrons

        if(nMuons == 0)
        {
          //get gen level information of leptons
          std::vector<int> genDecayPdgIdVec = tr.getVec<int>("genDecayPdgIdVec");
          std::vector<TLorentzVector> genDecayLVec = tr.getVec<TLorentzVector>("genDecayLVec");
          //get reco level information of muons
          std::vector<TLorentzVector> elesLVec = tr.getVec<TLorentzVector>("elesLVec");
          std::vector<double> elespfActivity = tr.getVec<double>("elespfActivity");
          std::vector<int> elesFlagVeto = tr.getVec<int>("elesFlagVeto");
          std::vector<double> elesMiniIso = tr.getVec<double>("elesMiniIso");

          for(int gen_emus_i = 0 ; gen_emus_i < gen_emus_count ; gen_emus_i++)
          {
            int genId;
            genId = emuVec_merge.at ( gen_emus_i );
            double genelectronpfActivity = emu_pfActivityVec_merge.at(gen_emus_i);

            LostLeptonObj myLostElectronObj;

            myLostElectronObj.SetMyLL(
                                       //LL variables that will be always useful 
                                       genDecayPdgIdVec.at ( genId ),
                                       ( genDecayLVec.at ( genId ) ),
                                       genelectronpfActivity,
                                       elesFlagVeto,
                                       elesLVec,
                                       elespfActivity,
                                       elesMiniIso
                                     );

            if( myLostElectronObj.isEl )
            { 
              ngenel++;
              if( myLostElectronObj.passAcc )
              { 
                if( myLostElectronObj.passId )
                { 
                  if( myLostElectronObj.passIso ){ }
                  else ngenelnotiso++;
                }
                else ngenelnotid++;
              }
              else ngeneloutacc++;
            }
          }//loop gen electrons/muons
        }//if no muons
      
        ///////////////////////////
        // expectation computation
        ///////////////////////////

        // exp 1 muon not iso
        if ( nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmunotiso==1) || ( ngenmu==2 && ngenmunotiso==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_iso_mus+=thisweight; 
    
	  if (!applyisotrkveto || nIsoTrks==0)
  	  {
	    (myClosureHistgram.h_exp_mu_iso_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_mu_iso_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_mu_iso_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_mu_iso_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_mu_iso_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_mu_iso_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 muon not id
        if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmunotid==1) || ( ngenmu==2 && ngenmunotid==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_id_mus+=thisweight;

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_id_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_mu_id_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_mu_id_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_mu_id_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_mu_id_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_mu_id_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 muon out acc
        if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmuoutacc==1) || ( ngenmu==2 && ngenmuoutacc==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_acc_mus+=thisweight;

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_acc_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_mu_acc_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_mu_acc_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_mu_acc_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_mu_acc_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_mu_acc_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 muon tot
        //if (nElectrons == 0 && nMuons==0 && ngenmu==1 && (ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1))
        if (nElectrons == 0 && nMuons==0 && ngenmu==1)
        {
	  if (!(ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1)) std::cout << "ngenmuoutacc = " << ngenmuoutacc << " , ngenmunotid = " << ngenmunotid << " , ngenmunotiso = " << ngenmunotiso << std::endl;
          myAccRecoIsoEffs.nevents_exp_all_mus+=thisweight;
          //myAccRecoIsoEffs.nevents_single_mus+=thisweight;

          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]+=thisweight*thisweight;
            myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[searchbin_id]+=thisweight;
          }
	
	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_musingle_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_musingle_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_musingle_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_musingle_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_musingle_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_musingle_all_ntopjets)->Fill(ntopjets,thisweight);

	    (myClosureHistgram.h_exp_mu_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_mu_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_mu_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_mu_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_mu_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 2 muons tot
        //if ( nElectrons == 0 && nMuons==0 && ngenmu==2 && ( ngenmuoutacc==2 || ngenmunotid==2 || ngenmunotiso==2 || ( ngenmuoutacc==1 && ngenmunotid==1 ) || (ngenmuoutacc==1 && ngenmunotiso==1 ) || ( ngenmunotiso==1 && ngenmunotid==1 ) ) )
        if ( nElectrons == 0 && nMuons==0 && ngenmu==2 )
        {
	  if (!( ngenmuoutacc==2 || ngenmunotid==2 || ngenmunotiso==2 || ( ngenmuoutacc==1 && ngenmunotid==1 ) || (ngenmuoutacc==1 && ngenmunotiso==1 ) || ( ngenmunotiso==1 && ngenmunotid==1 ) )) std::cout << "Warning in nElectrons == 0 && nMuons==0 && ngenmu==2" << std::endl;
          //myAccRecoIsoEffs.nevents_di_mus+=thisweight;
          myAccRecoIsoEffs.nevents_exp_all_mus+=thisweight;

          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]+=thisweight*thisweight;
            myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[searchbin_id]+=thisweight;
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_mu_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_mu_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_mu_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_mu_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp mu+ele
        if ( nElectrons == 0 && nMuons==0 && (ngenmu==1 || ngenmu==2 || ngenel==1 || ngenel==2) )
        {
	  //n_exp_lep_noitv+=thisweight;
	  //if (ngenmu==1 || ngenmu==2) ++n_exp_mu_noitv;
	  //if (ngenel==1 || ngenel==2) ++n_exp_ele_noitv;
	  //if (nIsoTrks==0) n_exp_lep_itv+=thisweight;

          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_lept_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_lept_exp_SB_MC[searchbin_id]+=thisweight*thisweight;
            myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[searchbin_id]+=thisweight;
            if (nIsoTrks==0) { myAccRecoIsoEffs.nevents_lept_exp_SB_MC_isotrk[searchbin_id]+=thisweight*thisweight; myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized_isotrk[searchbin_id]+=thisweight; }
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_lept_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_lept_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_lept_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_lept_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_lept_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_lept_all_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 electron not iso
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotiso==1) || ( ngenel==2 && ngenelnotiso==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_iso_els+=thisweight;

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_iso_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_el_iso_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_el_iso_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_el_iso_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_el_iso_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_el_iso_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 electron not id
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotid==1) || ( ngenel==2 && ngenelnotid==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_id_els+=thisweight;

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_id_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_el_id_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_el_id_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_el_id_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_el_id_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_el_id_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 electron not acc
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngeneloutacc==1) || ( ngenel==2 && ngeneloutacc==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_acc_els+=thisweight;

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_acc_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_el_acc_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_el_acc_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_el_acc_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_el_acc_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_el_acc_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }

        // exp 1 electron tot
        //if (nElectrons == 0 && nMuons==0 && ngenel==1 && (ngeneloutacc==1 || ngenelnotid==1 || ngenelnotiso==1))
        if (nElectrons == 0 && nMuons==0 && ngenel==1)
        {
	  if (!(ngeneloutacc==1 || ngenelnotid==1 || ngenelnotiso==1)) std::cout << "FSL: 1 ele tot warning" << std::endl;

          myAccRecoIsoEffs.nevents_exp_all_els+=thisweight;
          //myAccRecoIsoEffs.nevents_single_els+=thisweight;

          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]+=thisweight*thisweight;
            myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[searchbin_id]+=thisweight;
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_elsingle_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_elsingle_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_elsingle_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_elsingle_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_elsingle_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_elsingle_all_ntopjets)->Fill(ntopjets,thisweight);  

	    (myClosureHistgram.h_exp_el_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_el_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_el_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_el_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_el_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }  

        //if ( nElectrons == 0 && nMuons==0 && ngenel==2 && ( ngeneloutacc==2 || ngenelnotid==2 || ngenelnotiso==2 || ( ngeneloutacc==1 && ngenelnotid==1 ) || (ngeneloutacc==1 && ngenelnotiso==1 ) || ( ngenelnotiso==1 && ngenelnotid==1 ) ) )
        if (nElectrons == 0 && nMuons==0 && ngenel==2)
        {
          myAccRecoIsoEffs.nevents_exp_all_els+=thisweight;
          //myAccRecoIsoEffs.nevents_di_els+=thisweight;

          int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]+=thisweight*thisweight;
            myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[searchbin_id]+=thisweight;
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_all_met)->Fill(met,thisweight);
	    (myClosureHistgram.h_exp_el_all_njets)->Fill(njets30,thisweight);
	    (myClosureHistgram.h_exp_el_all_mt2)->Fill(MT2,thisweight);
	    (myClosureHistgram.h_exp_el_all_ht)->Fill(ht,thisweight);
	    (myClosureHistgram.h_exp_el_all_mht)->Fill(mht,thisweight);
	    (myClosureHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets,thisweight);
	  }
        }
      }//baseline, nolepveto
    }//TTjets samples class
  }//end of first loop

  // printSearchBin
//  //std::cout << "Muon CS # in Search Bins: " << std::endl;  
//  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
//  {
//    //myAccRecoIsoEffs.nevents_lept_exp_SB_MC[i_cal] = nevents_mus_exp_SB_MC[i_cal] + nevents_els_exp_SB_MC[i_cal];
//    //nevents_lept_pred_SB_MC[i_cal] = nevents_mus_pred_SB_MC[i_cal] + nevents_els_pred_SB_MC[i_cal];
//    //myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[i_cal] = nevents_mus_exp_SB_Normalized[i_cal] + nevents_els_exp_SB_Normalized[i_cal];
//    //myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] = nevents_mus_pred_SB_Normalized[i_cal] + nevents_els_pred_SB_Normalized[i_cal];
//    //nevents_lept_pred_SB_Normalized[i_cal] = nevents_lept_pred_SB_MC[i_cal]*scale;
//
//    //std::cout << "idx: " << i_cal << "; MC Numbers: " << nevents_mus_CS_SB_MC[i_cal] << "; Normalized Numbers: " << nevents_mus_CS_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Mus Exp MC Numbers: " << nevents_mus_exp_SB_MC[i_cal] << "; Mus Exp Normalized Numbers: " << nevents_mus_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Mus Pred Normalized Numbers: " << nevents_mus_pred_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Els Exp MC Numbers: " << nevents_els_exp_SB_MC[i_cal] << "; Els Exp Normalized Numbers: " << nevents_els_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Els Pred Normalized Numbers: " << nevents_els_pred_SB_Normalized[i_cal] << std::endl;
//    std::cout << "Lept Exp MC Numbers: " << myAccRecoIsoEffs.nevents_lept_exp_SB_MC[i_cal] << "; Lept Exp Normalized Numbers: " << myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Lept Pred Normalized Numbers: " << nevents_lept_pred_SB_Normalized[i_cal] << std::endl;
//  }
  
//  TH1D * h_cs_mus_sb = new TH1D("h_cs_mus_sb","",NSEARCH_BINS+1,0,NSEARCH_BINS+1);
//
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //double tmpscale = 0.01;
//    h_cs_mus_sb->SetBinContent( i_cal+1 , nevents_mus_CS_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_mu_sb->SetBinError(i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_exp_SB_MC[i_cal]));
    //myClosureHistgram.h_pred_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_mu_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC_err[i_cal])*tmpscale);
    myClosureHistgram.h_exp_el_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_el_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_els_exp_SB_MC[i_cal]));
    //myClosureHistgram.h_pred_el_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_el_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal])*tmpscale );
    myClosureHistgram.h_exp_lept_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_lept_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_lept_exp_SB_MC[i_cal]));
    myClosureHistgram.h_exp_lept_sb_isotrk->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized_isotrk[i_cal] );
    myClosureHistgram.h_exp_lept_sb_isotrk->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_lept_exp_SB_MC_isotrk[i_cal]));
    //myClosureHistgram.h_pred_lept_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_lept_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_lept_pred_SB_MC[i_cal])*tmpscale );
    //myClosureHistgram.h_pred_lept_sb_isotrk->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isotrackvetoeff );
    //myClosureHistgram.h_pred_lept_sb_isotrk->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_lept_pred_SB_MC[i_cal]*isotrackvetoeff)*tmpscale );
  }
//
//  //cmusCS
//  TCanvas *cmusCS = new TCanvas("cmusCS","A Simple Graph Example",200,10,700,500);
//  gStyle->SetOptStat(0);
//
//  h_cs_mus_sb->SetLineColor(1);
//  h_cs_mus_sb->SetLineWidth(3);
//  h_cs_mus_sb->Draw();
//
//  const std::string titre_musCS="CMS Simulation 2015, 10 fb^{-1}, #sqrt{s} = 13 TeV";
//  TLatex *title_musCS = new TLatex(0.09770115,0.9194915,titre_musCS.c_str());
//  title_musCS->SetNDC();
//  title_musCS->SetTextSize(0.045);
//  title_musCS->Draw("same");
//
//  //TLegend* leg_musCS = new TLegend(0.6,0.75,0.85,0.85);
//  //leg_musCS->SetBorderSize(0);
//  //leg_musCS->SetTextFont(42);
//  //leg_musCS->SetFillColor(0);
//  //leg_musCS->AddEntry(h_cs_mus_sb,"Number of Muon CS","l");
//  //leg_musCS->Draw("same");
//
//  cmusCS->SaveAs( "searchbin_mus_CS.png" );
//  cmusCS->SaveAs( "searchbin_mus_CS.C" );

  (myClosureHistgram.oFile)->Write();
  (myClosureHistgram.oFile)->Close();

  return ;
}

void LoopLLPred( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight )
{
  ClosureHistgram myClosureHistgram;
  myClosureHistgram.BookHistgram("PredLL.root");
  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;
   
  std::cout << "Prediction: " << std::endl;

  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  { 
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader trCS((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    AnaFunctions::prepareTopTagger();
    //The passBaseline is registered here
    trCS.registerFunction(&mypassBaselineFunc);

    double thisweight = (*iter_TTJetsSampleInfos).weight;
    std::cout << "Weight " << thisweight << std::endl;
    int neventc=0;

    //while(trCS.getNextEvent() && neventc<10000)
    while(trCS.getNextEvent())
    {
      ++neventc;
      if(trCS.getEvtNum()%20000 == 0) std::cout << trCS.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      bool passBaselinelostlept = trCS.getVar<bool>("passBaseline"+spec);
 
      if(
         passBaselinelostlept
        )
      {
        int nElectrons = trCS.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = trCS.getVar<int>("nMuons_CUT"+spec);

        double met = trCS.getVar<double>("met");
        double metphi = trCS.getVar<double>("metphi");

        int njets30 = trCS.getVar<int>("cntNJetsPt30Eta24"+spec);
        int ntopjets = trCS.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = trCS.getVar<int>("cntCSVS"+spec);
        double MT2 = trCS.getVar<double>("best_had_brJet_MT2"+spec);
        double ht = trCS.getVar<double>("ht");
        double mht = trCS.getVar<double>("mht");

        //muon CS
        if (nElectrons == 0 && nMuons == 1)
        {
          //counting the events for muon control sample
          myAccRecoIsoEffs.nevents_cs_mus++;
          //get reco muon variables
  	  std::vector<TLorentzVector> muonsLVec = trCS.getVec<TLorentzVector>("muonsLVec");
          std::vector<double> muonspfActivity = trCS.getVec<double>("muonspfActivity");
          std::vector<int> muonsFlagMedium = trCS.getVec<int>("muonsFlagMedium");
          std::vector<double> muonsMiniIso = trCS.getVec<double>("muonsMiniIso");

          double reco_mus_pt = -1, reco_mus_eta = 0, reco_mus_phi = 0;
          double activity = -1;
	  int nisomuons=0;
          for(unsigned int im = 0 ; im < muonsLVec.size() ; im++)
          {
            if(muonsFlagMedium[im] && muonsLVec[im].Pt()>(AnaConsts::muonsMiniIsoArr).minPt && fabs(muonsLVec[im].Eta()) < (AnaConsts::muonsMiniIsoArr).maxAbsEta && muonsMiniIso[im] < (AnaConsts::muonsMiniIsoArr).maxIso )
	    {
              reco_mus_pt  = ( muonsLVec.at(im) ).Pt();
              reco_mus_eta = ( muonsLVec.at(im) ).Eta();
              reco_mus_phi = ( muonsLVec.at(im) ).Phi();
              activity = muonspfActivity.at(im);
	      ++nisomuons;
	    }
	  }
	  if (nisomuons!=1) std::cout << "Error: nisomuons!=1: nisomuons = " << nisomuons << std::endl;
          //if ( reco_mus_pt < 0 ) continue;

          double deltaphi_mus = DeltaPhi( reco_mus_phi , metphi );
          double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * met * ( 1.0 - cos(deltaphi_mus) ) );

          if ( mtW_mus < 100.0 )
          {
	    //////////////////////////
	    // prediction computation
	    //////////////////////////
	    double EventWeight_mus = thisweight;
            int njetsbin_number = Set_njetsbin_number(njets30);
            int ptbin_number = Set_ptbin_number(reco_mus_pt);
            int acbin_number = Set_acbin_number(activity);
            int searchbin_id = find_Binning_Index( nbottomjets , ntopjets , MT2, met );
	    myAccRecoIsoEffs.EffstoWeights_fromH();

	    //std::cout << "ttbar_mtwcorrfactor[0] = " << ttbar_mtwcorrfactor[0] << std::endl;
	    //mtwcorrfactor
	    EventWeight_mus = EventWeight_mus * ttbar_mtwcorrfactor[ptbin_number];
	    //dimuon correction factor
            EventWeight_mus = EventWeight_mus * ttbar_corrfactor_di_mus;

	    if (applyisotrkveto) EventWeight_mus *= isotrackvetoeff;

            //muon prediction from muon CS
	    //Fill muon iso closure plots
	    (myClosureHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    //Fill muon id closure plots
	    (myClosureHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    //Fill muon acc closure plots
	    (myClosureHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus);
	    //Fill all muon closure plots
	    double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number];

	    (myClosureHistgram.h_pred_mu_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

	    (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

            //total events flow for muons, prediction
            myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;

            if( searchbin_id >= 0 )
            {
              myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[searchbin_id] += EventWeight_all_mus*EventWeight_mus;
              //myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_mus;

            //begin to predict lost electrons from muon CS
            double EventWeight_els = thisweight;
            //mtwcorrfactor
            EventWeight_els = EventWeight_els * ttbar_mtwcorrfactor[ptbin_number];
            //dielectron correction factor
            EventWeight_els = EventWeight_els * ttbar_corrfactor_di_els;
  
	    if (applyisotrkveto) EventWeight_els *= isotrackvetoeff;

            //electron prediction from muon CS
            //Fill electron iso closure plots
            (myClosureHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            //Fill electron id closure plots
            (myClosureHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            //Fill electron acc closure plots
            (myClosureHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els);
            //Fill all electron closure plots
            double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number] + myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number];

            (myClosureHistgram.h_pred_el_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);

            (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);

            //total events flow for electrons, prediction
            myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;

            if( searchbin_id >= 0 )
            {
              myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[searchbin_id] += EventWeight_all_els*EventWeight_els;
              //myAccRecoIsoEffs.nevents_els_pred_SB_MC[searchbin_id] += EventWeight_all_els*EventWeight_els/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_els_pred_SB_MC[searchbin_id] += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_els_err += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number]*EventWeight_els;
          }//mtW5_els<125 (muon CS)
        }//nElectrons == 0 && nMuons == 1 (muon CS)
      }//baseline_nolepveto
    }//TTJets samples class
  }//end of second loop

  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
//    h_cs_mus_sb->SetBinContent( i_cal+1 , nevents_mus_CS_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_mu_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal]));
    myClosureHistgram.h_pred_el_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_el_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]));

    myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] = myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] + myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal];

    myClosureHistgram.h_pred_lept_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_lept_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]));
    myClosureHistgram.h_pred_lept_sb_isotrk->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isotrackvetoeff );
    myClosureHistgram.h_pred_lept_sb_isotrk->SetBinError( i_cal+1 , (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isotrackvetoeff);
  }

  (myClosureHistgram.oFile)->Write();
  (myClosureHistgram.oFile)->Close();
  return ;
}


int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr <<"Please give 2 arguments " << "runListCal " << " " << "runListExpPred" << std::endl;
    std::cerr <<" Valid configurations are " << std::endl;
    std::cerr <<" ./LostLepton_MuCS_TTbar runlist_ttjets_cal.txt runlist_sample_exp_pred.txt" << std::endl;
    return -1;
  }
  const char *inputFileList_Cal = argv[1];
  const char *inputFileList_Exp_Pred = argv[2];

  //define my AccRecoIsoEffs class to stroe counts and efficiencies
  AccRecoIsoEffs myAccRecoIsoEffs;

  TTJetsSampleWeight myTTJetsSampleWeight;
  double W_Lept_BR = 0.1086*3;
  double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
  double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2
  //TTJets nominal
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_", 831.76, 11339232, LUMI, inputFileList_Cal );
  //TTJets single lepton and di-lepton
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromT_", 831.76*0.5*TTbar_SingleLept_BR, 60144642, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromTbar_", 831.76*0.5*TTbar_SingleLept_BR, 59816364, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_DiLept_", 831.76*TTbar_DiLept_BR, 30498962, LUMI, inputFileList_Cal );

  //TTJetsSampleWeight myExpPredSampleWeight;
  //myExpPredSampleWeight.TTJetsSampleInfo_push_back( "TTJets_", 831.76, 11339232, LUMI, inputFileList_Exp_Pred );

  LoopLLCal( myAccRecoIsoEffs, myTTJetsSampleWeight );
  //LoopLLExp( myAccRecoIsoEffs, myTTJetsSampleWeight );
  //LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight );


  //std::cout << "main: printOverview" << std::endl;
  //myAccRecoIsoEffs.printOverview();
  //myAccRecoIsoEffs.NormalizeFlowNumber();
  // myAccRecoIsoEffs.printNormalizeFlowNumber();

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
    for (int htbinc=0;htbinc<NHT_BINS;++htbinc)
    {
      mus_acc[i_cal][htbinc] = nmus_acc[i_cal][htbinc]/nmus[i_cal][htbinc];
      mus_acc_err[i_cal][htbinc] = get_stat_Error(nmus_acc_MC[i_cal][htbinc],nmus_MC[i_cal][htbinc]);
    }
    // mus_acc_err[i_cal] = std::sqrt( get_stat_Error(nmus_acc_MC[i_cal],nmus_MC[i_cal])*get_stat_Error(nmus_acc_MC[i_cal],nmus_MC[i_cal]) + get_sys_Error(mus_acc[i_cal],0.09)*get_sys_Error(mus_acc[i_cal],0.09) );

    els_acc[i_cal] = nels_acc[i_cal]/nels[i_cal];
    //els_acc_err[i_cal] = std::sqrt( get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal])*get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal]) + get_sys_Error(els_acc[i_cal],0.09)*get_sys_Error(els_acc[i_cal],0.09) );
    els_acc_err[i_cal] = get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal]);
  }

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      mus_recoeff[i_cal][j_cal]     = nmus_reco[i_cal][j_cal]/nmus_acc_bin[i_cal][j_cal];
      //here the error calculation is wrong?

      mus_recoeff_err[i_cal][j_cal] = get_stat_Error(nmus_reco_MC[i_cal][j_cal],nmus_acc_bin_MC[i_cal][j_cal]);
      els_recoeff[i_cal][j_cal]     = nels_reco[i_cal][j_cal]/nels_acc_bin[i_cal][j_cal];
      els_recoeff_err[i_cal][j_cal] = get_stat_Error(nels_reco_MC[i_cal][j_cal],nels_acc_bin_MC[i_cal][j_cal]);

      mus_isoeff[i_cal][j_cal]     = nmus_iso[i_cal][j_cal]/nmus_reco[i_cal][j_cal];
      mus_isoeff_err[i_cal][j_cal] = get_stat_Error(nmus_iso_MC[i_cal][j_cal],nmus_reco_MC[i_cal][j_cal]);
      els_isoeff[i_cal][j_cal]     = nels_iso[i_cal][j_cal]/nels_reco[i_cal][j_cal];
      els_isoeff_err[i_cal][j_cal] = get_stat_Error(nels_iso_MC[i_cal][j_cal],nels_reco_MC[i_cal][j_cal]);

      mus_isoeff_allreco[i_cal][j_cal]     = nmus_iso_allreco[i_cal][j_cal]/nmus_reco_allreco[i_cal][j_cal];
      mus_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nmus_iso_MC_allreco[i_cal][j_cal],nmus_reco_MC_allreco[i_cal][j_cal]);
      els_isoeff_allreco[i_cal][j_cal]     = nels_iso_allreco[i_cal][j_cal]/nels_reco_allreco[i_cal][j_cal];
      els_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nels_iso_MC_allreco[i_cal][j_cal],nels_reco_MC_allreco[i_cal][j_cal]);
    }
  }
  
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    mtwcorrfactor[i_cal] = mtwall[i_cal]/mtw100[i_cal];
    mtwcorrfactor_err[i_cal] = get_stat_Error_APNOA(mtw100[i_cal],mtwall[i_cal],std::sqrt(mtw100_MC[i_cal]),std::sqrt(mtwall_MC[i_cal]-mtw100_MC[i_cal]));
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
      mus_recoeffs2d->SetBinContent( i_cal+1 , j_cal+2, mus_recoeff[i_cal][j_cal] );
      mus_isoeffs2d->SetBinContent( i_cal+1 , j_cal+2, mus_isoeff_allreco[i_cal][j_cal] );
      els_recoeffs2d->SetBinContent( i_cal+1 , j_cal+2, els_recoeff[i_cal][j_cal] );
      els_isoeffs2d->SetBinContent( i_cal+1 , j_cal+2, els_isoeff_allreco[i_cal][j_cal] );

      mus_recoeffs2d->SetBinError( i_cal+1 , j_cal+2, mus_recoeff_err[i_cal][j_cal] );
      mus_isoeffs2d->SetBinError( i_cal+1 , j_cal+2, mus_isoeff_err_allreco[i_cal][j_cal] );
      els_recoeffs2d->SetBinError( i_cal+1 , j_cal+2, els_recoeff_err[i_cal][j_cal] );
      els_isoeffs2d->SetBinError( i_cal+1 , j_cal+2, els_isoeff_err_allreco[i_cal][j_cal] );
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

  for( i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    for (int htbinc=0;htbinc<NHT_BINS;++htbinc)
    {
      mus_acc2d->SetBinContent( i_cal+1 , htbinc+1, mus_acc[i_cal][htbinc] );
      //els_recoeffs2d->SetBinContent( i_cal+1 , j_cal+2, els_recoeff[i_cal][j_cal] );

      mus_acc2d->SetBinError( i_cal+1 , htbinc+1, mus_acc_err[i_cal][htbinc] );
      //els_recoeffs2d->SetBinError( i_cal+1 , j_cal+2, els_recoeff_err[i_cal][j_cal] );
    }       
  }
  mus_acc2d->GetXaxis()->SetTitle("N jets");
  mus_acc2d->GetYaxis()->SetTitle("HT [GeV]");

  Effs2dPlots->Write();
}

void AccRecoIsoEffs::EffstoWeights()
{
//FSLtoadd//  int i_cal;
//FSLtoadd//  int j_cal;
//FSLtoadd//  int k_cal;
//FSLtoadd//
//FSLtoadd//  for(i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++)
//FSLtoadd//  {
//FSLtoadd//    for(j_cal = 0 ; j_cal < PT_BINS ; j_cal++)
//FSLtoadd//    {
//FSLtoadd//      for(k_cal = 0 ; k_cal < AC_BINS ; k_cal++)
//FSLtoadd//      {
//FSLtoadd//        //mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - mus_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal];
//FSLtoadd//        //mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1.0 - mus_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] ); 
//FSLtoadd//        //mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - mus_acc[i_cal])/mus_acc[i_cal] );
//FSLtoadd//    
//FSLtoadd//        //els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - els_acc[i_cal])/mus_acc[i_cal] );
//FSLtoadd//        //els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1 - els_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] )* (els_acc[i_cal]/mus_acc[i_cal]); 
//FSLtoadd//        //els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - els_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal] ) * (els_recoeff[j_cal][k_cal]/mus_recoeff[j_cal][k_cal])* (els_acc[i_cal]/mus_acc[i_cal]);
//FSLtoadd//
//FSLtoadd//        mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - mus_isoeff_allreco[j_cal][k_cal])/mus_isoeff_allreco[j_cal][k_cal];
//FSLtoadd//        mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * ( (1.0 - mus_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] );
//FSLtoadd//        mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - mus_acc[i_cal])/mus_acc[i_cal] );
//FSLtoadd//
//FSLtoadd//        els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - els_acc[i_cal])/mus_acc[i_cal] );
//FSLtoadd//        els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff_allreco[j_cal][k_cal]) * ( (1 - els_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] )* (els_acc[i_cal]/mus_acc[i_cal]);
//FSLtoadd//        els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - els_isoeff_allreco[j_cal][k_cal])/mus_isoeff_allreco[j_cal][k_cal] ) * (els_recoeff[j_cal][k_cal]/mus_recoeff[j_cal][k_cal])* (els_acc[i_cal]/mus_acc[i_cal]);
//FSLtoadd//      }
//FSLtoadd//    }
//FSLtoadd//  }

  return ;
}

void AccRecoIsoEffs::EffstoWeights_fromH()
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

        mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - ttbar_mus_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal];
        mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] );
        mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_acc[i_cal])/ttbar_mus_acc[i_cal] );

        els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_els_acc[i_cal])/ttbar_mus_acc[i_cal] );
        els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1 - ttbar_els_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] )* (ttbar_els_acc[i_cal]/ttbar_mus_acc[i_cal]);
        els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - ttbar_els_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal] ) * (ttbar_els_recoeff[j_cal][k_cal]/ttbar_mus_recoeff[j_cal][k_cal])* (ttbar_els_acc[i_cal]/ttbar_mus_acc[i_cal]);
      }
    }
  }

  return ;
}

void AccRecoIsoEffs::GetDiLeptonFactor()
{
  corrfactor_di_mus = ( nevents_single_mus + nevents_di_mus ) / ( nevents_single_mus + 2 * nevents_di_mus );
  corrfactor_di_els = ( nevents_single_els + nevents_di_els ) / ( nevents_single_els + 2 * nevents_di_els );

  const double alpha = 1-0.6827;

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

  std::cout << "mtW correction factor & ";
  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
  {
    std::cout << mtwcorrfactor[i_cal] << "$\\pm$" << mtwcorrfactor_err[i_cal] << " & ";
    if( i_cal == PT_BINS-1 )
    {
      std::cout <<  " \\\\" << std::endl;
    }
  }

  std::cout << std::endl << "Muon information: " << std::endl;

//  std::cout << "number of muons from top (njets bins): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << nmus[i_cal] << " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, accepted (njets bins): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << nmus_acc[i_cal] << " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, accepted, bins (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nmus_acc_bin[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal == PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, reconstructed (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nmus_reco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, reconstructed (allreco) (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nmus_reco_allreco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, isolated (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nmus_iso[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of muons from top, isolated (allreco) (PT_BINS,AC_BINS): " << std::endl;
//
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nmus_iso_allreco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "muons from top, acceptance (NJETS_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << mus_acc[i_cal] << "(" << mus_acc_err[i_cal] << ")"<< " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  } 
//
//  std::cout << "muons from top, reconstruction efficiency (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << mus_recoeff[i_cal][j_cal] << "(" << mus_recoeff_err[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout<<"muons from top, isolation efficiency (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << mus_isoeff[i_cal][j_cal] << "(" << mus_isoeff_err[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout<<"muons from top, isolation efficiency (allreco) (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << mus_isoeff_allreco[i_cal][j_cal] << "(" << mus_isoeff_err_allreco[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
  std::cout << "correction factor from di muons: "<< corrfactor_di_mus << "(" << corrfactor_di_mus_err << ")" << std::endl;
//
//  std::cout << std::endl << "Electron information: " << std::endl;
//
//  std::cout << "number of electrons from top (NJETS_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << nels[i_cal] << " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of electrons from top, accepted (NJETS_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << nels_acc[i_cal] << " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//
//  std::cout << "number of electrons from top, accepted, bins (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nels_acc_bin[i_cal][j_cal] <<" ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout<<std::endl;
//    }
//  }
//
//  std::cout << "number of electrons from top, reconstructed (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nels_reco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of electrons from top, reconstructed (allreco) (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nels_reco_allreco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//
//  std::cout << "number of electrons from top, isolated (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nels_iso[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "number of electrons from top, isolated (allreco) (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << nels_iso_allreco[i_cal][j_cal] << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "electrons from top, acceptance (NJETS_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal < NJETS_BINS ; i_cal++ )
//  {
//    std::cout << els_acc[i_cal] << "(" << els_acc_err[i_cal] << ")"<< " ";
//    if( i_cal == NJETS_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//
//  std::cout << "electrons from top, reconstruction efficiency (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << els_recoeff[i_cal][j_cal] << "(" << els_recoeff_err[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal == PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "electrons from top, isolation efficiency (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << els_isoeff[i_cal][j_cal] << "(" << els_isoeff_err[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }
//
//  std::cout << "electrons from top, isolation efficiency (allreco) (PT_BINS,AC_BINS): " << std::endl;
//  for( i_cal=0 ; i_cal<PT_BINS ; i_cal++ )
//  {
//    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
//    {
//      std::cout << els_isoeff_allreco[i_cal][j_cal] << "(" << els_isoeff_err_allreco[i_cal][j_cal] << ")" << " ";
//      if( j_cal == AC_BINS-1 )
//      {
//        std::cout << std::endl;
//      }
//    }
//    if( i_cal==PT_BINS-1 )
//    {
//      std::cout << std::endl;
//    }
//  }

  std::cout << "correction factor from di electrons: " << corrfactor_di_els << "(" << corrfactor_di_els_err << ")" <<std::endl;

  return ;
}


void AccRecoIsoEffs::printEffsHeader()
{
  std::ofstream EffsHeader;
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

  EffsHeader << "  const double ttbar_mus_acc[" << NJETS_BINS << "][" << NHT_BINS << "] = "; 
  for( i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++ )
  {
    for (int htbinc=0;htbinc<NHT_BINS;++htbinc)
    {
      if( i_cal == 0 && htbinc == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && htbinc == 0 ) { EffsHeader << "{"; }
      EffsHeader << mus_acc[i_cal][htbinc];
      if( htbinc != NHT_BINS-1 ) { EffsHeader << ","; }
      if( i_cal != NJETS_BINS-1 && htbinc == NHT_BINS-1) { EffsHeader << "},"; }
      if( i_cal == NJETS_BINS-1 && htbinc == NHT_BINS-1) { EffsHeader << "}};" << std::endl; }
    }
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


void LostLeptonObj::SetMyLL( 
                             int pid,
                             TLorentzVector onegenlept,
                             double genpfActivity,
                             std::vector<int> IdFlag,
                             std::vector<TLorentzVector> recoleptLVec,
                             std::vector<double> recoleptactivityVec,
                             std::vector<double> MiniIso
                           )
{


  SetFlavor( pid );
  genLeptonSetup( onegenlept, genpfActivity );
  gogoAcc();
  gogoId( IdFlag, recoleptLVec );
 
  //here is just a temporary protect against segmentation error we may have during the run, will get rid of this once we go to use the actvity variable stored in then flattrees
  if( reco_index >= 0 )
  {
    recoLeptonSetup( recoleptLVec.at(reco_index), recoleptactivityVec.at(reco_index) );
  }

  gogoIso( MiniIso );

  return ;
}

void LostLeptonObj::SetFlavor(int pid)
{
  isMu = (pid == 13) || (pid == -13);
  isEl = (pid == 11) || (pid == -11);

  if( (!isMu) && (!isEl) ) { std::cout << "The particle in the emuVec is not mu nor el, what the fuck is going on ?!" << std::endl; }

  return ;
}

void LostLeptonObj::genLeptonSetup(TLorentzVector onegenlept, double activity)
{
  gen_eta = onegenlept.Eta(); 
  gen_phi = onegenlept.Phi(); 
  gen_pt = onegenlept.Pt(); 
  gen_activity = activity;
}

void LostLeptonObj::gogoAcc()
{
  if( isMu ) passAcc = (std::abs(gen_eta)) < (AnaConsts::muonsMiniIsoArr).maxAbsEta && gen_pt > (AnaConsts::muonsMiniIsoArr).minPt;
  if( isEl ) passAcc = (std::abs(gen_eta)) < (AnaConsts::elesMiniIsoArr).maxAbsEta && gen_pt > (AnaConsts::elesMiniIsoArr).minPt;

  doneAcc = true;
  return ;
}

void LostLeptonObj::gogoId( std::vector<int> IdFlag, std::vector<TLorentzVector> recoleptLVec )
{
  if( (doneAcc && !passAcc) ){ return ; }

  //loop over reco lepton information to find out smallest deltaR value, but those muon must be medium or tight muon first
  int reco_count = recoleptLVec.size();
  std::vector<double> deltar_pool;
  for(int reco_i = 0 ; reco_i < reco_count ; reco_i++)
  {
    //make sure the reco lepton also pass the acc requirement
    if( isMu && !( (recoleptLVec.at(reco_i)).Pt()>(AnaConsts::muonsMiniIsoArr).minPt && std::abs((recoleptLVec.at(reco_i)).Eta())<(AnaConsts::muonsMiniIsoArr).maxAbsEta ) )  continue;
    if( isEl && !( (recoleptLVec.at(reco_i)).Pt()>(AnaConsts::elesMiniIsoArr).minPt && std::abs((recoleptLVec.at(reco_i)).Eta())<(AnaConsts::elesMiniIsoArr).maxAbsEta ) )  continue;

    if ( IdFlag.at(reco_i) )
    {
      double deltar_media;
      deltar_media = DeltaR(
                            gen_eta,
		            gen_phi,
			    (recoleptLVec.at(reco_i)).Eta(),
			    (recoleptLVec.at(reco_i)).Phi()
			   );

      deltar_pool.push_back(deltar_media);
    }
  }

  double deltar;
  deltar = 1000;

  if(deltar_pool.size() > 0)
  {
    deltar = *(std::min_element(deltar_pool.begin(), deltar_pool.end()));
    reco_index = std::min_element(deltar_pool.begin(), deltar_pool.end()) - deltar_pool.begin();
  }

  deltar_pool.clear();

  passId = (deltar < 0.2);

  doneId = true;
  return ;
}

void LostLeptonObj::recoLeptonSetup(TLorentzVector onerecolept, double activity)
{
  reco_eta = onerecolept.Eta(); 
  reco_phi = onerecolept.Phi(); 
  reco_pt = onerecolept.Pt(); 
  reco_activity = activity;
}

void LostLeptonObj::gogoIso(std::vector<double> iso)
{
  if( (doneAcc && !passAcc) || (doneId && !passId) ){ return ; }

  if( reco_index < 0) std::cout << "Segmentation voilation protection do not work! Please think hard!" << std::endl;   

  if( isMu ) passIso = ( iso.at(reco_index) < (AnaConsts::muonsMiniIsoArr).maxIso );
  if( isEl ) passIso = ( iso.at(reco_index) < (AnaConsts::elesMiniIsoArr).maxIsoEB );

  doneIso = true;
  return ; 
}
