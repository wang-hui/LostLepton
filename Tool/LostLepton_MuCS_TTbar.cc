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
//#include "v151203_EffsHeader_MuCS.h"
//#include "v151209_EffsHeader_MuCS.h"
//#include "v160105_EffsHeader_MuCS.h"
//#include "v160114_EffsHeader_MuCS.h"
//#include "v160114_emuaccSB_EffsHeader_MuCS.h"
//#include "v160114_dimuonStatUp_EffsHeader_MuCS.h"
//#include "v160114_dimuonSystUp_EffsHeader_MuCS.h"
//#include "v160114_dieleStatUp_EffsHeader_MuCS.h"
//#include "v160114_dieleSystUp_EffsHeader_MuCS.h"
//#include "v160208_muIsoSystUp_EffsHeader_MuCS.h"
//#include "v160208_muIsoStatUp_EffsHeader_MuCS.h"
//#include "v160208_muIsoStatDown_EffsHeader_MuCS.h"
//#include "v160208_muIsoSystDown_EffsHeader_MuCS.h"
//#include "v160208_muRecoSystUp_EffsHeader_MuCS.h"
//#include "v160208_muRecoStatUp_EffsHeader_MuCS.h"
//#include "v160208_muRecoStatDown_EffsHeader_MuCS.h"
//#include "v160208_muRecoSystDown_EffsHeader_MuCS.h"
//#include "v160208_elIsoStatUp_EffsHeader_MuCS.h"
//#include "v160208_elIsoStatDown_EffsHeader_MuCS.h"
//#include "v160208_elIsoSystUp_EffsHeader_MuCS.h"
//#include "v160208_elIsoSystDown_EffsHeader_MuCS.h"
//#include "v160208_elIdStatUp_EffsHeader_MuCS.h"
//#include "v160208_elIdStatDown_EffsHeader_MuCS.h"
//#include "v160208_elIdSystUp_EffsHeader_MuCS.h"
//#include "v160208_elIdSystDown_EffsHeader_MuCS.h"
//#include "v160208_istrkUp_EffsHeader_MuCS.h"
//#include "v160122_accfromInvertedDPhi_EffsHeader_MuCS.h"
//#include "v160122_invertedDPhi_EffsHeader_MuCS.h"
//#include "v160217_ttbar_EffsHeader_MuCS.h"
//#include "v160217_ttbar_v3_EffsHeader_MuCS.h"
//#include "v160217_ttbar_v3_invertedDPhi_EffsHeader_MuCS.h"
//#include "v160224_ttbarSingleTopWjets_EffsHeader_MuCS.h"
//#include "v160302_ttbarv6_EffsHeader_MuCS.h"
//#include "v160302_ttbarSTWv6_EffsHeader_MuCS.h"
//#include "v160302_muIsoStatUp_EffsHeader_MuCS.h"
//#include "v160302_muIsoStatDown_EffsHeader_MuCS.h"
//#include "v160302_muIsoSystUp_EffsHeader_MuCS.h"
//#include "v160302_muIsoSystDown_EffsHeader_MuCS.h"
//#include "v160302_eIsoStatUp_EffsHeader_MuCS.h"
//#include "v160302_eIsoStatDown_EffsHeader_MuCS.h"
//#include "v160302_eIsoSystUp_EffsHeader_MuCS.h"
//#include "v160302_eIsoSystDown_EffsHeader_MuCS.h"
//#include "v160302_muIDStatUp_EffsHeader_MuCS.h"
//#include "v160302_muIDStatDown_EffsHeader_MuCS.h"
//#include "v160302_muIDSystUp_EffsHeader_MuCS.h"
//#include "v160302_muIDSystDown_EffsHeader_MuCS.h"
//#include "v160302_eIDStatUp_EffsHeader_MuCS.h"
//#include "v160302_eIDStatDown_EffsHeader_MuCS.h"
//#include "v160302_eIDSystUp_EffsHeader_MuCS.h"
//#include "v160302_eIDSystDown_EffsHeader_MuCS.h"
//#include "v160302_istrkUp_EffsHeader_MuCS.h"
//#include "v160309_ttbarSingletopW_45bins_EffsHeader_MuCS.h"
//#include "v160309_ttbarSingletopW_45bins_invertedDphi_EffsHeader_MuCS.h"
//#include "v2_EffsHeader_MuCS.h"
//#include "v4_59_EffsHeader_MuCS.h"
//#include "v5_59_EffsHeader_MuCS.h"
//#include "v4_59_ttbarW_EffsHeader_MuCS.h"
//#include "v7_EffsHeader_MuCS.h"
//#include "invertedDPhi_v160706_EffsHeader_MuCS.h"
//#include "v160707_invertedDPhi_EffsHeader_MuCS.h"
#include "v160707c_EffsHeader_MuCS.h"
// v7 is ttbar+Wjets after bug fix
#include "TriggerEff.h"

//#include "SusyAnaTools/Tools/PDFUncertainty.h"

//const double isotrackvetoeff = 0.563499421;
const bool applyisotrkveto = false; // should be false
//const double isotrackvetoeff = 1;
double isotrkeff[NSEARCH_BINS];

void LoopLLCal( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight )
{
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram("test.root");

  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  SearchBins theSearchBins("SB_59_2016");

  std::cout << "Efficiencies Calculation: " << std::endl;

  double neventsSB[NSEARCH_BINS];
  double neventsSB_afterITV[NSEARCH_BINS];
  double neventsSB_MC[NSEARCH_BINS];
  double neventsSB_afterITV_MC[NSEARCH_BINS];
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    neventsSB[searchbinc]=0.0;
    neventsSB_afterITV[searchbinc]=0.0;
    neventsSB_MC[searchbinc]=0.0;
    neventsSB_afterITV_MC[searchbinc]=0.0;
    isotrkeff[searchbinc]=0.0;
  }

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

      bool passLeptVeto = tr.getVar<bool>("passLeptVeto"+spec);
      bool passnJets = tr.getVar<bool>("passnJets"+spec);
      bool passMET = tr.getVar<bool>("passMET"+spec);
      bool passHT = tr.getVar<bool>("passHT"+spec);
      bool passMT2 = tr.getVar<bool>("passMT2"+spec);
      bool passTagger = tr.getVar<bool>("passTagger"+spec);
      bool passBJets = tr.getVar<bool>("passBJets"+spec);
      bool passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter"+spec);
      bool passdPhis = tr.getVar<bool>("passdPhis"+spec);

      bool passInvertedBaseline = 
	//passLeptVeto &&
                          passnJets
                          && passHT
                          && passMT2
                          && passTagger
                          && passBJets
                          && passNoiseEventFilter
	                  && passMET
	                  && !passdPhis;

      //baseline cut without lepton veto
      bool passBaselinelostlept = tr.getVar<bool>("passBaseline"+spec);

      if (passBaselinelostlept)
      //if (passInvertedBaseline)
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
	      //const int HTbin_number = Set_MT2bin_number(MT2);
	      //std::cout << "HTbin_number = " << HTbin_number << std::endl;            

              myAccRecoIsoEffs.nmus[njetsbin_number][HTbin_number]+=thisweight;
              //myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
              myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]++;

	      const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );            
	      myAccRecoIsoEffs.nmus_sb[searchbin_id]+=thisweight;
	      myAccRecoIsoEffs.nmus_MC_sb[searchbin_id]++;

	      //debug
	      //if (searchbin_id==21) std::cout << "thisweight = " << thisweight << std::endl;

              if( myLostMuonObj.passAcc )
              {
                myAccRecoIsoEffs.nmus_acc[njetsbin_number][HTbin_number]+=thisweight;
                //myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
                myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]++;

                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);

		myAccRecoIsoEffs.nmus_acc_sb[searchbin_id]+=thisweight;
		myAccRecoIsoEffs.nmus_acc_MC_sb[searchbin_id]++;

		//debug
		//if (searchbin_id==21) std::cout << "acc thisweight = " << thisweight << std::endl;

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
	      const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );   
            
              //myAccRecoIsoEffs.nels[njetsbin_number]+=thisweight;
              //myAccRecoIsoEffs.nels_MC[njetsbin_number]++;
	      myAccRecoIsoEffs.nels[searchbin_id]+=thisweight;
	      myAccRecoIsoEffs.nels_MC[searchbin_id]++;

              if( myLostElectronObj.passAcc )
              {
                //myAccRecoIsoEffs.nels_acc[njetsbin_number]+=thisweight;
                //myAccRecoIsoEffs.nels_acc_MC[njetsbin_number]++;
                myAccRecoIsoEffs.nels_acc[searchbin_id]+=thisweight;
                myAccRecoIsoEffs.nels_acc_MC[searchbin_id]++;

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
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbin_id]+=thisweight;
          }
        }

	if (nElectrons == 0 && nMuons == 0 && (ngenmu==1 || ngenmu==2 || ngenel==1 || ngenel==2))
	{
          const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
	    neventsSB[searchbin_id]+=thisweight;
	    if (nIsoTrks==0) neventsSB_afterITV[searchbin_id]+=thisweight;

	    if (neventsSB_afterITV[searchbin_id]>neventsSB[searchbin_id]) std::cout << "neventsSB_afterITV[searchbin_id]>neventsSB[searchbin_id]" << std::endl;

	    ++neventsSB_MC[searchbin_id];
	    if (nIsoTrks==0) ++neventsSB_afterITV_MC[searchbin_id];

            if (neventsSB_afterITV[searchbin_id]>neventsSB[searchbin_id]) std::cout << "neventsSB_afterITV[searchbin_id]>neventsSB[searchbin_id]" << std::endl;

	  }
	}
      }//baseline, nolepveto
    }//TTjets samples class
  }//end of first loop

  // printSearchBin
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    if (neventsSB[searchbinc]==0)
    {
      std::cout << "neventsSB[" << searchbinc << "]==0, isotrkeff[" << searchbinc << "] set to 0.5." << std::endl;
      isotrkeff[searchbinc]=0.5;
    }
    else
    {
      isotrkeff[searchbinc]=neventsSB_afterITV[searchbinc]/neventsSB[searchbinc];
      //std::cout << "neventsSB[" << searchbinc << "] = " << neventsSB[searchbinc] << " , neventsSB_afterITV = " << neventsSB_afterITV[searchbinc] << " , ratio = " << isotrkeff[searchbinc] << std::endl;
      std::cout << "isotrackeff_SB[" << searchbinc << "] = " << isotrkeff[searchbinc] << std::endl;
    }
  }
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    std::cout << "isoTrackErr[" << searchbinc << "] = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[searchbinc],neventsSB_MC[searchbinc]) << ";" << std::endl;
  }

//  // debug
//  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
//  {
//    std::cout << "nmus_acc_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_acc_sb[searchbinc] << std::endl;
//    std::cout << "nmus_acc_MC_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_acc_MC_sb[searchbinc] << std::endl;
//    std::cout << "nmus_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_sb[searchbinc] << std::endl;
//    std::cout << "nmus_MC_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_MC_sb[searchbinc] << std::endl;
//  }
//  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
//  {
//    std::cout << "nels_acc[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_acc[searchbinc] << std::endl;
//    std::cout << "nels_acc_MC[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_acc_MC[searchbinc] << std::endl;
//    std::cout << "nels[" << searchbinc << "] = " << myAccRecoIsoEffs.nels[searchbinc] << std::endl;
//    std::cout << "nels_MC[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_MC[searchbinc] << std::endl;
//  }

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

  SearchBins theSearchBins("SB_59_2016");

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


      bool passLeptVeto = tr.getVar<bool>("passLeptVeto"+spec);
      bool passnJets = tr.getVar<bool>("passnJets"+spec);
      bool passMET = tr.getVar<bool>("passMET"+spec);
      bool passHT = tr.getVar<bool>("passHT"+spec);
      bool passMT2 = tr.getVar<bool>("passMT2"+spec);
      bool passTagger = tr.getVar<bool>("passTagger"+spec);
      bool passBJets = tr.getVar<bool>("passBJets"+spec);
      bool passNoiseEventFilter = tr.getVar<bool>("passNoiseEventFilter"+spec);
      bool passdPhis = tr.getVar<bool>("passdPhis"+spec);

      //normal baseline without dPhis cut
      bool passInvertedBaseline = 
	//passLeptVeto &&
                          passnJets
                          && passHT
                          && passMT2
                          && passTagger
                          && passBJets
                          && passNoiseEventFilter
	                  && passMET
	                  && !passdPhis;


      if ( 
          passBaselinelostlept 
	  //passInvertedBaseline
         )
      {
        myAccRecoIsoEffs.nevents_sel_base+=thisweight;

        int nElectrons = tr.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = tr.getVar<int>("nMuons_CUT"+spec);

        double met = tr.getVar<double>("met");
        double metphi = tr.getVar<double>("metphi");
        int njets30 = tr.getVar<int>("cntNJetsPt30Eta24"+spec);
        int njets50 = tr.getVar<int>("cntNJetsPt50Eta24"+spec);
        //const double ht = tr.getVar<double>("ht");
        const double ht = tr.getVar<double>("HT"+spec);
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

	//const double metEff = QCDGetTriggerEff( "TTJets", met );
	const double metEff = 1.0;

        if(nElectrons == 0)
        {
          myAccRecoIsoEffs.nevents_sel_mus+=(thisweight*metEff);

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
          myAccRecoIsoEffs.nevents_exp_iso_mus+=(thisweight*metEff); 
    
	  if (!applyisotrkveto || nIsoTrks==0)
  	  {
	    (myClosureHistgram.h_exp_mu_iso_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_iso_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_iso_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_iso_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_iso_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_iso_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //if( searchbin_id >= 0 && ngenmu==1)
	  if( searchbin_id >= 0)
          {
            myAccRecoIsoEffs.nevents_mus_exp_iso_SB_Normalized[searchbin_id]+=thisweight*metEff;
            if (nIsoTrks==0) myAccRecoIsoEffs.nevents_mus_exp_iso_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
            if (nIsoTrks==0) myAccRecoIsoEffs.w2_mus_exp_iso_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }


        }

        // exp 1 muon not id
        if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmunotid==1) || ( ngenmu==2 && ngenmunotid==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_id_mus+=(thisweight*metEff);

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_id_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_id_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_id_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_id_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_id_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_id_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //if( searchbin_id >= 0  && ngenmu==1)
	  if( searchbin_id >= 0)
          {
            myAccRecoIsoEffs.nevents_mus_exp_reco_SB_Normalized[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.nevents_mus_exp_reco_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.w2_mus_exp_reco_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }
        }

        // exp 1 muon out acc
        if (nElectrons == 0 && nMuons==0 && ( (ngenmu==1 && ngenmuoutacc==1) || ( ngenmu==2 && ngenmuoutacc==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_acc_mus+=(thisweight*metEff);

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_acc_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_acc_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_acc_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_acc_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_acc_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_acc_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //if( searchbin_id >= 0  && ngenmu==1)
	  if( searchbin_id >= 0)
          {
            myAccRecoIsoEffs.nevents_mus_exp_acc_SB_Normalized[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.nevents_mus_exp_acc_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.w2_mus_exp_acc_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }
        }

        // exp 1 muon tot
        //if (nElectrons == 0 && nMuons==0 && ngenmu==1 && (ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1))
        if (nElectrons == 0 && nMuons==0 && ngenmu==1)
        {
	  if (!(ngenmuoutacc==1 || ngenmunotid==1 || ngenmunotiso==1)) std::cout << "ngenmuoutacc = " << ngenmuoutacc << " , ngenmunotid = " << ngenmunotid << " , ngenmunotiso = " << ngenmunotiso << std::endl;
          myAccRecoIsoEffs.nevents_exp_all_mus+=(thisweight*metEff);
          //myAccRecoIsoEffs.nevents_single_mus+=(thisweight*metEff);

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff);
            myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[searchbin_id]+=(thisweight*metEff);
          }
	
	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_musingle_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_musingle_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_musingle_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_musingle_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_musingle_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_musingle_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));

	    (myClosureHistgram.h_exp_mu_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }
        }

        // exp 2 muons tot
        //if ( nElectrons == 0 && nMuons==0 && ngenmu==2 && ( ngenmuoutacc==2 || ngenmunotid==2 || ngenmunotiso==2 || ( ngenmuoutacc==1 && ngenmunotid==1 ) || (ngenmuoutacc==1 && ngenmunotiso==1 ) || ( ngenmunotiso==1 && ngenmunotid==1 ) ) )
        if ( nElectrons == 0 && nMuons==0 && ngenmu==2 )
        {
	  if (!( ngenmuoutacc==2 || ngenmunotid==2 || ngenmunotiso==2 || ( ngenmuoutacc==1 && ngenmunotid==1 ) || (ngenmuoutacc==1 && ngenmunotiso==1 ) || ( ngenmunotiso==1 && ngenmunotid==1 ) )) std::cout << "Warning in nElectrons == 0 && nMuons==0 && ngenmu==2" << std::endl;
          //myAccRecoIsoEffs.nevents_di_mus+=(thisweight*metEff);
          myAccRecoIsoEffs.nevents_exp_all_mus+=(thisweight*metEff);

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_exp_SB_MC[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff);
            myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[searchbin_id]+=(thisweight*metEff);
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_mu_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_mu_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }
        }

        // exp mu+ele
        if ( nElectrons == 0 && nMuons==0 && (ngenmu==1 || ngenmu==2 || ngenel==1 || ngenel==2) )
        {
	  //n_exp_lep_noitv+=(thisweight*metEff);
	  //if (ngenmu==1 || ngenmu==2) ++n_exp_mu_noitv;
	  //if (ngenel==1 || ngenel==2) ++n_exp_ele_noitv;
	  //if (nIsoTrks==0) n_exp_lep_itv+=(thisweight*metEff);

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_lept_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_lept_exp_SB_MC[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff);
            myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[searchbin_id]+=(thisweight*metEff);
            if (nIsoTrks==0) { myAccRecoIsoEffs.nevents_lept_exp_SB_MC_isotrk[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff); myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized_isotrk[searchbin_id]+=(thisweight*metEff); }
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_lept_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_lept_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_lept_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_lept_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_lept_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_lept_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }
	  if (nIsoTrks==0)
	  {	  
  	    (myClosureHistgram.h_exp_lept_all_Z_sb)->Fill(searchbin_id,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_njets30)->Fill(njets30,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_njets50)->Fill(njets50,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_ntops)->Fill(ntopjets,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_nbjets)->Fill(nbottomjets,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_MET)->Fill(met,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_MT2)->Fill(MT2,(thisweight*metEff));
  	    (myClosureHistgram.h_exp_lept_all_Z_HT)->Fill(ht,(thisweight*metEff));
	  }
        }

        // exp 1 electron not iso
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotiso==1) || ( ngenel==2 && ngenelnotiso==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_iso_els+=(thisweight*metEff);

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_iso_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_iso_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_iso_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_iso_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_iso_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_iso_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
	  if( searchbin_id >= 0)
          {
            //myAccRecoIsoEffs.nevents_mus_exp_iso_SB_Normalized[searchbin_id]+=thisweight*metEff;
            if (nIsoTrks==0) myAccRecoIsoEffs.nevents_els_exp_iso_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
            if (nIsoTrks==0) myAccRecoIsoEffs.w2_els_exp_iso_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }


        }

        // exp 1 electron not id
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngenelnotid==1) || ( ngenel==2 && ngenelnotid==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_id_els+=(thisweight*metEff);

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_id_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_id_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_id_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_id_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_id_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_id_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
	  if( searchbin_id >= 0)
          {
            //myAccRecoIsoEffs.nevents_mus_exp_reco_SB_Normalized[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.nevents_els_exp_reco_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.w2_els_exp_reco_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }


        }

        // exp 1 electron not acc
        if (nElectrons == 0 && nMuons==0 && ( (ngenel==1 && ngeneloutacc==1) || ( ngenel==2 && ngeneloutacc==2) ) )
        {
          myAccRecoIsoEffs.nevents_exp_acc_els+=(thisweight*metEff);

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_acc_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_acc_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_acc_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_acc_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_acc_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_acc_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
	  if( searchbin_id >= 0)
          {
            //myAccRecoIsoEffs.nevents_mus_exp_acc_SB_Normalized[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.nevents_els_exp_acc_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff;
	    if (nIsoTrks==0) myAccRecoIsoEffs.w2_els_exp_acc_SB_Normalized_isotrk[searchbin_id]+=thisweight*metEff*thisweight*metEff;
          }
        }

        // exp 1 electron tot
        //if (nElectrons == 0 && nMuons==0 && ngenel==1 && (ngeneloutacc==1 || ngenelnotid==1 || ngenelnotiso==1))
        if (nElectrons == 0 && nMuons==0 && ngenel==1)
        {
	  if (!(ngeneloutacc==1 || ngenelnotid==1 || ngenelnotiso==1)) std::cout << "FSL: 1 ele tot warning" << std::endl;

          myAccRecoIsoEffs.nevents_exp_all_els+=(thisweight*metEff);
          //myAccRecoIsoEffs.nevents_single_els+=(thisweight*metEff);

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff);
            myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[searchbin_id]+=(thisweight*metEff);
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_elsingle_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_elsingle_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_elsingle_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_elsingle_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_elsingle_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_elsingle_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));  

	    (myClosureHistgram.h_exp_el_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));
	  }
        }  

        //if ( nElectrons == 0 && nMuons==0 && ngenel==2 && ( ngeneloutacc==2 || ngenelnotid==2 || ngenelnotiso==2 || ( ngeneloutacc==1 && ngenelnotid==1 ) || (ngeneloutacc==1 && ngenelnotiso==1 ) || ( ngenelnotiso==1 && ngenelnotid==1 ) ) )
        if (nElectrons == 0 && nMuons==0 && ngenel==2)
        {
          myAccRecoIsoEffs.nevents_exp_all_els+=(thisweight*metEff);
          //myAccRecoIsoEffs.nevents_di_els+=(thisweight*metEff);

          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          if( searchbin_id >= 0 )
          {
            //myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_els_exp_SB_MC[searchbin_id]+=(thisweight*metEff)*(thisweight*metEff);
            myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[searchbin_id]+=(thisweight*metEff);
          }

	  if (!applyisotrkveto || nIsoTrks==0)
	  {
	    (myClosureHistgram.h_exp_el_all_met)->Fill(met,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_njets)->Fill(njets30,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_mt2)->Fill(MT2,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_ht)->Fill(ht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_mht)->Fill(mht,(thisweight*metEff));
	    (myClosureHistgram.h_exp_el_all_ntopjets)->Fill(ntopjets,(thisweight*metEff));
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

  double nmuacc=0;
  double nmuiso=0;
  double nmuid=0;
  double neacc=0;
  double neiso=0;
  double neid=0;
  double ntot=0;
  double nmuacc_err=0;
  double nmuiso_err=0;
  double nmuid_err=0;
  double neacc_err=0;
  double neiso_err=0;
  double neid_err=0;
  double ntot_err=0;

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

    myClosureHistgram.h_exp_iso_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_exp_iso_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_reco_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_exp_reco_SB_Normalized[i_cal] );
    myClosureHistgram.h_exp_acc_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_exp_acc_SB_Normalized[i_cal] );

    nmuacc+=myAccRecoIsoEffs.nevents_mus_exp_acc_SB_Normalized_isotrk[i_cal];
    nmuiso+=myAccRecoIsoEffs.nevents_mus_exp_iso_SB_Normalized_isotrk[i_cal];
    nmuid+=myAccRecoIsoEffs.nevents_mus_exp_reco_SB_Normalized_isotrk[i_cal];
    neacc+=myAccRecoIsoEffs.nevents_els_exp_acc_SB_Normalized_isotrk[i_cal];
    neiso+=myAccRecoIsoEffs.nevents_els_exp_iso_SB_Normalized_isotrk[i_cal];
    neid+=myAccRecoIsoEffs.nevents_els_exp_reco_SB_Normalized_isotrk[i_cal];
    ntot+=myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized_isotrk[i_cal];
    nmuacc_err+=std::sqrt(myAccRecoIsoEffs.w2_mus_exp_acc_SB_Normalized_isotrk[i_cal]);
    nmuiso_err+=std::sqrt(myAccRecoIsoEffs.w2_mus_exp_iso_SB_Normalized_isotrk[i_cal]);
    nmuid_err+=std::sqrt(myAccRecoIsoEffs.w2_mus_exp_reco_SB_Normalized_isotrk[i_cal]);
    neacc_err+=std::sqrt(myAccRecoIsoEffs.w2_els_exp_acc_SB_Normalized_isotrk[i_cal]);
    neiso_err+=std::sqrt(myAccRecoIsoEffs.w2_els_exp_iso_SB_Normalized_isotrk[i_cal]);
    neid_err+=std::sqrt(myAccRecoIsoEffs.w2_els_exp_reco_SB_Normalized_isotrk[i_cal]);
    ntot_err+=std::sqrt(myAccRecoIsoEffs.nevents_lept_exp_SB_MC_isotrk[i_cal]);
  }

  std::cout << "nmuiso = " << nmuiso << " pm " << nmuiso_err << std::endl;
  std::cout << "nmuid = " << nmuid << " pm " << nmuid_err << std::endl;
  std::cout << "nmuacc = " << nmuacc << " pm " << nmuacc_err << std::endl;
  std::cout << "neiso = " << neiso << " pm " << neiso_err << std::endl;
  std::cout << "neid = " << neid << " pm " << neid_err << std::endl;
  std::cout << "neacc = " << neacc << " pm " << neacc_err << std::endl;
  std::cout << "ntot = " << ntot << " pm " << ntot_err << std::endl;

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
  //myBaselineVessel = new BaselineVessel(spec, "fastsim");

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  SearchBins theSearchBins("SB_59_2016");
   
  std::cout << "Prediction: " << std::endl;
  int sidebandLowMT2=0;
  int sidebandHighMT2=0;
  std::vector< std::pair<int, int> > viMP;
  std::map<std::pair<int, int>, std::pair<std::array<double,NSEARCH_BINS>, std::array<double,NSEARCH_BINS> > > mapSignalpred;

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

    //while(trCS.getNextEvent() && neventc<100000)
    while(trCS.getNextEvent())
    {
      ++neventc;
      if(trCS.getEvtNum()%20000 == 0) std::cout << trCS.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      bool passBaselinelostlept = trCS.getVar<bool>("passBaseline"+spec);
 
      //const double& SusyMotherMass  = trCS.getVar<double>("SusyMotherMass");
      //const double& SusyLSPMass     = trCS.getVar<double>("SusyLSPMass");

      std::vector<std::string> TriggerNames = trCS.getVec<std::string>("TriggerNames");
      std::vector<int> PassTrigger = trCS.getVec<int>("PassTrigger");
      bool foundTrigger = false;

      for(unsigned it=0; it<TriggerNames.size(); it++)
      {
	if( TriggerNames[it].find("HLT_PFHT350_PFMET100_JetIdCleaned_v") || TriggerNames[it].find("HLT_PFHT350_PFMET100_NoiseCleaned_v") || TriggerNames[it].find("HLT_PFHT350_PFMET100_v*") )
	{
	  if( PassTrigger[it] ) foundTrigger = true;
	}
      }

      //HLT_PFHT300_PFMET100_v
      //if( !foundTrigger ) std::cout << "FL: trigger not found" << std::endl;
      //FSLICHEP//if( !foundTrigger ) continue;

      //std::cout << "FSL: test" << std::endl;

      bool passLeptVeto = trCS.getVar<bool>("passLeptVeto"+spec);
      bool passnJets = trCS.getVar<bool>("passnJets"+spec);
      bool passMET = trCS.getVar<bool>("passMET"+spec);
      bool passHT = trCS.getVar<bool>("passHT"+spec);
      bool passMT2 = trCS.getVar<bool>("passMT2"+spec);
      bool passTagger = trCS.getVar<bool>("passTagger"+spec);
      bool passBJets = trCS.getVar<bool>("passBJets"+spec);
      bool passNoiseEventFilter = trCS.getVar<bool>("passNoiseEventFilter"+spec);
      bool passdPhis = trCS.getVar<bool>("passdPhis"+spec);

      //trCS.registerFunction(&mypdf);
      //const double pdfCentral = trCS.getVar<double>("NNPDF_From_Median_Central");
      //const double pdfUp = trCS.getVar<double>("NNPDF_From_Median_Up");
      //const double pdfDown = trCS.getVar<double>("NNPDF_From_Median_Down");
      //const double scaleUp = trCS.getVar<double>("Scaled_Variations_Up");
      //const double scaleDown = trCS.getVar<double>("Scaled_Variations_Down");
      //int iSR = theSearchBins.find_Binning_Index(cnt1CSVS, nTopCandSortedCnt_acc, MT2_acc, combmet);

      //std::cout << "pdfCentral = " << pdfCentral << ", pdfUp = " << pdfUp << std::endl;
      //std::cout << "scaleUp = " << scaleUp << std::endl;

      bool passInvertedBaseline = 
	//passLeptVeto &&
                          passnJets
                          && passHT
                          && passMT2
                          && passTagger
                          && passBJets
                          && passNoiseEventFilter
	                  && passMET
	                  && !passdPhis;


      if(
         passBaselinelostlept
         //passInvertedBaseline
        )
      {
        int nElectrons = trCS.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = trCS.getVar<int>("nMuons_CUT"+spec);

        double met = trCS.getVar<double>("met");
        double metphi = trCS.getVar<double>("metphi");

        int njets30 = trCS.getVar<int>("cntNJetsPt30Eta24"+spec);
        int njets50 = trCS.getVar<int>("cntNJetsPt50Eta24"+spec);
        int ntopjets = trCS.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = trCS.getVar<int>("cntCSVS"+spec);
        double MT2 = trCS.getVar<double>("best_had_brJet_MT2"+spec);
        //double ht = trCS.getVar<double>("ht");
        double ht = trCS.getVar<double>("HT"+spec);
        double mht = trCS.getVar<double>("mht");

	//const std::vector<double> &metMagUp = trCS.getVec<double>("metMagUp");
	//const std::vector<double> &metMagDown = trCS.getVec<double>("metMagDown");

	//double metjecUp = metMagUp[1];//met Uncertainty
	//double metjecLow = metMagDown[1];//met Uncertainty
	//double metjerUp = metMagUp[0];//met Uncertainty
	//double metjerLow = metMagDown[0];//met Uncertainty

	//std::cout << "met = " << met << " , metjecUp = " << metjecUp << " , metjecLow = " << metjecLow << " , metjerUp = " << metjerUp << " , metjerLow = " << metjerLow << std::endl;

	if (met<200) std::cout << "met<200" << std::endl;

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
          //double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * metjecUp * ( 1.0 - cos(deltaphi_mus) ) );
          //double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * metjerUp * ( 1.0 - cos(deltaphi_mus) ) );
          //double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * metjecLow * ( 1.0 - cos(deltaphi_mus) ) );
          //double mtW_mus = std::sqrt( 2.0 * reco_mus_pt * metjerLow * ( 1.0 - cos(deltaphi_mus) ) );

	  //if (met>450.0)
	  //{
	  //  std::cout << "met = " << met << std::endl;
	  //  std::cout << "mtW_mus = " << mtW_mus << std::endl;
	  //  std::cout << "muon pt = " << reco_mus_pt << std::endl;
	  //}

          if ( mtW_mus < 100.0 )
          {
	    //////////////////////////
	    // prediction computation
	    //////////////////////////
	    double EventWeight_mus = thisweight;
            const int njetsbin_number = Set_njetsbin_number(njets30);
            const int ptbin_number = Set_ptbin_number(reco_mus_pt);
            const int acbin_number = Set_acbin_number(activity);
	    const int htbin_number = Set_HTbin_number(ht);
            //const int htbin_number = Set_MT2bin_number(MT2);
            int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );

	    //// for signal scan
	    ////std::cout << "SusyMotherMass = " << SusyMotherMass << std::endl;
	    ////std::cout << "SusyLSPMass = " << SusyLSPMass << std::endl;
	    //std::pair<int, int> iMP((int)SusyMotherMass, (int)SusyLSPMass);
	    //std::vector<std::pair<int, int> >::iterator it;
	    //it = find (viMP.begin(), viMP.end(), iMP);
	    //if (it == viMP.end()) viMP.push_back(iMP);

	    //if (met<200.0) std::cout << "searchbin_id = " << searchbin_id << std::endl;

	    myAccRecoIsoEffs.EffstoWeights_fromH();

	    ++myAccRecoIsoEffs.nevents_cs_mus_sb[searchbin_id];

	    //std::cout << "ttbar_mtwcorrfactor[0] = " << ttbar_mtwcorrfactor[0] << std::endl;
	    //mtwcorrfactor
	    EventWeight_mus = EventWeight_mus * ttbar_mtwcorrfactor[ptbin_number];
	    //dimuon correction factor
            EventWeight_mus = EventWeight_mus * ttbar_corrfactor_di_mus;

	    //if (applyisotrkveto) EventWeight_mus *= isotrackvetoeff;

            //muon prediction from muon CS
	    //Fill muon iso closure plots
	    (myClosureHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    //Fill muon id closure plots
	    (myClosureHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus);
	    //Fill muon acc closure plots
	    (myClosureHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus);
	    //Fill all muon closure plots
	    double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number] + myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number] + myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id];

	    (myClosureHistgram.h_pred_mu_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

	    (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_2d_met_mupt)->Fill(deltaphi_mus, reco_mus_pt, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_lept_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);
	  
  	    (myClosureHistgram.h_exp_lept_all_Z_sb)->Fill(searchbin_id,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_njets30)->Fill(njets30,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_njets50)->Fill(njets50,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_ntops)->Fill(ntopjets,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_nbjets)->Fill(nbottomjets,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_MET)->Fill(met,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_MT2)->Fill(MT2,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_HT)->Fill(ht,EventWeight_all_mus*EventWeight_mus*isoTrackEff_SB[searchbin_id]);


            //total events flow for muons, prediction
            myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;

            if( searchbin_id >= 0 )
            {
              myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[searchbin_id] += EventWeight_all_mus*EventWeight_mus;

	      //// signal scan
	      //std::array<double,NSEARCH_BINS> vspred;
	      //std::array<double,NSEARCH_BINS> vsuncprec; 
	      //std::pair<std::array<double,NSEARCH_BINS>, std::array<double,NSEARCH_BINS> > signalPredAndUnc(vspred,vsuncprec);
	      //auto iter = mapSignalpred.find(iMP);
	      //if(iter == mapSignalpred.end()) iter = mapSignalpred.emplace(iMP, signalPredAndUnc).first;
	      //((iter->second).first)[searchbin_id]+=EventWeight_all_mus*EventWeight_mus;
	      //((iter->second).second)[searchbin_id]+=EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;

              myAccRecoIsoEffs.nevents_mus_pred_iso_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_reco_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_acc_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus;

              //myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number]*EventWeight_mus;

	    if (met>175.0 && met<=200.0)
	    {
	      //std::cout << "met = " << met << std::endl;
	      //std::cout << "MT2 = " << MT2 << std::endl;
	      if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_mus*EventWeight_mus;
	      if (MT2>300.0) sidebandHighMT2 += EventWeight_all_mus*EventWeight_mus;
	    }

            //begin to predict lost electrons from muon CS
            double EventWeight_els = thisweight;
            //mtwcorrfactor
            EventWeight_els = EventWeight_els * ttbar_mtwcorrfactor[ptbin_number];
            //dielectron correction factor
            EventWeight_els = EventWeight_els * ttbar_corrfactor_di_els;
  
	    //if (applyisotrkveto) EventWeight_els *= isotrackvetoeff;

            //electron prediction from muon CS
            //Fill electron iso closure plots
            (myClosureHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            //Fill electron id closure plots
            (myClosureHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            //Fill electron acc closure plots
            (myClosureHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els);
            //Fill all electron closure plots
            double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id];

            (myClosureHistgram.h_pred_el_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);

            (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_2d_met_mupt)->Fill(deltaphi_mus, reco_mus_pt, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);

  	    (myClosureHistgram.h_exp_lept_all_Z_sb)->Fill(searchbin_id,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_njets30)->Fill(njets30,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_njets50)->Fill(njets50,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_ntops)->Fill(ntopjets,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_nbjets)->Fill(nbottomjets,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_MET)->Fill(met,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_MT2)->Fill(MT2,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);
  	    (myClosureHistgram.h_exp_lept_all_Z_HT)->Fill(ht,EventWeight_all_els*EventWeight_els*isoTrackEff_SB[searchbin_id]);

            //total events flow for electrons, prediction
            myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;

	    if (met>175.0 && met<=200.0)
	    {
	      if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_els*EventWeight_els;
	      if (MT2>300.0) sidebandHighMT2 += EventWeight_all_els*EventWeight_els;
	    }


            if( searchbin_id >= 0 )
            {
              myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[searchbin_id] += EventWeight_all_els*EventWeight_els;
              //myAccRecoIsoEffs.nevents_els_pred_SB_MC[searchbin_id] += EventWeight_all_els*EventWeight_els/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_els_pred_SB_MC[searchbin_id] += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;

	      // signal scan
	      //auto iter = mapSignalpred.find(iMP);
	      //if(iter == mapSignalpred.end()) std::cout << "Error: iter == mapSignalpred.end()" << std::endl;
	      //else
	      //{
	      //	((iter->second).first)[searchbin_id]+=EventWeight_all_els*EventWeight_els;
	      //	((iter->second).second)[searchbin_id]+=EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
	      //}

            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_els_err += EventWeight_all_els*EventWeight_els*EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][htbin_number][searchbin_id]*EventWeight_els;
          }//mtW5_els<125 (muon CS)
        }//nElectrons == 0 && nMuons == 1 (muon CS)
      }//baseline_nolepveto
    }//TTJets samples class
  }//end of second loop

  //std::cout << "sidebandLowMT2 = " << sidebandLowMT2 << std::endl;
  //std::cout << "sidebandHighMT2 = " << sidebandHighMT2 << std::endl;

  //double aveTFfromMC[37]={0.50938294641, 0.54479857283, 0.55589368664, 0.47562009798, 0.47313954851, 0.55585283862, 0.63137075827, 0.49694571031, 0.45654695127, 0.73899860882, 0.77450214989, 0.51858315555, 0.58676123185, 0.61502503663, 0.65623508004, 0.49805541749, 0.62407702117, 0.6823095162, 0.61661615032, 0.60866461656, 0.91230657595, 0.30047597469, 0.31295993515, 0.26905854413, 0.45670624417, 0.47015903259, 0.35177654579, 0.50399066718, 0.44055971649, 0.32361672859, 0.3155307976, 0.30018569397, 0.45286988726, 0.50602656801, 0.41621381186, 0.44509945257, 0.61782010353};


  //double aveTFfromMC[45]={0.50938294641, 0.54479857283, 0.55589368664, 0.47562009798, 0.47313954851, 0.55585283862, 0.63137075827, 0.49694571031, 0.45654695127, 0.73899860882, 0.77450214989, 0.52104365541, 0.60239121649, 0.61841946255, 0.68602977735, 0.47571369344, 0.61390093265, 0.70282326462, 0.59551049235, 0.62121857692, 0.79823741008, 0.52378337381, 0.54342635586, 0.70728870479, 0.29863341115, 0.30230795586, 0.25385280375, 0.45733333925, 0.46685659786, 0.35692544256, 0.50504359785, 0.43853793836, 0.32081313052, 0.30348002588, 0.28963527215, 0.44438981801, 0.51552256953, 0.40326995862, 0.36276343558, 0.64257071159, 0.39887754975, 0.38596378496, 0.42245089568, 0.38479498238, 0.37833687369};

  //double aveTFfromMC[59]={0.55347, 0.67387, 0.61141, 0.54937, 0.54725, 0.63727, 0.55182, 0.62907, 0.30486, 0.58223, 0.57274, 0.71943, 0.57290, 0.72634, 0.72020, 0.33790, 0.53360, 0.72222, 1.21865, 1.02533, 0.74698, 0.71353, 0.66977, 0.52476, 0.50922, 0.89214, 0.74934, 0.68842, 0.36153, 0.34438, 0.52775, 0.56833, 0.60254, 0.60309, 0.76332, 0.46730, 0.46887, 0.45869, 0.62370, 0.54502, 0.34302, 0.37405, 0.40613, 0.58999, 0.55086, 0.59882, 0.39324, 0.90137, 0.90115, 0.37601, 0.24105, 0.57766, 0.34594, 0.42217, 0.33566, 0.37051, 0.19585, 0.36191, 0.49791};

  double aveTFfromMC[59]={0.521, 0.689, 0.644, 0.717, 0.504, 0.613, 0.561, 0.871, 0.238, 0.574, 0.505, 0.715, 0.533, 0.652, 0.697, 0.259, 0.467, 0.469, 0.851, 0.779, 0.615, 0.878, 0.335, 0.532, 0.624, 0.737, 0.703, 0.838, 0.352, 0.334, 0.373, 0.246, 0.582, 0.573, 0.758, 0.523, 0.476, 0.452, 0.464, 0.900, 0.326, 0.318, 0.431, 0.541, 0.498, 0.558, 0.460, 0.692, 0.730, 0.387, 0.269, 0.596, 0.493, 0.366, 0.358, 0.374, 0.248, 0.439, 0.178};

  //double psystup[37]={0.149294 ,  0.156028 ,  0.170782 ,  0.251557 ,  0.157696 ,  0.159797 ,  0.196816 ,  0.269567 ,  0.217605 ,  0.188411 ,  0.205831 ,  0.152414 ,  0.1556 ,  0.184521 ,  0.271894 ,  0.184819 ,  0.173421 ,  0.191517 ,  0.278766 ,  0.245826 ,  0.296231 ,  0.176127 ,  0.235712 ,  0.228626 ,  0.16598 ,  0.177104 ,  0.222109 ,  0.251645 ,  0.208664 ,  0.210948 ,  0.279213 ,  0.231705 ,  0.175164 ,  0.186562 ,  0.28631 ,  0.260183 ,  0.239561};

  //double psystup[45]={0.149706 ,  0.156635 ,  0.170685 ,  0.251444 ,  0.157666 ,  0.160303 ,  0.197513 ,  0.268814 ,  0.217869 ,  0.188927 ,  0.206648 ,  0.149137 ,  0.146228 ,  0.161018 ,  0.300102 ,  0.192106 ,  0.174856 ,  0.191143 ,  0.265615 ,  0.273642 ,  0.333063 ,  0.160709 ,  0.186936 ,  0.281573 ,  0.17605 ,  0.233689 ,  0.221836 ,  0.159932 ,  0.169461 ,  0.220388 ,  0.246237 ,  0.200758 ,  0.211591 ,  0.27611 ,  0.232983 ,  0.162081 ,  0.173356 ,  0.31402 ,  0.281801 ,  0.235023 ,  0.172506 ,  0.211313 ,  0.210235 ,  0.21777 ,  0.430502};

// double psystup[59]={0.153826 ,  0.165889 ,  0.247046 ,  1.66503 ,  0.178513 ,  0.177058 ,  0.318357 ,  0.554478 ,  0.332809 ,  0.243916 ,  0.249555 ,  0.661046 ,  0.15448 ,  0.172 ,  0.428146 ,  1.32462 ,  0.209322 ,  0.360871 ,  0.371814 ,  1.03167 ,  0.257624 ,  0.355837 ,  0.763525 ,  0.172662 ,  0.434775 ,  0.677375 ,  0.273899 ,  0.470478 ,  0.155295 ,  0.205738 ,  0.372261 ,  1.64157 ,  0.174096 ,  0.212279 ,  0.705538 ,  0.555891 ,  0.347277 ,  0.29998 ,  0.377704 ,  0.873694 ,  0.158215 ,  0.211727 ,  0.355306 ,  0.186035 ,  0.294232 ,  0.363652 ,  0.696344 ,  0.430562 ,  0.47057 ,  0.189598 ,  0.388406 ,  0.322059 ,  0.437904 ,  0.207241 ,  0.359394 ,  0.246785 ,  0.386214 ,  0.545055 ,  0.605631};

 double psystup[59]={0.153892 ,  0.164843 ,  0.222729 ,  0.407309 ,  0.173378 ,  0.177311 ,  0.265726 ,  0.328192 ,  0.387449 ,  0.210767 ,  0.247455 ,  0.261955 ,  0.153933 ,  0.166551 ,  0.285281 ,  1.83659 ,  0.183942 ,  0.20439 ,  0.286498 ,  0.486872 ,  0.222592 ,  0.338231 ,  0.422755 ,  0.162146 ,  0.210147 ,  0.368106 ,  0.207237 ,  0.24439 ,  0.157696 ,  0.232213 ,  0.408331 ,  3.73454 ,  0.175242 ,  0.215195 ,  0.425469 ,  0.694079 ,  0.344309 ,  0.253567 ,  0.442248 ,  0.449888 ,  0.156293 ,  0.184127 ,  0.289806 ,  0.177576 ,  0.260624 ,  0.334429 ,  0.511034 ,  0.253931 ,  0.336598 ,  0.176379 ,  0.291908 ,  0.231118 ,  0.312253 ,  0.224609 ,  0.422877 ,  0.219063 ,  0.366674 ,  0.367327 ,  0.616121};

 //double psystdown[37]={0.146536 ,  0.154 ,  0.168976 ,  0.249199 ,  0.153894 ,  0.157871 ,  0.195473 ,  0.26838 ,  0.213521 ,  0.187269 ,  0.204796 ,  0.15015 ,  0.153952 ,  0.183137 ,  0.270346 ,  0.182874 ,  0.172492 ,  0.190701 ,  0.278411 ,  0.244047 ,  0.295598 ,  0.169856 ,  0.231512 ,  0.219935 ,  0.162314 ,  0.173871 ,  0.218384 ,  0.249355 ,  0.206481 ,  0.205355 ,  0.274888 ,  0.223034 ,  0.170817 ,  0.182979 ,  0.282588 ,  0.257035 ,  0.236604 };

 // double psystdown[45]={0.146472 ,  0.154139 ,  0.168741 ,  0.249032 ,  0.153711 ,  0.157852 ,  0.195798 ,  0.267831 ,  0.213881 ,  0.187428 ,  0.205209 ,  0.146121 ,  0.143684 ,  0.15986 ,  0.299084 ,  0.189316 ,  0.173324 ,  0.190605 ,  0.265193 ,  0.271275 ,  0.332337 ,  0.158828 ,  0.18535 ,  0.280073 ,  0.169179 ,  0.228928 ,  0.213637 ,  0.156528 ,  0.166648 ,  0.217643 ,  0.244577 ,  0.19968 ,  0.206139 ,  0.272343 ,  0.223643 ,  0.159098 ,  0.170654 ,  0.311411 ,  0.278421 ,  0.233373 ,  0.169348 ,  0.207698 ,  0.205487 ,  0.211313 ,  0.427667};

// double psystdown[59]={0.155121 ,  0.167091 ,  0.247854 ,  1.66515 ,  0.17963 ,  0.178184 ,  0.318984 ,  0.554839 ,  0.33341 ,  0.244734 ,  0.250356 ,  0.661349 ,  0.155769 ,  0.173159 ,  0.428613 ,  1.32477 ,  0.210275 ,  0.361425 ,  0.372352 ,  1.03187 ,  0.258399 ,  0.356399 ,  0.763787 ,  0.173817 ,  0.435234 ,  0.67767 ,  0.274628 ,  0.470903 ,  0.156577 ,  0.206708 ,  0.372798 ,  1.64169 ,  0.175241 ,  0.213219 ,  0.705822 ,  0.556251 ,  0.347852 ,  0.300646 ,  0.378233 ,  0.873923 ,  0.159474 ,  0.21267 ,  0.355868 ,  0.187107 ,  0.294911 ,  0.364202 ,  0.696631 ,  0.431027 ,  0.470995 ,  0.19065 ,  0.38892 ,  0.322679 ,  0.438361 ,  0.208204 ,  0.35995 ,  0.247594 ,  0.386731 ,  0.545422 ,  0.605961};

 double psystdown[59]={0.155186 ,  0.166052 ,  0.223625 ,  0.4078 ,  0.174528 ,  0.178435 ,  0.266478 ,  0.328801 ,  0.387965 ,  0.211713 ,  0.248262 ,  0.262718 ,  0.155227 ,  0.167747 ,  0.285981 ,  1.8367 ,  0.185026 ,  0.205367 ,  0.287195 ,  0.487282 ,  0.223489 ,  0.338821 ,  0.423228 ,  0.163374 ,  0.211096 ,  0.368649 ,  0.2082 ,  0.245207 ,  0.158959 ,  0.233073 ,  0.40882 ,  3.7346 ,  0.176379 ,  0.216122 ,  0.425939 ,  0.694367 ,  0.34489 ,  0.254355 ,  0.4427 ,  0.450332 ,  0.157567 ,  0.18521 ,  0.290495 ,  0.178699 ,  0.26139 ,  0.335027 ,  0.511425 ,  0.254717 ,  0.337191 ,  0.177509 ,  0.292592 ,  0.231981 ,  0.312893 ,  0.225498 ,  0.42335 ,  0.219974 ,  0.367219 ,  0.367871 ,  0.616445};

 std::cout.precision(1);
 //std::cout.precision(3);
 //std::cout.precision(5);
 std::cout << std::fixed;
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
    myClosureHistgram.h_pred_lept_sb_isotrk->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] );
    myClosureHistgram.h_pred_lept_sb_isotrk->SetBinError( i_cal+1 , (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]);

    myClosureHistgram.h_pred_mu_iso_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_iso_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_mu_iso_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_iso_SB_MC[i_cal]));
    myClosureHistgram.h_pred_mu_reco_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_reco_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_mu_reco_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_reco_SB_MC[i_cal]));
    myClosureHistgram.h_pred_mu_acc_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_acc_SB_Normalized[i_cal] );
    //myClosureHistgram.h_pred_mu_acc_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_acc_SB_MC[i_cal]));

    //std::cout << "N events [" << i_cal << "] = " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal] << " , pred = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << std::endl;
 
    // data card
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal];
    //std::cout << " " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(60144642+59816364+30498962)/2262.946/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/4004.345/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/8000.0/831.76;
    // 1.83333 c'est 1/(0.43930872+0.10614564)
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
    //std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal];
    //if (myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]==0.0) std::cout << " 0";
    //else std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))/myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal];
    //if (myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]==0.0) std::cout << " 0";
    //else std::cout << " " << std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]);

    //std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal];
       //                //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ",";

    // syst
    //std::cout << "Pred[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredSystUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredEUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredESystUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;

    //std::cout << "uncUpPred[" << i_cal << "] = " << std::sqrt(((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])*((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])+1.84*aveTFfromMC[i_cal]*1.84*aveTFfromMC[i_cal]) << ";" << std::endl;
    //std::cout << "uncDownPred[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal] << ";" << std::endl;

    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(60144642+59816364+30498962)/2262.946/831.76 << ";" << std::endl;
    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(59654914+51873969+30587326)/2262.946/831.76 << ";" << std::endl;
    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/8000.0/831.76 << ";" << std::endl;




    if (true)
    {
    if (i_cal == 0) std::cout  <<  "0 & 1 & 1 & 200-350 & 200-350 & ";
    if (i_cal == 1) std::cout  <<  "1 & 1 & 1 & 200-350 & 350-500 & ";
    if (i_cal == 2) std::cout  <<  "2 & 1 & 1 & 200-350 & 500-650 & ";
    if (i_cal == 3) std::cout  <<  "3 & 1 & 1 & 200-350 & 650+ & ";
    if (i_cal == 4) std::cout  <<  "4 & 1 & 1 & 350-450 & 200-350 & ";
    if (i_cal == 5) std::cout  <<  "5 & 1 & 1 & 350-450 & 350-500 & ";
    if (i_cal == 6) std::cout  <<  "6 & 1 & 1 & 350-450 & 500-650 & ";
    if (i_cal == 7) std::cout  <<  "7 & 1 & 1 & 350-450 & 650+ & ";
    if (i_cal == 8) std::cout  <<  "8 & 1 & 1 & 450+ & 200-350 & ";
    if (i_cal == 9) std::cout  <<  "9 & 1 & 1 & 450+ & 350-500 & ";
    if (i_cal == 10) std::cout << "10 & 1 & 1 & 450+ & 500-650 & ";
    if (i_cal == 11) std::cout << "11 & 1 & 1 & 450+ & 650+ & ";

    if (i_cal == 12) std::cout << "12 & 1 & 2 & 200-350 & 200-350 & ";
    if (i_cal == 13) std::cout << "13 & 1 & 2 & 200-350 & 350-500 & ";
    if (i_cal == 14) std::cout << "14 & 1 & 2 & 200-350 & 500-650 & ";
    if (i_cal == 15) std::cout << "15 & 1 & 2 & 200-350 & 650+ & ";
    if (i_cal == 16) std::cout << "16 & 1 & 2 & 350-450 & 200-350 & ";
    if (i_cal == 17) std::cout << "17 & 1 & 2 & 350-450 & 350-500 & ";
    if (i_cal == 18) std::cout << "18 & 1 & 2 & 350-450 & 500-650 & ";
    if (i_cal == 19) std::cout << "19 & 1 & 2 & 350-450 & 650+ & ";
    if (i_cal == 20) std::cout << "20 & 1 & 2 & 450+ & 200-500 & ";
    if (i_cal == 21) std::cout << "21 & 1 & 2 & 450+ & 500-650 & ";
    if (i_cal == 22) std::cout << "22 & 1 & 2 & 450+ & 650+ & ";

    if (i_cal == 23) std::cout << "23 & 1 & 3+ & 200-350 & 200-350 & ";
    if (i_cal == 24) std::cout << "24 & 1 & 3+ & 200-350 & 350-500 & ";
    if (i_cal == 25) std::cout << "25 & 1 & 3+ & 200-350 & 500+ & ";
    if (i_cal == 26) std::cout << "26 & 1 & 3+ & 350+ & 200-350 & ";
    if (i_cal == 27) std::cout << "27 & 1 & 3+ & 350+ & 350+ & ";

    if (i_cal == 28) std::cout << "28 & 2 & 1 & 200-350 & 200-350 & ";
    if (i_cal == 29) std::cout << "29 & 2 & 1 & 200-350 & 350-500 & ";
    if (i_cal == 30) std::cout << "30 & 2 & 1 & 200-350 & 500-650 & ";
    if (i_cal == 31) std::cout << "31 & 2 & 1 & 200-350 & 650+ & ";
    if (i_cal == 32) std::cout << "32 & 2 & 1 & 350-450 & 200-350 & ";
    if (i_cal == 33) std::cout << "33 & 2 & 1 & 350-450 & 350-500 & ";
    if (i_cal == 34) std::cout << "34 & 2 & 1 & 350-450 & 500-650 & ";
    if (i_cal == 35) std::cout << "35 & 2 & 1 & 350-450 & 650+ & ";
    if (i_cal == 36) std::cout << "36 & 2 & 1 & 450+ & 200-350 & ";
    if (i_cal == 37) std::cout << "37 & 2 & 1 & 450+ & 350-500 & ";
    if (i_cal == 38) std::cout << "38 & 2 & 1 & 450+ & 500-650 & ";
    if (i_cal == 39) std::cout << "39 & 2 & 1 & 450+ & 650+ & ";

    if (i_cal == 40) std::cout << "40 & 2 & 2 & 200-350 & 200-350 & ";
    if (i_cal == 41) std::cout << "41 & 2 & 2 & 200-350 & 350-500 & ";
    if (i_cal == 42) std::cout << "42 & 2 & 2 & 200-350 & 500+ & ";
    if (i_cal == 43) std::cout << "43 & 2 & 2 & 350-450 & 200-350 & ";
    if (i_cal == 44) std::cout << "44 & 2 & 2 & 350-450 & 350-500 & ";
    if (i_cal == 45) std::cout << "45 & 2 & 2 & 350-450 & 500+ & ";
    if (i_cal == 46) std::cout << "46 & 2 & 2 & 450+ & 200-350 & ";
    if (i_cal == 47) std::cout << "47 & 2 & 2 & 450+ & 350-500 & ";
    if (i_cal == 48) std::cout << "48 & 2 & 2 & 450+ & 500+ & ";

    if (i_cal == 49) std::cout << "49 & 2 & 3+ & 200-350 & 200-350 & ";
    if (i_cal == 50) std::cout << "50 & 2 & 3+ & 200-350 & 350+ & ";
    if (i_cal == 51) std::cout << "51 & 2 & 3+ & 350+ & 200-350 & ";
    if (i_cal == 52) std::cout << "52 & 2 & 3+ & 350+ & 350+ & ";

    if (i_cal == 53) std::cout << "53 & 3+ & 1 & 200+ & 200-350 & ";
    if (i_cal == 54) std::cout << "54 & 3+ & 1 & 200+ & 350+ & ";

    if (i_cal == 55) std::cout << "55 & 3+ & 2 & 200+ & 200-350 & ";
    if (i_cal == 56) std::cout << "56 & 3+ & 2 & 200+ & 350+ & ";

    if (i_cal == 57) std::cout << "57 & 3+ & 3+ & 200+ & 200-350 & ";
    if (i_cal == 58) std::cout << "58 & 3+ & 3+ & 200+ & 350+ & ";

//    if (i_cal == 11) std::cout << "11 & 1 & 2+ & 200-300 & 200-275 & ";
//    if (i_cal == 12) std::cout << "12 & 1 & 2+ & 200-300 & 275-350 & ";
//    if (i_cal == 13) std::cout << "13 & 1 & 2+ & 200-300 & 350-450 & ";
//    if (i_cal == 14) std::cout << "14 & 1 & 2+ & 200-300 & 450+ & ";
//    if (i_cal == 15) std::cout << "15 & 1 & 2+ & 300-400 & 200-275 & ";
//    if (i_cal == 16) std::cout << "16 & 1 & 2+ & 300-400 & 275-350 & ";
//    if (i_cal == 17) std::cout << "17 & 1 & 2+ & 300-400 & 350-450 & ";
//    if (i_cal == 18) std::cout << "18 & 1 & 2+ & 300-400 & 450+ & ";
//    if (i_cal == 19) std::cout << "19 & 1 & 2+ & 400+ & 200-450 & ";
//    if (i_cal == 20) std::cout << "20 & 1 & 2+ & 400+ & 450+ & ";
//    if (i_cal == 21) std::cout << "21 & 2+ & 1 & 200-300 & 200-275 & ";
//    if (i_cal == 22) std::cout << "22 & 2+ & 1 & 200-300 & 275-350 & ";
//    if (i_cal == 23) std::cout << "23 & 2+ & 1 & 200-300 & 350+ & ";
//    if (i_cal == 24) std::cout << "24 & 2+ & 1 & 300-400 & 200-275 & ";
//    if (i_cal == 25) std::cout << "25 & 2+ & 1 & 300-400 & 275-350 & ";
//    if (i_cal == 26) std::cout << "26 & 2+ & 1 & 300-400 & 350+ & ";
//    if (i_cal == 27) std::cout << "27 & 2+ & 1 & 400+ & 200-350 & ";
//    if (i_cal == 28) std::cout << "28 & 2+ & 1 & 400+ & 350+ & ";
//    if (i_cal == 29) std::cout << "29 & 2+ & 2+ & 200-300 & 200-275 & ";
//    if (i_cal == 30) std::cout << "30 & 2+ & 2+ & 200-300 & 275-350 & ";
//    if (i_cal == 31) std::cout << "31 & 2+ & 2+ & 200-300 & 350+ & ";
//    if (i_cal == 32) std::cout << "32 & 2+ & 2+ & 300-400 & 200-275 & ";
//    if (i_cal == 33) std::cout << "33 & 2+ & 2+ & 300-400 & 275-350 & ";
//    if (i_cal == 34) std::cout << "34 & 2+ & 2+ & 300-400 & 350+ & ";
//    if (i_cal == 35) std::cout << "35 & 2+ & 2+ & 400+ & 200-350 & ";
//    if (i_cal == 36) std::cout << "36 & 2+ & 2+ & 400+ & 350+ & ";

    std::cout << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << " $^{+" << std::sqrt(((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])*((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])+1.84*aveTFfromMC[i_cal]*1.84*aveTFfromMC[i_cal]) << "}_{-" << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal] << "}$ $^{+" << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]*psystup[i_cal] << "}_{-" << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]*psystdown[i_cal] << "}$" << " \\\\ " << 
std::endl;
    std::cout << "\\hline" << std::endl;
    }

  }
  std::cout  << std::endl;

//  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
//  {
//    std::cout << " " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
//  }
//  std::cout  << std::endl;
//
//  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
//  {
//    std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/4004.345/831.76;
//  }
//  std::cout  << std::endl;
//  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
//  {
//    std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
//  }
//  std::cout  << std::endl;



//  // signal scan
//  for (int signalPointc=0;signalPointc<viMP.size();++signalPointc)
//  {
//    char hname[128];
//    sprintf(hname, "%s_%d_%d", "h_signal_pred", viMP[signalPointc].first, viMP[signalPointc].second);
//
//    //std::cout << "viMP[signalPointc].first = " << viMP[signalPointc].first << " , viMP[signalPointc].second = " << viMP[signalPointc].second << std::endl;
//
//    TH1D *h_signal_pred_100_0;
//    h_signal_pred_100_0 = new TH1D(hname,hname,NSEARCH_BINS+1,0,NSEARCH_BINS+1);
//
//    auto iter = mapSignalpred.find(viMP[signalPointc]);
//    for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
//    {
//      h_signal_pred_100_0->SetBinContent( i_cal+1 , ((iter->second).first)[i_cal]*isoTrackEff_SB[i_cal] );
//      h_signal_pred_100_0->SetBinError( i_cal+1 , ((iter->second).second)[i_cal]*isoTrackEff_SB[i_cal]);
//    }
//  }
  std::cout << "pred is done" << std::endl;
  (myClosureHistgram.oFile)->Write();
  std::cout << "write is done" << std::endl;
  (myClosureHistgram.oFile)->Close();
  std::cout << "close is done" << std::endl;
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
  //const char *inputFileList_Exp_Pred = argv[2];

  //define my AccRecoIsoEffs class to stroe counts and efficiencies
  AccRecoIsoEffs myAccRecoIsoEffs;

  TTJetsSampleWeight myTTJetsSampleWeight;
  double W_Lept_BR = 0.1086*3;
  double TTbar_SingleLept_BR = 0.43930872; // 2*W_Lept_BR*(1-W_Lept_BR)
  double TTbar_DiLept_BR = 0.10614564; // W_Lept_BR^2
  // https://github.com/susy2015/SusyAnaTools/blob/master/Tools/samples.cc
  //TTJets nominal
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_", 831.76, 11339232, LUMI, inputFileList_Cal );
  //TTJets single lepton and di-lepton
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromT_", 831.76*0.5*TTbar_SingleLept_BR, 59654914, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromTbar", 831.76*0.5*TTbar_SingleLept_BR, 51873969, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_DiLept", 831.76*TTbar_DiLept_BR, 30587326, LUMI, inputFileList_Cal );

  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "stop", 0.0189612, 240685, LUMI, inputFileList_Cal ); // 850,100
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "HTMHT" , 1, 1, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "HTMHT" , 1, 1, 1.0, inputFileList_Cal );

  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "tW_top" , 35.6, 998400, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "tW_antitop" , 35.6, 967600, LUMI, inputFileList_Cal );

  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-400To600" , 48.91, 7299788, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-600To800" , 12.05, 3723054, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-800To1200" , 5.501, 7480017, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-1200To2500" , 1.329, 7042924, LUMI, inputFileList_Cal );
  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-2500ToInf" , 0.03216, 252809, LUMI, inputFileList_Cal );

  //myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "SMS-T2tt_mStop" , 1.0, 1.0, 1.0, inputFileList_Cal );

  //TTJetsSampleWeight myExpPredSampleWeight;
  //myExpPredSampleWeight.TTJetsSampleInfo_push_back( "TTJets_", 831.76, 11339232, LUMI, inputFileList_Exp_Pred );

  //LoopLLCal( myAccRecoIsoEffs, myTTJetsSampleWeight );
  //LoopLLExp( myAccRecoIsoEffs, myTTJetsSampleWeight );
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight );

  std::cout << "done" << std::endl;
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

    //els_acc[i_cal] = nels_acc[i_cal]/nels[i_cal];
    //els_acc_err[i_cal] = std::sqrt( get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal])*get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal]) + get_sys_Error(els_acc[i_cal],0.09)*get_sys_Error(els_acc[i_cal],0.09) );
    //els_acc_err[i_cal] = get_stat_Error(nels_acc_MC[i_cal],nels_MC[i_cal]);
  }


  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    //std::cout << "acc[" << searchbinc << "] = " << nmus_acc_sb[searchbinc]/nmus_sb[searchbinc] << ";" << std::endl;
    std::cout << "err_acc[" << searchbinc << "] = " << get_stat_Error(nmus_acc_MC_sb[searchbinc],nmus_MC_sb[searchbinc]) << ";" << std::endl;
    els_acc[searchbinc] = nels_acc[searchbinc]/nels[searchbinc];
    els_acc_err[searchbinc] = get_stat_Error(nels_acc_MC[searchbinc],nels_MC[searchbinc]);
  }


  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    std::cout << "err_acc_el[" << searchbinc << "] = " << els_acc_err[searchbinc] << ";" << std::endl;
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

  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
  for(i_cal = 0 ; i_cal < NJETS_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < PT_BINS ; j_cal++)
    {
      for(k_cal = 0 ; k_cal < AC_BINS ; k_cal++)
      {
	for (int htbinc=0;htbinc<NHT_BINS;++htbinc)
	{
	  //mus_EventWeight_iso[i_cal][j_cal][k_cal]  = (1.0 - mus_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal];
	  //mus_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1.0 - mus_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] ); 
	  //mus_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - mus_acc[i_cal])/mus_acc[i_cal] );
    
	  //els_EventWeight_acc[i_cal][j_cal][k_cal]  = (1.0/mus_isoeff[j_cal][k_cal]) * (1.0/mus_recoeff[j_cal][k_cal]) * ( (1.0 - els_acc[i_cal])/mus_acc[i_cal] );
	  //els_EventWeight_reco[i_cal][j_cal][k_cal] = (1.0/mus_isoeff[j_cal][k_cal]) * ( (1 - els_recoeff[j_cal][k_cal])/mus_recoeff[j_cal][k_cal] )* (els_acc[i_cal]/mus_acc[i_cal]); 
	  //els_EventWeight_iso[i_cal][j_cal][k_cal]  = ( (1.0 - els_isoeff[j_cal][k_cal])/mus_isoeff[j_cal][k_cal] ) * (els_recoeff[j_cal][k_cal]/mus_recoeff[j_cal][k_cal])* (els_acc[i_cal]/mus_acc[i_cal]);

	  mus_EventWeight_iso[i_cal][j_cal][k_cal][htbinc]  = (1.0 - ttbar_mus_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal];
	  mus_EventWeight_reco[i_cal][j_cal][k_cal][htbinc] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] );
	  //mus_EventWeight_acc[i_cal][j_cal][k_cal][htbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_acc[i_cal][htbinc])/ttbar_mus_acc[i_cal][htbinc] );
	  mus_EventWeight_acc[i_cal][j_cal][k_cal][htbinc][searchbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_acc[searchbinc])/ttbar_mus_acc[searchbinc] );

	  //els_EventWeight_acc[i_cal][j_cal][k_cal][htbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_els_acc[i_cal])/ttbar_mus_acc[i_cal][htbinc] );
	  //els_EventWeight_reco[i_cal][j_cal][k_cal][htbinc] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1 - ttbar_els_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] )* (ttbar_els_acc[i_cal]/ttbar_mus_acc[i_cal][htbinc]);
	  //els_EventWeight_iso[i_cal][j_cal][k_cal][htbinc]  = ( (1.0 - ttbar_els_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal] ) * (ttbar_els_recoeff[j_cal][k_cal]/ttbar_mus_recoeff[j_cal][k_cal])* (ttbar_els_acc[i_cal]/ttbar_mus_acc[i_cal][htbinc]);
	  els_EventWeight_acc[i_cal][j_cal][k_cal][htbinc][searchbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_els_acc[searchbinc])/ttbar_mus_acc[searchbinc] );
	  els_EventWeight_reco[i_cal][j_cal][k_cal][htbinc][searchbinc] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1 - ttbar_els_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] )* (ttbar_els_acc[searchbinc]/ttbar_mus_acc[searchbinc]);
	  els_EventWeight_iso[i_cal][j_cal][k_cal][htbinc][searchbinc]  = ( (1.0 - ttbar_els_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal] ) * (ttbar_els_recoeff[j_cal][k_cal]/ttbar_mus_recoeff[j_cal][k_cal])* (ttbar_els_acc[searchbinc]/ttbar_mus_acc[searchbinc]);
	}
      }
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

  EffsHeader << "  const double ttbar_mus_acc[" << NSEARCH_BINS << "] = "; 
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    if( searchbinc == 0 ) { EffsHeader << "{"; }
    EffsHeader << nmus_acc_sb[searchbinc]/nmus_sb[searchbinc];
    if( searchbinc != NSEARCH_BINS-1 ) { EffsHeader << ","; }
    if( searchbinc == NSEARCH_BINS-1 ) { EffsHeader << "};" << std::endl; }
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

  EffsHeader << "  const double ttbar_els_acc[" << NSEARCH_BINS << "] = ";                            
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    if( searchbinc == 0 ) { EffsHeader << "{"; }
    EffsHeader << els_acc[searchbinc];
    if( searchbinc != NSEARCH_BINS-1 ) { EffsHeader << ","; }
    if( searchbinc == NSEARCH_BINS-1 ) { EffsHeader << "};" << std::endl; }
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

  EffsHeader << "  double isoTrackEff_SB[" << NSEARCH_BINS << "] = {" << std::endl;
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    EffsHeader << isotrkeff[searchbinc] << ", ";
  }
  EffsHeader << "};" << std::endl;

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
