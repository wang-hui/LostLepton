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
#include "TGraphAsymmErrors.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "TROOT.h"
//#include "TInterpreter.h"

#include "Activity.h"
#include "LostLepton_MuCS_TTbar.h"
#include "TTJetsReWeighting.h"
//#include "v160714_newMuonID_accNoSingleTop_bin7f6_trkSF_EffsHeader_MuCS.h"
#include "v3_EffsHeader_MuCS_data.h"
#include "TriggerEff.h"
//#include "SusyAnaTools/Tools/PDFUncertainty.h"

//const double isotrackvetoeff = 0.563499421;
const bool applyisotrkveto = false; // should be false
//const double isotrackvetoeff = 1;
double isotrkeff[NSEARCH_BINS];

const bool use_muon_control_sample = false;
const bool use_electron_control_sample = true;
const double electron_purity = 0.96;
const bool isData = true;

void LoopLLCal( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight )
{
  BaseHistgram myBaseHistgram;
  myBaseHistgram.BookHistgram("CalLL.root");

  NTupleReader *tr =0;

  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(*tr, spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  SearchBins theSearchBins("SB_v1_2017");

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

	TH2D * mu_mediumID_SF = 0, * mu_miniISO_SF = 0;
	TGraphAsymmErrors * mu_trkptGT10_SF = 0, * mu_trkptLT10_SF = 0;

	TH2D * ele_VetoID_SF = 0, * ele_miniISO_SF = 0;
	TH2D * ele_trkpt_SF = 0;
		
	if( isData ){
    	// Lepton SF
	TFile * allINone_leptonSF_file = new TFile("allINone_leptonSF_Moriond17.root");
    	mu_mediumID_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0");
      	mu_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("pt_abseta_PLOT_pair_probeMultiplicity_bin0_&_Medium2016_pass");
	// SF for muon tracks
      	mu_trkptGT10_SF = (TGraphAsymmErrors*) allINone_leptonSF_file->Get("ratio_eff_eta3_dr030e030_corr");
      	mu_trkptLT10_SF = (TGraphAsymmErrors*) allINone_leptonSF_file->Get("ratio_eff_eta3_tk0_dr030e030_corr");
  
     	ele_VetoID_SF = (TH2D*) allINone_leptonSF_file->Get("GsfElectronToVeto");
      	ele_miniISO_SF = (TH2D*) allINone_leptonSF_file->Get("MVAVLooseElectronToMini");
      	ele_trkpt_SF = (TH2D*) allINone_leptonSF_file->Get("EGamma_SF2D");
      	}

  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  {  
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file

    NTupleReader tr((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    //AnaFunctions::prepareTopTagger();

    tr.registerFunction(&mypassBaselineFunc);    
 
    //PDFUncertainty pdf;
    //tr.registerFunction(pdf);

    double thisweight = (*iter_TTJetsSampleInfos).weight;

    std::cout << (*iter_TTJetsSampleInfos).TTJetsTag << " Weight: " << thisweight << std::endl;

    int neventc=0;

    //std::cout << "tr.getNextEvent()=" << tr.getNextEvent() <<std::endl;

    while(tr.getNextEvent())
    //while(tr.getNextEvent() && neventc<1000)
    {
      ++neventc;
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      const double genHT = tr.hasVar("genHT") ? tr.getVar<double>("genHT") : -999;
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromT_" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromTbar" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_DiLept") && genHT >= 600) continue; 
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-600to800" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-800to1200" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-1200to2500" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-2500toInf") && genHT < 600) continue;

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


      //if (passBaselinelostlept)
      if (true)		//for iso and reco only!!!
      {
	//const double pdf_Unc_Up = tr.getVar<double>("NNPDF_From_Median_Up");
	//const double pdf_Unc_Down = tr.getVar<double>("NNPDF_From_Median_Down");
        //const double scaleUp = tr.getVar<double>("Scaled_Variations_Up");
        //const double scaleDown = tr.getVar<double>("Scaled_Variations_Down");
 
        //std::cout << "pdfUp = " << pdf_Unc_Up << std::endl;
	const double systWeight=1.0;
	//double systWeight=pdf_Unc_Up;
	//const double systWeight=pdf_Unc_Down;
	//const double systWeight=scaleUp;
	//const double systWeight=scaleDown;

	//if (pdf_Unc_Up>2.0) systWeight=0.0;

        myAccRecoIsoEffs.nevents_sel_base+=thisweight;

        int nElectrons = tr.getVar<int>("nElectrons_CUT"+spec);
        int nMuons = tr.getVar<int>("nMuons_CUT"+spec);


        //searchbin variables
        int ntopjets = tr.getVar<int>("nTopCandSortedCnt"+spec);
        int nbottomjets = tr.getVar<int>("cntCSVS"+spec);
        double MT2 = tr.getVar<double>("best_had_brJet_MT2"+spec);
        double met = tr.getVar<double>("met");
        double metphi = tr.getVar<double>("metphi");
	//const double ht = tr.getVar<double>("ht");
	const double ht = tr.getVar<double>("HT"+spec);
	//std::cout << "ht = " << ht << std::endl;

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
              //const int njetsbin_number = Set_njetsbin_number(njets30);
	      //const int HTbin_number = Set_HTbin_number(ht);
	      //const int HTbin_number = Set_MT2bin_number(MT2);
	      //std::cout << "HTbin_number = " << HTbin_number << std::endl;            

              //myAccRecoIsoEffs.nmus[njetsbin_number][HTbin_number]+=thisweight;
              //myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
              //myAccRecoIsoEffs.nmus_MC[njetsbin_number][HTbin_number]++;

	      //const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
	      //59 bins
	      const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

	      myAccRecoIsoEffs.nmus_sb[searchbin_id]+=thisweight*systWeight;

	      myAccRecoIsoEffs.nmus_MC_sb[searchbin_id]++;

		double mu_id_SF = 1.0, mu_iso_SF = 1.0, mu_trk_SF = 1.0;
		double mu_id_SF_allreco = 1.0, mu_iso_SF_allreco = 1.0, mu_trk_SF_allreco = 1.0;	
	
		if (isData)
		{	
		mu_id_SF = mu_mediumID_SF->GetBinContent(mu_mediumID_SF->FindBin(myLostMuonObj.gen_pt, fabs(myLostMuonObj.gen_eta)));
		if( mu_id_SF == 0 ) mu_id_SF = 1.0; // very simple way dealing with out of range issue of the TH2D

      		mu_iso_SF = mu_miniISO_SF->GetBinContent(mu_miniISO_SF->FindBin(myLostMuonObj.gen_pt, fabs(myLostMuonObj.gen_eta)));
		if( mu_iso_SF == 0 ) mu_iso_SF = 1.0;

       		if( myLostMuonObj.gen_pt < 10 ){ mu_trk_SF = mu_trkptLT10_SF->Eval(myLostMuonObj.gen_eta); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0;}
    		if( myLostMuonObj.gen_pt >= 10 ){ mu_trk_SF = mu_trkptGT10_SF->Eval(myLostMuonObj.gen_eta); if( mu_trk_SF == 0 ) mu_trk_SF = 1.0;}
        	//mu_id_SF_err = mu_id_SF * 0.03;
		
		mu_id_SF_allreco = mu_mediumID_SF->GetBinContent(mu_mediumID_SF->FindBin(myLostMuonObj.reco_pt, fabs(myLostMuonObj.reco_eta)));
		if( mu_id_SF_allreco == 0 ) mu_id_SF_allreco = 1.0; // very simple way dealing with out of range issue of the TH2D

      		mu_iso_SF_allreco = mu_miniISO_SF->GetBinContent(mu_miniISO_SF->FindBin(myLostMuonObj.reco_pt, fabs(myLostMuonObj.reco_eta)));
		if( mu_iso_SF_allreco == 0 ) mu_iso_SF_allreco = 1.0;

       		if( myLostMuonObj.reco_pt < 10 )
		{mu_trk_SF_allreco = mu_trkptLT10_SF->Eval(myLostMuonObj.reco_eta);
		if( mu_trk_SF_allreco == 0 ) mu_trk_SF_allreco = 1.0;}

    		if( myLostMuonObj.reco_pt >= 10 )
		{mu_trk_SF_allreco = mu_trkptGT10_SF->Eval(myLostMuonObj.reco_eta);
		if( mu_trk_SF_allreco == 0 ) mu_trk_SF_allreco = 1.0;}
		}

              if( myLostMuonObj.passAcc )
              {
                //myAccRecoIsoEffs.nmus_acc[njetsbin_number][HTbin_number]+=thisweight;
                //myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]+=thisweight*thisweight;
                //myAccRecoIsoEffs.nmus_acc_MC[njetsbin_number][HTbin_number]++;

                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);
                int Etabin_number = Set_Etabin_number(fabs(myLostMuonObj.gen_eta));
		//std::cout << "myLostMuonObj.gen_eta = " << abs(myLostMuonObj.gen_eta) << " Etabin_number = " << Etabin_number << std::endl;

		myAccRecoIsoEffs.nmus_acc_bin[ptbin_number][Etabin_number]+=thisweight;
		myAccRecoIsoEffs.nmus_acc_bin_MC[ptbin_number][Etabin_number]++;

		myAccRecoIsoEffs.nmus_acc_sb[searchbin_id]+=thisweight*systWeight;
		myAccRecoIsoEffs.nmus_acc_MC_sb[searchbin_id]++;

		myAccRecoIsoEffs.nmus_acc_weight_squar[searchbin_id]+=thisweight*thisweight;
              }

	      else myAccRecoIsoEffs.nmus_not_acc_weight_squar[searchbin_id]+=thisweight*thisweight;

		// apply muon tracking SF:
		double TrkSF=1.0;
		/*if (myLostMuonObj.reco_eta<-2.1) TrkSF=0.982399;
		else if (myLostMuonObj.reco_eta<-1.6) TrkSF=0.9917468;
		else if (myLostMuonObj.reco_eta<-1.1) TrkSF=0.995945;
		else if (myLostMuonObj.reco_eta<-0.6) TrkSF=0.9934131;
		else if (myLostMuonObj.reco_eta<0.0) TrkSF=0.9914607;
		else if (myLostMuonObj.reco_eta<0.6) TrkSF=0.9946801;
		else if (myLostMuonObj.reco_eta<1.1) TrkSF=0.9966664;
		else if (myLostMuonObj.reco_eta<1.6) TrkSF=0.9949339;
		else if (myLostMuonObj.reco_eta<2.1) TrkSF=0.9911866;
		else TrkSF=0.9768119;
		*/
              //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt

	      if( myLostMuonObj.passId )
              {
                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);
                int Etabin_number = Set_Etabin_number(fabs(myLostMuonObj.gen_eta));

                myAccRecoIsoEffs.nmus_reco[ptbin_number][Etabin_number]+=thisweight*TrkSF*mu_trk_SF*mu_id_SF;
                myAccRecoIsoEffs.nmus_reco_MC[ptbin_number][Etabin_number]++;

                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostMuonObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostMuonObj.reco_activity);

                myAccRecoIsoEffs.nmus_reco_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight*TrkSF*mu_trk_SF_allreco*mu_id_SF_allreco;
                myAccRecoIsoEffs.nmus_reco_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //std::cout << myLostMuonObj.gen_activity <<"," << myLostMuonObj.reco_activity << std::endl;
                (myBaseHistgram.h_id_genactivity_mus)->Fill(myLostMuonObj.gen_activity,thisweight*TrkSF*mu_trk_SF*mu_id_SF);
                (myBaseHistgram.h_id_recoactivity_mus)->Fill(myLostMuonObj.reco_activity,thisweight*TrkSF*mu_trk_SF_allreco*mu_id_SF_allreco);
              }

              if( myLostMuonObj.passIso )
              {
                int ptbin_number = Set_ptbin_number(myLostMuonObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostMuonObj.gen_activity);

                myAccRecoIsoEffs.nmus_iso[ptbin_number][acbin_number]+=thisweight*TrkSF*mu_trk_SF*mu_id_SF*mu_iso_SF;
                myAccRecoIsoEffs.nmus_iso_MC[ptbin_number][acbin_number]++;
              
                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostMuonObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostMuonObj.reco_activity);
                
                myAccRecoIsoEffs.nmus_iso_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight*TrkSF*mu_trk_SF_allreco*mu_id_SF_allreco*mu_iso_SF_allreco;
                myAccRecoIsoEffs.nmus_iso_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
                //std::cout << muonsMiniIso.size() << "," << myLostMuonObj.reco_index << "ISO:" << muonsMiniIso.at(myLostMuonObj.reco_index) << std::endl;
              }

              //check warning function when we calculate the efficienies!
              if ( nMuons == 0 && myLostMuonObj.passIso )
              { 
		std::cout << "Warning: mu pass iso but nMuons=0!" << std::endl;
		std::cout << "reco_mus_pt = " << myLostMuonObj.reco_pt << std::endl;
		std::cout << "reco_mus_eta = " << myLostMuonObj.reco_eta << std::endl;

		std::cout << " isMu = " << myLostMuonObj.isMu << " isEl = " << myLostMuonObj.isEl << std::endl;
		std::cout << " gen_eta = " << myLostMuonObj.gen_eta << " gen_phi = " << myLostMuonObj.gen_phi << " gen_pt = " << myLostMuonObj.gen_pt << " gen_activity = " << myLostMuonObj.gen_activity << std::endl;
		std::cout << " reco_eta = " << myLostMuonObj.reco_eta << " reco_phi = " << myLostMuonObj.reco_phi << " reco_pt = " << myLostMuonObj.reco_pt << "reco_activity = " << myLostMuonObj.reco_activity << std::endl;
		std::cout << " passAcc = " << myLostMuonObj.passAcc << " passId = " << myLostMuonObj.passId << " passIso = " << myLostMuonObj.passIso << std::endl;
		std::cout << " doneAcc = " << myLostMuonObj.doneAcc << " doneId = " << myLostMuonObj.doneId << " doneIso = " << myLostMuonObj.doneIso << std::endl;
		std::cout << " reco_index = " << myLostMuonObj.reco_index << std::endl;

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

                //std::cout << "myLostElectronObj.isEl=" << myLostElectronObj.isEl << std::endl;

            if( myLostElectronObj.isEl )
            {
	      ngenel++;
              int njetsbin_number = Set_njetsbin_number(njets30);
	      //const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );   
	      //59 bins
	      const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);
            
              //myAccRecoIsoEffs.nels[njetsbin_number]+=thisweight;
              //myAccRecoIsoEffs.nels_MC[njetsbin_number]++;
	      myAccRecoIsoEffs.nels[searchbin_id]+=thisweight*systWeight;
	      myAccRecoIsoEffs.nels_MC[searchbin_id]++;

	      double ele_id_SF = 1.0, ele_iso_SF = 1.0, ele_trk_SF = 1.0;
      	      double ele_id_SF_err = 0.0, ele_iso_SF_err = 0.0, ele_trk_SF_err = 0.0;
	      
	      double ele_id_SF_allreco = 1.0, ele_iso_SF_allreco = 1.0, ele_trk_SF_allreco = 1.0;
      	      double ele_id_SF_err_allreco = 0.0, ele_iso_SF_err_allreco = 0.0, ele_trk_SF_err_allreco = 0.0;

	      if( isData )
      	      {
          	ele_id_SF = ele_VetoID_SF->GetBinContent(ele_VetoID_SF->FindBin(myLostElectronObj.gen_pt, fabs(myLostElectronObj.gen_eta)));
          	ele_id_SF_err = ele_VetoID_SF->GetBinError(ele_VetoID_SF->FindBin(myLostElectronObj.gen_pt, fabs(myLostElectronObj.gen_eta)));
          	if( ele_id_SF == 0 ){ ele_id_SF = 1.0; ele_id_SF_err = 0.0; }
        	// very simple way dealing with out of range issue of the TH2D
          	
		ele_iso_SF = ele_miniISO_SF->GetBinContent(ele_miniISO_SF->FindBin(myLostElectronObj.gen_pt, fabs(myLostElectronObj.gen_eta)));
          	ele_iso_SF_err = ele_miniISO_SF->GetBinError(ele_miniISO_SF->FindBin(myLostElectronObj.gen_pt, fabs(myLostElectronObj.gen_eta)));
          	if( ele_iso_SF == 0 ){ ele_iso_SF = 1.0; ele_iso_SF_err = 0.0; }
        	
        	
          	ele_trk_SF = ele_trkpt_SF->GetBinContent(ele_trkpt_SF->FindBin(myLostElectronObj.gen_eta, myLostElectronObj.gen_pt));
          	ele_trk_SF_err = myLostElectronObj.gen_pt<20 || myLostElectronObj.gen_pt>80? 0.01: 0.00;
          	if( ele_trk_SF == 0 ){ ele_trk_SF = 1.0; ele_trk_SF_err = 0.0; }
          	
		ele_id_SF_allreco = ele_VetoID_SF->GetBinContent(ele_VetoID_SF->FindBin(myLostElectronObj.reco_pt, fabs(myLostElectronObj.reco_eta)));
          	ele_id_SF_err_allreco = ele_VetoID_SF->GetBinError(ele_VetoID_SF->FindBin(myLostElectronObj.reco_pt, fabs(myLostElectronObj.reco_eta)));
          	if( ele_id_SF_allreco == 0 ){ ele_id_SF_allreco = 1.0; ele_id_SF_err_allreco = 0.0; }
        	// very simple way dealing with out of range issue of the TH2D
          	
		ele_iso_SF_allreco = ele_miniISO_SF->GetBinContent(ele_miniISO_SF->FindBin(myLostElectronObj.reco_pt, fabs(myLostElectronObj.reco_eta)));
          	ele_iso_SF_err_allreco = ele_miniISO_SF->GetBinError(ele_miniISO_SF->FindBin(myLostElectronObj.reco_pt, fabs(myLostElectronObj.reco_eta)));
          	if( ele_iso_SF_allreco == 0 ){ ele_iso_SF_allreco = 1.0; ele_iso_SF_err_allreco = 0.0; }
        	
        	
          	ele_trk_SF_allreco = ele_trkpt_SF->GetBinContent(ele_trkpt_SF->FindBin(myLostElectronObj.reco_eta, myLostElectronObj.reco_pt));
          	ele_trk_SF_err_allreco = myLostElectronObj.reco_pt<20 || myLostElectronObj.reco_pt>80? 0.01: 0.00;
          	if( ele_trk_SF_allreco == 0 ){ ele_trk_SF_allreco = 1.0; ele_trk_SF_err_allreco = 0.0; }
      	      }
              if( myLostElectronObj.passAcc )
              {
                //myAccRecoIsoEffs.nels_acc[njetsbin_number]+=thisweight;
                //myAccRecoIsoEffs.nels_acc_MC[njetsbin_number]++;
                myAccRecoIsoEffs.nels_acc[searchbin_id]+=thisweight*systWeight;
                myAccRecoIsoEffs.nels_acc_MC[searchbin_id]++;

		myAccRecoIsoEffs.nels_acc_weight_squar[searchbin_id]+=thisweight*thisweight;

                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);
                int Etabin_number = Set_Etabin_number(fabs(myLostElectronObj.gen_eta));

                myAccRecoIsoEffs.nels_acc_bin[ptbin_number][Etabin_number]+=thisweight;
                myAccRecoIsoEffs.nels_acc_bin_MC[ptbin_number][Etabin_number]++;
              }

	      else myAccRecoIsoEffs.nels_not_acc_weight_squar[searchbin_id]+=thisweight*thisweight;

		// apply electron tracking SF:
		// check with POG
		double TrkSF=1.0;
		/*if (myLostElectronObj.reco_eta<-2.4) TrkSF=1.170338;
		else if (myLostElectronObj.reco_eta<-2.3) TrkSF=1.00852;
		else if (myLostElectronObj.reco_eta<-2.2) TrkSF=1.010471;
		else if (myLostElectronObj.reco_eta<-2.0) TrkSF=1.005187;
		else if (myLostElectronObj.reco_eta<-1.8) TrkSF=0.9979317;
		else if (myLostElectronObj.reco_eta<-1.63) TrkSF=0.9917012;
		else if (myLostElectronObj.reco_eta<-1.566) TrkSF=0.9864865;
		else if (myLostElectronObj.reco_eta<-1.444) TrkSF=0.9615819;
		else if (myLostElectronObj.reco_eta<-1.2) TrkSF=0.9866667;
		else if (myLostElectronObj.reco_eta<-1.0) TrkSF=0.9775051;
		else if (myLostElectronObj.reco_eta<-0.6) TrkSF=0.9693878;
		else if (myLostElectronObj.reco_eta<-0.4) TrkSF=0.9663609;
		else if (myLostElectronObj.reco_eta<-0.2) TrkSF=0.9633027;
		else if (myLostElectronObj.reco_eta<0) TrkSF=0.96;
		else if (myLostElectronObj.reco_eta<0.2) TrkSF=0.9661885;
		else if (myLostElectronObj.reco_eta<0.4) TrkSF=0.9796334;
		else if (myLostElectronObj.reco_eta<0.6) TrkSF=0.9765784;
		else if (myLostElectronObj.reco_eta<1.0) TrkSF=0.9806517;
		else if (myLostElectronObj.reco_eta<1.2) TrkSF=0.9867347;
		else if (myLostElectronObj.reco_eta<1.444) TrkSF=0.9866803;
		else if (myLostElectronObj.reco_eta<1.566) TrkSF=0.9707207;
		else if (myLostElectronObj.reco_eta<1.63) TrkSF=0.9896694;
		else if (myLostElectronObj.reco_eta<1.8) TrkSF=0.995872;
		else if (myLostElectronObj.reco_eta<2.0) TrkSF=0.989733;
		else if (myLostElectronObj.reco_eta<2.2) TrkSF=0.9948612;
		else if (myLostElectronObj.reco_eta<2.3) TrkSF=0.9927686;
		else if (myLostElectronObj.reco_eta<2.4) TrkSF=0.9666319;
		else TrkSF=0.8840206;
		*/
              //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
              if( myLostElectronObj.passId )
              {

                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);
                int Etabin_number = Set_Etabin_number(fabs(myLostElectronObj.gen_eta));

                myAccRecoIsoEffs.nels_reco[ptbin_number][Etabin_number]+=thisweight*TrkSF*ele_trk_SF*ele_id_SF;
                myAccRecoIsoEffs.nels_reco_MC[ptbin_number][Etabin_number]++;

                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostElectronObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostElectronObj.reco_activity);

                myAccRecoIsoEffs.nels_reco_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight*TrkSF*ele_trk_SF_allreco*ele_id_SF_allreco;
                myAccRecoIsoEffs.nels_reco_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;

                (myBaseHistgram.h_id_genactivity_els)->Fill(myLostElectronObj.gen_activity,thisweight*TrkSF*ele_trk_SF*ele_id_SF);
                (myBaseHistgram.h_id_recoactivity_els)->Fill(myLostElectronObj.reco_activity,thisweight*TrkSF*ele_trk_SF_allreco*ele_id_SF_allreco);
              }

              if( myLostElectronObj.passIso )
              {
                int ptbin_number = Set_ptbin_number(myLostElectronObj.gen_pt);
                int acbin_number = Set_acbin_number(myLostElectronObj.gen_activity);

                myAccRecoIsoEffs.nels_iso[ptbin_number][acbin_number]+=thisweight*TrkSF*ele_trk_SF*ele_id_SF*ele_iso_SF;
                myAccRecoIsoEffs.nels_iso_MC[ptbin_number][acbin_number]++;
              
                //call another process for iso eff calculation, reset pt bin number for iso efficiency, as reco_pt
                int ptbin_number_allreco = Set_ptbin_number(myLostElectronObj.reco_pt);
                int acbin_number_allreco = Set_acbin_number(myLostElectronObj.reco_activity);
                
                myAccRecoIsoEffs.nels_iso_allreco[ptbin_number_allreco][acbin_number_allreco]+=thisweight*TrkSF*ele_trk_SF_allreco*ele_id_SF_allreco*ele_iso_SF_allreco;
                myAccRecoIsoEffs.nels_iso_MC_allreco[ptbin_number_allreco][acbin_number_allreco]++;
              }
              //check warning function when we calculate the efficienies!
              if ( nElectrons == 0 && myLostElectronObj.passIso )
              { 
		std::cout << "Warning: el pass iso but nElectrons=0!" << std::endl;
		std::cout << "reco_els_pt = " << myLostElectronObj.reco_pt << std::endl;
		std::cout << "reco_els_eta = " << myLostElectronObj.reco_eta << std::endl;

		std::cout << " isMu = " << myLostElectronObj.isMu << " isEl = " << myLostElectronObj.isEl << std::endl;
		std::cout << " gen_eta = " << myLostElectronObj.gen_eta << " gen_phi = " << myLostElectronObj.gen_phi << " gen_pt = " << myLostElectronObj.gen_pt << " gen_activity = " << myLostElectronObj.gen_activity << std::endl;
		std::cout << " reco_eta = " << myLostElectronObj.reco_eta << " reco_phi = " << myLostElectronObj.reco_phi << " reco_pt = " << myLostElectronObj.reco_pt << "reco_activity = " << myLostElectronObj.reco_activity << std::endl;
		std::cout << " passAcc = " << myLostElectronObj.passAcc << " passId = " << myLostElectronObj.passId << " passIso = " << myLostElectronObj.passIso << std::endl;
		std::cout << " doneAcc = " << myLostElectronObj.doneAcc << " doneId = " << myLostElectronObj.doneId << " doneIso = " << myLostElectronObj.doneIso << std::endl;
		std::cout << " reco_index = " << myLostElectronObj.reco_index << std::endl;

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

          //muon CS statistics
          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

          if( searchbin_id >= 0 )
          {
            myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbin_id]++;
            myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbin_id]+=thisweight;
          }
	  }
        }

	if (nElectrons == 0 && nMuons == 0 && (ngenmu==1 || ngenmu==2 || ngenel==1 || ngenel==2))
	{
          //const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          const int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);
          if( searchbin_id >= 0 )
          {
	    neventsSB[searchbin_id]+=thisweight;
	    (myBaseHistgram.h_iso_track_veto_nb_all)->Fill(nbottomjets,thisweight);
	    (myBaseHistgram.h_iso_track_veto_nt_all)->Fill(ntopjets,thisweight);
	    (myBaseHistgram.h_iso_track_veto_mt2_all)->Fill(MT2,thisweight);
	    (myBaseHistgram.h_iso_track_veto_met_all)->Fill(met,thisweight);
	    (myBaseHistgram.h_iso_track_veto_ht_all)->Fill(ht,thisweight);

	    if (nIsoTrks==0)
		{neventsSB_afterITV[searchbin_id]+=thisweight;
		 (myBaseHistgram.h_iso_track_veto_nb)->Fill(nbottomjets,thisweight);
                 (myBaseHistgram.h_iso_track_veto_nt)->Fill(ntopjets,thisweight);
                 (myBaseHistgram.h_iso_track_veto_mt2)->Fill(MT2,thisweight);
                 (myBaseHistgram.h_iso_track_veto_met)->Fill(met,thisweight);
                 (myBaseHistgram.h_iso_track_veto_ht)->Fill(ht,thisweight);
		}

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
      std::cout << "isotrack bin " << searchbinc << " pass " << neventsSB_afterITV[searchbinc] << " all " << neventsSB[searchbinc] << " eff " << isotrkeff[searchbinc] << " raw pass " << neventsSB_afterITV_MC[searchbinc] << " raw all " << neventsSB_MC[searchbinc] << std::endl;
      //std::cout << "isotrackeff_SB[" << searchbinc << "] = " << isotrkeff[searchbinc] << std::endl;
    }
  }

/*   std::cout << "13 () = " <<(neventsSB_afterITV[13]+neventsSB_afterITV[14])/(neventsSB[13]+neventsSB[14])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[13] + neventsSB_afterITV_MC[14] ,neventsSB_MC[13] + neventsSB_MC[14]) << std::endl;
   std::cout << "18 () = " <<(neventsSB_afterITV[18]+neventsSB_afterITV[19]+neventsSB_afterITV[20])/(neventsSB[18]+neventsSB[19]+neventsSB[20])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[18] + neventsSB_afterITV_MC[19] + neventsSB_afterITV_MC[20] ,neventsSB_MC[18] + neventsSB_MC[19] + neventsSB_MC[20]) << std::endl;
   std::cout << "19 () = " <<(neventsSB_afterITV[18]+neventsSB_afterITV[19]+neventsSB_afterITV[20])/(neventsSB[18]+neventsSB[19]+neventsSB[20])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[18] + neventsSB_afterITV_MC[19] + neventsSB_afterITV_MC[20] ,neventsSB_MC[18] + neventsSB_MC[19] + neventsSB_MC[20]) << std::endl;
   std::cout << "20 () = " <<(neventsSB_afterITV[18]+neventsSB_afterITV[19]+neventsSB_afterITV[20])/(neventsSB[18]+neventsSB[19]+neventsSB[20])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[18] + neventsSB_afterITV_MC[19] + neventsSB_afterITV_MC[20] ,neventsSB_MC[18] + neventsSB_MC[19] + neventsSB_MC[20]) << std::endl;
   std::cout << "29 () = " <<(neventsSB_afterITV[29]+neventsSB_afterITV[28])/(neventsSB[28]+neventsSB[29])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[29] + neventsSB_afterITV_MC[28] ,neventsSB_MC[29] + neventsSB_MC[28]) << std::endl;
   std::cout << "34 () = " <<(neventsSB_afterITV[34]+neventsSB_afterITV[35]+neventsSB_afterITV[36])/(neventsSB[34]+neventsSB[35]+neventsSB[36])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[34] + neventsSB_afterITV_MC[35] + neventsSB_afterITV_MC[36] ,neventsSB_MC[34] + neventsSB_MC[35] + neventsSB_MC[36]) << std::endl;
   std::cout << "35 () = " <<(neventsSB_afterITV[34]+neventsSB_afterITV[35]+neventsSB_afterITV[36])/(neventsSB[34]+neventsSB[35]+neventsSB[36])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[34] + neventsSB_afterITV_MC[35] + neventsSB_afterITV_MC[36] ,neventsSB_MC[34] + neventsSB_MC[35] + neventsSB_MC[36]) << std::endl;
   std::cout << "36 () = " <<(neventsSB_afterITV[34]+neventsSB_afterITV[35]+neventsSB_afterITV[36])/(neventsSB[34]+neventsSB[35]+neventsSB[36])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[34] + neventsSB_afterITV_MC[35] + neventsSB_afterITV_MC[36] ,neventsSB_MC[34] + neventsSB_MC[35] + neventsSB_MC[36]) << std::endl;
   std::cout << "40 () = " <<(neventsSB_afterITV[39]+neventsSB_afterITV[40])/(neventsSB[39]+neventsSB[40])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[39] + neventsSB_afterITV_MC[40] ,neventsSB_MC[39] + neventsSB_MC[40]) << std::endl;
   std::cout << "51 () = " <<(neventsSB_afterITV[51]+neventsSB_afterITV[50])/(neventsSB[51]+neventsSB[50])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[51] + neventsSB_afterITV_MC[50] ,neventsSB_MC[51] + neventsSB_MC[50]) << std::endl;
   std::cout << "57 () = " <<(neventsSB_afterITV[57]+neventsSB_afterITV[56])/(neventsSB[57]+neventsSB[56])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[56] + neventsSB_afterITV_MC[57] ,neventsSB_MC[56] + neventsSB_MC[57]) << std::endl;
   std::cout << "61 () = " <<(neventsSB_afterITV[61]+neventsSB_afterITV[60])/(neventsSB[61]+neventsSB[60])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[60] + neventsSB_afterITV_MC[61] ,neventsSB_MC[60] + neventsSB_MC[61]) << std::endl;
   std::cout << "67 () = " <<(neventsSB_afterITV[67]+neventsSB_afterITV[69])/(neventsSB[67]+neventsSB[69])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[67] + neventsSB_afterITV_MC[69] ,neventsSB_MC[67] + neventsSB_MC[69]) << std::endl;
   std::cout << "77 () = " <<(neventsSB_afterITV[77]+neventsSB_afterITV[76])/(neventsSB[77]+neventsSB[76])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[76] + neventsSB_afterITV_MC[77] ,neventsSB_MC[76] + neventsSB_MC[77]) << std::endl;
   std::cout << "79 () = " <<(neventsSB_afterITV[79]+neventsSB_afterITV[78])/(neventsSB[79]+neventsSB[78])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[79] + neventsSB_afterITV_MC[78] ,neventsSB_MC[79] + neventsSB_MC[78]) << std::endl;
   std::cout << "81 () = " <<(neventsSB_afterITV[81]+neventsSB_afterITV[80])/(neventsSB[81]+neventsSB[80])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[81] + neventsSB_afterITV_MC[80] ,neventsSB_MC[81] + neventsSB_MC[80]) << std::endl;
   std::cout << "83 () = " <<(neventsSB_afterITV[83]+neventsSB_afterITV[82])/(neventsSB[83]+neventsSB[82])<< " err = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[83] + neventsSB_afterITV_MC[82] ,neventsSB_MC[83] + neventsSB_MC[82]) << std::endl;
*/
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    std::cout << "isoTrackErr[" << searchbinc << "] = " << myAccRecoIsoEffs.get_stat_Error(neventsSB_afterITV_MC[searchbinc],neventsSB_MC[searchbinc]) << ";" << std::endl;
  }

  std::cout << "Muon CS # in Search Bins: " << std::endl;  
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
  std::cout << "bin " << searchbinc << " raw events " << myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbinc] << " weighted events " << myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbinc] << std::endl;
  } 

// debug
  /*for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    std::cout << "nmus_acc_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_acc_sb[searchbinc] << std::endl;
    std::cout << "nmus_acc_MC_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_acc_MC_sb[searchbinc] << std::endl;
    std::cout << "nmus_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_sb[searchbinc] << std::endl;
    std::cout << "nmus_MC_sb[" << searchbinc << "] = " << myAccRecoIsoEffs.nmus_MC_sb[searchbinc] << std::endl;
  }
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    std::cout << "nels_acc[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_acc[searchbinc] << std::endl;
    std::cout << "nels_acc_MC[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_acc_MC[searchbinc] << std::endl;
    std::cout << "nels[" << searchbinc << "] = " << myAccRecoIsoEffs.nels[searchbinc] << std::endl;
    std::cout << "nels_MC[" << searchbinc << "] = " << myAccRecoIsoEffs.nels_MC[searchbinc] << std::endl;
  }*/

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

  NTupleReader *tr =0;

  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(*tr, spec);

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  SearchBins theSearchBins("SB_v1_2017");

  std::cout << "Expectation: " << std::endl;
  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  {  
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader tr((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    //AnaFunctions::prepareTopTagger();
    //The passBaseline is registered here
    tr.registerFunction(&mypassBaselineFunc);    
 
    double thisweight = (*iter_TTJetsSampleInfos).weight;
    std::cout << (*iter_TTJetsSampleInfos).TTJetsTag << " Weight: " << thisweight << std::endl;
    int neventc=0;

    //while(tr.getNextEvent() && neventc<1000)
    while(tr.getNextEvent())
    {
      ++neventc;
      if(tr.getEvtNum()%20000 == 0) std::cout << tr.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

      const double genHT = tr.hasVar("genHT") ? tr.getVar<double>("genHT") : -999;
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromT_" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromTbar" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_DiLept") && genHT >= 600) continue; 
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-600to800" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-800to1200" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-1200to2500" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-2500toInf") && genHT < 600) continue; 
      double HT_sample_weight = thisweight;
      //0.43930872 is single lept BR from ttbar and 0.10614564 is dilepton BR from ttbar
      //if ((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-600to800" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-800to1200" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-1200to2500" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-2500toInf") HT_sample_weight = thisweight * (0.43930872+0.10614564);
      (myClosureHistgram.h_genHT)->Fill(genHT,HT_sample_weight);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
	  int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met,ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

          //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
          //59 bins
          int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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
  //std::cout << "Muon CS # in Search Bins: " << std::endl;
  //for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  //{
//    //myAccRecoIsoEffs.nevents_lept_exp_SB_MC[i_cal] = nevents_mus_exp_SB_MC[i_cal] + nevents_els_exp_SB_MC[i_cal];
//    //nevents_lept_pred_SB_MC[i_cal] = nevents_mus_pred_SB_MC[i_cal] + nevents_els_pred_SB_MC[i_cal];
//    //myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[i_cal] = nevents_mus_exp_SB_Normalized[i_cal] + nevents_els_exp_SB_Normalized[i_cal];
//    //myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] = nevents_mus_pred_SB_Normalized[i_cal] + nevents_els_pred_SB_Normalized[i_cal];
//    //nevents_lept_pred_SB_Normalized[i_cal] = nevents_lept_pred_SB_MC[i_cal]*scale;
//
  //std::cout << "idx: " << i_cal << "; MC Numbers: " << myAccRecoIsoEffs.nevents_mus_CS_SB_MC[i_cal] << "; Normalized Numbers: " << myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[i_cal] << std::endl;
  //std::cout << "Mus Exp MC Numbers: " << myAccRecoIsoEffs.nevents_mus_exp_SB_MC[i_cal] << "; Mus Exp Normalized Numbers: " << myAccRecoIsoEffs.nevents_mus_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Mus Pred Normalized Numbers: " << nevents_mus_pred_SB_Normalized[i_cal] << std::endl;
  //std::cout << "Els Exp MC Numbers: " << myAccRecoIsoEffs.nevents_els_exp_SB_MC[i_cal] << "; Els Exp Normalized Numbers: " << myAccRecoIsoEffs.nevents_els_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Els Pred Normalized Numbers: " << nevents_els_pred_SB_Normalized[i_cal] << std::endl;
//    std::cout << "Lept Exp MC Numbers: " << myAccRecoIsoEffs.nevents_lept_exp_SB_MC[i_cal] << "; Lept Exp Normalized Numbers: " << myAccRecoIsoEffs.nevents_lept_exp_SB_Normalized[i_cal] << std::endl;
//    //std::cout << "Lept Pred Normalized Numbers: " << nevents_lept_pred_SB_Normalized[i_cal] << std::endl;
  //}
  
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

void LoopLLPred( AccRecoIsoEffs& myAccRecoIsoEffs, TTJetsSampleWeight& myTTJetsSampleWeight, double resultspred[NSEARCH_BINS] )
{
  const bool storePlots = false;  //set to false when running LLsyst
  ClosureHistgram myClosureHistgram;
  if (storePlots) myClosureHistgram.BookHistgram("PredLL_mu_cs.root");
  //if (storePlots) myClosureHistgram.BookHistgram("v2_PredLL_data.root");

  NTupleReader *tr =0;

  //use class BaselineVessel in the SusyAnaTools/Tools/baselineDef.h file
  std::string spec = "lostlept";
  myBaselineVessel = new BaselineVessel(*tr, spec);
  //myBaselineVessel = new BaselineVessel(spec, "fastsim");

  size_t t0 = clock();
  std::vector<TTJetsSampleInfo>::iterator iter_TTJetsSampleInfos;

  SearchBins theSearchBins("SB_v1_2017");
   
  std::cout << "Prediction: " << std::endl;
  //int sidebandLowMT2=0;
  //int sidebandHighMT2=0;
  //std::vector< std::pair<int, int> > viMP;
  //std::map<std::pair<int, int>, std::pair<std::array<double,NSEARCH_BINS>, std::array<double,NSEARCH_BINS> > > mapSignalpred;

  for(iter_TTJetsSampleInfos = myTTJetsSampleWeight.TTJetsSampleInfos.begin(); iter_TTJetsSampleInfos != myTTJetsSampleWeight.TTJetsSampleInfos.end(); iter_TTJetsSampleInfos++)
  { 
    //use class NTupleReader in the SusyAnaTools/Tools/NTupleReader.h file
    NTupleReader trCS((*iter_TTJetsSampleInfos).chain);
    //initialize the type3Ptr defined in the customize.h
    //AnaFunctions::prepareTopTagger();
    //The passBaseline is registered here
    trCS.registerFunction(&mypassBaselineFunc);

    double thisweight = (*iter_TTJetsSampleInfos).weight;

    std::cout << (*iter_TTJetsSampleInfos).TTJetsTag << " Weight: " << thisweight << std::endl;

    int neventc=0;

    //while(trCS.getNextEvent() && neventc<1000)
    while(trCS.getNextEvent())
    {
      ++neventc;
      if(trCS.getEvtNum()%20000 == 0) std::cout << trCS.getEvtNum() << "\t" << ((clock() - t0)/1000000.0) << std::endl;

/*      //for MC only!!    start
      const double genHT = trCS.hasVar("genHT") ? trCS.getVar<double>("genHT") : -999;
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromT_" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_SingleLeptFromTbar" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_DiLept") && genHT >= 600) continue; 
      if (((*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-600to800" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-800to1200" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-1200to2500" || (*iter_TTJetsSampleInfos).TTJetsTag == "TTJets_HT-2500toInf") && genHT < 600) continue;
      //for MC only!!    end
*/
      bool passBaselinelostlept = trCS.getVar<bool>("passBaseline"+spec);
 
      //const double& SusyMotherMass  = trCS.getVar<double>("SusyMotherMass");
      //const double& SusyLSPMass     = trCS.getVar<double>("SusyLSPMass");

      //for data only!!  start
      std::vector<std::string> TriggerNames = trCS.getVec<std::string>("TriggerNames");
      std::vector<int> PassTrigger = trCS.getVec<int>("PassTrigger");
      bool foundTrigger = false;

      for(unsigned it=0; it<TriggerNames.size(); it++)
      {
	using namespace std;

	//if( TriggerNames[it].find("HLT_PFHT350_PFMET100_JetIdCleaned_v") || TriggerNames[it].find("HLT_PFHT350_PFMET100_NoiseCleaned_v") || TriggerNames[it].find("HLT_PFHT350_PFMET100_v*") )
	if( TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v")!= string::npos || TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET170_HBHECleaned_v") != string::npos || TriggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v") != string::npos || TriggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v")!= string::npos || TriggerNames[it].find("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v") != string::npos || TriggerNames[it].find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") != string::npos || TriggerNames[it].find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") != string::npos)

	{
	  if( PassTrigger[it] ) foundTrigger = true;
	}
      }
      if( !foundTrigger ) continue;
      //for data only!!  end

      bool passLeptVeto = trCS.getVar<bool>("passLeptVeto"+spec);
      bool passnJets = trCS.getVar<bool>("passnJets"+spec);
      bool passMET = trCS.getVar<bool>("passMET"+spec);
      bool passHT = trCS.getVar<bool>("passHT"+spec);
      bool passMT2 = trCS.getVar<bool>("passMT2"+spec);
      bool passTagger = trCS.getVar<bool>("passTagger"+spec);
      bool passBJets = trCS.getVar<bool>("passBJets"+spec);
      bool passNoiseEventFilter = trCS.getVar<bool>("passNoiseEventFilter"+spec);
      bool passdPhis = trCS.getVar<bool>("passdPhis"+spec);

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

	const std::vector<double> &metMagUp = trCS.getVec<double>("metMagUp");
	const std::vector<double> &metMagDown = trCS.getVec<double>("metMagDown");

	double metjecUp = metMagUp[1];//met Uncertainty
	double metjecLow = metMagDown[1];//met Uncertainty
	double metjerUp = metMagUp[0];//met Uncertainty
	double metjerLow = metMagDown[0];//met Uncertainty

	//std::cout << "met = " << met << " , metjecUp = " << metjecUp << " , metjecLow = " << metjecLow << " , metjerUp = " << metjerUp << " , metjerLow = " << metjerLow << std::endl;

	if (met<200) std::cout << "met<200" << std::endl;

        //muon CS
        
	if (use_muon_control_sample)
	{
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
	  //if (true)	//mtw test
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
            const int Etabin_number = Set_Etabin_number(fabs(reco_mus_eta));
            //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
            //59 bins
            int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

	    myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbin_id]++;
	    myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbin_id]+=thisweight;

	    //std::cout << "ttbar_mtwcorrfactor[0] = " << ttbar_mtwcorrfactor[0] << std::endl;
	    //mtwcorrfactor
	    EventWeight_mus = EventWeight_mus * ttbar_mtwcorrfactor[ptbin_number];	//mtw test
	    //dimuon correction factor
            EventWeight_mus = EventWeight_mus * ttbar_corrfactor_di_mus;

	    //if (applyisotrkveto) EventWeight_mus *= isotrackvetoeff;

            //muon prediction from muon CS
	    if (storePlots)
	    {
	    //Fill muon iso closure plots
	    (myClosureHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    //Fill muon id closure plots
	    (myClosureHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    //Fill muon acc closure plots
	    (myClosureHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    }
	    //Fill all muon closure plots
	    double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id];

	    if (storePlots)
	    {
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
	    }

            //total events flow for muons, prediction
            myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

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

              myAccRecoIsoEffs.nevents_mus_pred_iso_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_reco_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_acc_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

              //myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

	    //if (met>175.0 && met<=200.0)
	    //{
	      //std::cout << "met = " << met << std::endl;
	      //std::cout << "MT2 = " << MT2 << std::endl;
	      //if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_mus*EventWeight_mus;
	      //if (MT2>300.0) sidebandHighMT2 += EventWeight_all_mus*EventWeight_mus;
	    //}

            //begin to predict lost electrons from muon CS
            double EventWeight_els = thisweight;
            //mtwcorrfactor
            EventWeight_els = EventWeight_els * ttbar_mtwcorrfactor[ptbin_number];	//mtw test
            //dielectron correction factor
            EventWeight_els = EventWeight_els * ttbar_corrfactor_di_els;
  
	    //if (applyisotrkveto) EventWeight_els *= isotrackvetoeff;

            //electron prediction from muon CS
	    if (storePlots)
	    {
            //Fill electron iso closure plots
            (myClosureHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            //Fill electron id closure plots
            (myClosureHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            //Fill electron acc closure plots
            (myClosureHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
	    }
            //Fill all electron closure plots
            double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id];

	    if (storePlots)
	    {
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
	    }
            //total events flow for electrons, prediction
            myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;

	    //if (met>175.0 && met<=200.0)
	    //{
	    //  if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_els*EventWeight_els;
	    //  if (MT2>300.0) sidebandHighMT2 += EventWeight_all_els*EventWeight_els;
	    //}


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
            myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
          }//mtW5_els<125 (muon CS)
        }//nElectrons == 0 && nMuons == 1 (muon CS)
	}//end of bool: use_muon_control_sample
        
	//electron CS
        
	if (use_electron_control_sample)
	{
        if (nElectrons == 1 && nMuons == 0)
        {
          //counting the events for muon control sample
          myAccRecoIsoEffs.nevents_cs_mus++;
          //get reco muon variables
  	  std::vector<TLorentzVector> elesLVec = trCS.getVec<TLorentzVector>("elesLVec");
          std::vector<double> elespfActivity = trCS.getVec<double>("elespfActivity");
          std::vector<int> elesFlagVeto = trCS.getVec<int>("elesFlagVeto");
          std::vector<double> elesMiniIso = trCS.getVec<double>("elesMiniIso");

	//for electron purity
	//std::vector<int> PdgId = trCS.getVec<int>("genDecayPdgIdVec");
	//std::vector<TLorentzVector> gen_LVec = trCS.getVec<TLorentzVector>("genDecayLVec");

          double reco_eles_pt = -1, reco_eles_eta = 0, reco_eles_phi = 0;
          double activity = -1;
	  int nisoeles=0;
          for(unsigned int im = 0 ; im < elesLVec.size() ; im++)
          {
            if(elesFlagVeto[im] && elesLVec[im].Pt()>(AnaConsts::elesMiniIsoArr).minPt && fabs(elesLVec[im].Eta()) < (AnaConsts::elesMiniIsoArr).maxAbsEta && elesMiniIso[im] < (AnaConsts::elesMiniIsoArr).maxIsoEB )
	    {
              reco_eles_pt  = ( elesLVec.at(im) ).Pt();
              reco_eles_eta = ( elesLVec.at(im) ).Eta();
              reco_eles_phi = ( elesLVec.at(im) ).Phi();
              activity = elespfActivity.at(im);
	      ++nisoeles;
	    }
	  }
	  if (nisoeles!=1) std::cout << "Error: nisoeles!=1: nisoeles = " << nisoeles << std::endl;
          //if ( reco_mus_pt < 0 ) continue;

          double deltaphi_eles = DeltaPhi( reco_eles_phi , metphi );
          double mtW_eles = std::sqrt( 2.0 * reco_eles_pt * met * ( 1.0 - cos(deltaphi_eles) ) );
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

          if ( mtW_eles < 100.0 )
          {
	    //////////////////////////
	    // prediction computation
	    //////////////////////////
	    double EventWeight_mus = thisweight;
            const int njetsbin_number = Set_njetsbin_number(njets30);
            const int ptbin_number = Set_ptbin_number(reco_eles_pt);
            const int acbin_number = Set_acbin_number(activity);
	    const int htbin_number = Set_HTbin_number(ht);
	    const int Etabin_number = Set_Etabin_number(fabs(reco_eles_eta));
            //int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met );
            //59 bins
            int searchbin_id = theSearchBins.find_Binning_Index( nbottomjets , ntopjets , MT2, met, ht);

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

	    myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbin_id]++;
 	    myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbin_id]+=thisweight;

	    //std::cout << "ttbar_mtwcorrfactor[0] = " << ttbar_mtwcorrfactor[0] << std::endl;
	    //mtwcorrfactor
	    EventWeight_mus = EventWeight_mus * ttbar_mtwcorrfactor[ptbin_number];
	    //dimuon correction factor
            EventWeight_mus = EventWeight_mus * ttbar_corrfactor_di_mus;

	    //if (applyisotrkveto) EventWeight_mus *= isotrackvetoeff;
	//for electron purity	
	/*int temp=0;
        int n_gen_electron=0;
        for(int i=0; i < PdgId.size(); i++)
        {
        if(PdgId.at(i) == 11 || PdgId.at(i) == -11)
        {temp = i;
        n_gen_electron++;}
	}
        if (n_gen_electron=1 && elesLVec.at(0).DeltaR(gen_LVec.at(temp))< 0.4)
        {*/ 

        /*int n_gen_electron=0;
        for(int i=0; i < PdgId.size(); i++)
        {
        if(PdgId.at(i) == 11 || PdgId.at(i) == -11 && elesLVec.at(0).DeltaR(gen_LVec.at(i))< 0.2)
        {
        n_gen_electron++;}
        }
        if (n_gen_electron > 0)
	{*/

            //muon prediction from muon CS
	    if (storePlots)
	    {
	    //Fill muon iso closure plots
	    (myClosureHistgram.h_pred_mu_iso_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    //Fill muon id closure plots
	    (myClosureHistgram.h_pred_mu_id_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_id_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    //Fill muon acc closure plots
	    (myClosureHistgram.h_pred_mu_acc_met)->Fill(met, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_njets)->Fill(njets30, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ht)->Fill(ht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_mht)->Fill(mht, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus);
	    }
	    //Fill all muon closure plots
	    double EventWeight_all_mus = myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id];

	    if (storePlots)
	    {
		//EventWeight_mus=1;  //for purity study only!!!
		//EventWeight_all_mus=1;  //for purity study only!!!
	    (myClosureHistgram.h_pred_mu_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_njets)->Fill(njets30, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_mu_all_mt2)->Fill(MT2, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ht)->Fill(ht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_mht)->Fill(mht, EventWeight_all_mus*EventWeight_mus);
            (myClosureHistgram.h_pred_mu_all_ntopjets)->Fill(ntopjets, EventWeight_all_mus*EventWeight_mus);

	    (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_mus*EventWeight_mus);
	    (myClosureHistgram.h_pred_lept_all_2d_met_mupt)->Fill(deltaphi_eles, reco_eles_pt, EventWeight_all_mus*EventWeight_mus);
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
	    }

            //total events flow for muons, prediction
            myAccRecoIsoEffs.nevents_pred_all_mus += EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

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

              myAccRecoIsoEffs.nevents_mus_pred_iso_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_reco_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
              myAccRecoIsoEffs.nevents_mus_pred_acc_SB_Normalized[searchbin_id] += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

              //myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus/thisweight/thisweight;
              myAccRecoIsoEffs.nevents_mus_pred_SB_MC[searchbin_id] += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            }
            //here the error calculation is wrong...
            myAccRecoIsoEffs.nevents_pred_all_mus_err += EventWeight_all_mus*EventWeight_mus*EventWeight_all_mus*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_acc_mus_err += myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_id_mus_err += myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;
            myAccRecoIsoEffs.nevents_pred_iso_mus_err += myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus*myAccRecoIsoEffs.mus_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_mus;

	    //if (met>175.0 && met<=200.0)
	    //{
	      //std::cout << "met = " << met << std::endl;
	      //std::cout << "MT2 = " << MT2 << std::endl;
	      //if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_mus*EventWeight_mus;
	      //if (MT2>300.0) sidebandHighMT2 += EventWeight_all_mus*EventWeight_mus;
	    //}

            //begin to predict lost electrons from muon CS
            double EventWeight_els = thisweight;
            //mtwcorrfactor
            EventWeight_els = EventWeight_els * ttbar_mtwcorrfactor[ptbin_number];
            //dielectron correction factor
            EventWeight_els = EventWeight_els * ttbar_corrfactor_di_els;
  
	    //if (applyisotrkveto) EventWeight_els *= isotrackvetoeff;

            //electron prediction from muon CS
	    if (storePlots)
	    {
            //Fill electron iso closure plots
            (myClosureHistgram.h_pred_el_iso_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_iso_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            //Fill electron id closure plots
            (myClosureHistgram.h_pred_el_id_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_id_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            //Fill electron acc closure plots
            (myClosureHistgram.h_pred_el_acc_met)->Fill(met, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_njets)->Fill(njets30, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mt2)->Fill(MT2, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ht)->Fill(ht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_mht)->Fill(mht, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
            (myClosureHistgram.h_pred_el_acc_ntopjets)->Fill(ntopjets, myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els);
	    }
            //Fill all electron closure plots
            double EventWeight_all_els = myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id] + myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id];

	    if (storePlots)
	    {
		//EventWeight_els=1;  //for purity study only!!!
		//EventWeight_all_els=1;  //for purity study only!!!
            (myClosureHistgram.h_pred_el_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_njets)->Fill(njets30, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mt2)->Fill(MT2, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ht)->Fill(ht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_mht)->Fill(mht, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_el_all_ntopjets)->Fill(ntopjets, EventWeight_all_els*EventWeight_els);

            (myClosureHistgram.h_pred_lept_all_met)->Fill(met, EventWeight_all_els*EventWeight_els);
            (myClosureHistgram.h_pred_lept_all_2d_met_mupt)->Fill(deltaphi_eles, reco_eles_pt, EventWeight_all_els*EventWeight_els);
	    (myClosureHistgram.h_pred_el_all_2d_met_njets)->Fill(met, njets30);
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
	    }
            //total events flow for electrons, prediction
            myAccRecoIsoEffs.nevents_pred_all_els += EventWeight_all_els*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_acc_els += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;

	    //if (met>175.0 && met<=200.0)
	    //{
	    //  if (MT2>200.0 && MT2<300.0) sidebandLowMT2 += EventWeight_all_els*EventWeight_els;
	    //  if (MT2>300.0) sidebandHighMT2 += EventWeight_all_els*EventWeight_els;
	    //}


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
            myAccRecoIsoEffs.nevents_pred_acc_els_err += myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_acc[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_id_els_err += myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_reco[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
            myAccRecoIsoEffs.nevents_pred_iso_els_err += myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els*myAccRecoIsoEffs.els_EventWeight_iso[njetsbin_number][ptbin_number][acbin_number][Etabin_number][searchbin_id]*EventWeight_els;
         //}//end of electron purity
	 }//mtW5_els<125 (muon CS)
        }//nElectrons == 0 && nMuons == 1 (muon CS)
	}//end of bool: use_electron_control_sample
      }//baseline_nolepveto
    }//TTJets samples class
  }//end of second loop

  //std::cout << "sidebandLowMT2 = " << sidebandLowMT2 << std::endl;
  //std::cout << "sidebandHighMT2 = " << sidebandHighMT2 << std::endl;

  // TF muon CS
std::array <double, 84> aveTFfromMC_mu_cs={0.65452, 0.53775, 0.58909, 0.41141, 0.30614, 0.62969, 0.65089, 0.65674, 0.50386, 0.45485, 0.61707, 0.51231, 0.48089, 0.28145, 0.33864, 0.48679, 0.55109, 0.54004, 0.90584, 0.54099, 0.52136, 0.59685, 0.56194, 0.52321, 0.67406, 0.53257, 0.60314, 0.85751, 0.57650, 0.59478, 0.60984, 0.48227, 0.32457, 0.48768, 0.69924, 0.61275, 0.51995, 0.64178, 0.63747, 0.68836, 0.43362, 0.54726, 0.54114, 0.56058, 0.64455, 0.47306, 0.53588, 0.55759, 0.32667, 0.29986, 0.32745, 0.38446, 0.64409, 0.51746, 0.45185, 0.49802, 0.99901, 0.35219, 0.31627, 0.31500, 0.29196, 0.35471, 0.53624, 0.59835, 0.51133, 0.44067, 0.47647, 0.56281, 0.38735, 0.77547, 0.41253, 0.38302, 0.70333, 0.62000, 0.32163, 0.37719, 0.31216, 0.41283, 0.36693, 0.36078, 0.20553, 0.34576, 0.29329, 0.22609};

  // TF ele CS
std::array <double, 84> aveTFfromMC_el_cs={0.76111, 0.64307, 0.81815, 0.49388, 0.40288, 0.71423, 0.71977, 0.86432, 0.52828, 0.54922, 0.75247, 0.67931, 0.60780, 0.33887, 0.42287, 0.52772, 0.60815, 0.55598, 0.42422, 0.68920, 0.51671, 0.70016, 0.63413, 0.55442, 0.97103, 0.62896, 0.67011, 1.01647, 0.61673, 0.52730, 0.59042, 0.56515, 0.33646, 0.43478, 0.78056, 0.82275, 0.47275, 0.81923, 0.68846, 0.72558, 0.78916, 0.58734, 0.61198, 0.68876, 0.76163, 0.57087, 0.66599, 0.67437, 0.38373, 0.46969, 0.32441, 0.60692, 0.83025, 0.72363, 0.51538, 0.62168, 0.92336, 0.67919, 0.35023, 0.44408, 0.33970, 0.45095, 0.63113, 0.64244, 0.79089, 0.59042, 0.67512, 0.60941, 0.39625, 0.84927, 0.46653, 0.58492, 0.72392, 0.64514, 0.51673, 0.37487, 0.35070, 0.41767, 0.48407, 0.38646, 0.32989, 0.35159, 0.63034, 0.31140};

std::array <double, 84> aveTFfromMC;

if (use_muon_control_sample) aveTFfromMC = aveTFfromMC_mu_cs;
if (use_electron_control_sample) aveTFfromMC = aveTFfromMC_el_cs;

// sys unc on mu cs
std::array <double, 84> psystup_mu_cs = {0.141396, 0.164841, 0.281423, 0.550958, 0.648367, 0.141798, 0.168322, 0.327499, 0.41299, 0.234346, 0.206379, 0.244513, 0.343945, 1.03763, 0.427655, 0.315673, 0.294014, 0.670768, 0.987603, 0.8371, 0.54593, 0.135934, 0.192384, 0.311753, 0.404397, 0.871196, 0.277264, 0.295366, 0.492708, 0.642887, 0.545674, 0.299616, 0.353631, 0.729562, 0.927877, 1.04834, 8.01852, 0.186901, 0.229194, 0.421644, 0.750882, 0.223623, 0.31601, 0.526376, 1.2957, 0.336038, 0.577684, 0.513048, 0.176665, 0.346984, 1.01585, 0.5955, 0.187025, 0.30519, 0.606232, 0.528413, 0.420343, 1.14631, 0.183418, 0.40556, 0.552962, 0.744836, 0.201447, 0.355197, 0.499681, 0.301073, 0.657972, 1.54633, 0.836082, 0.684243, 0.297229, 0.677555, 1.10599, 0.382347, 0.444965, 0.405665, 0.822449, 1.02792, 0.734019, 0.656343, 0.594371, 0.651944, 1.13586, 1.30939};
std::array <double, 84> psystdown_el_cs = {0.149493, 0.172479, 0.286882, 0.552634, 0.65103, 0.15044, 0.175667, 0.332019, 0.415796, 0.240966, 0.212587, 0.25, 0.347237, 1.03983, 0.432305, 0.321864, 0.299818, 0.674733, 0.98935, 0.838822, 0.548159, 0.144552, 0.198913, 0.31716, 0.40666, 0.878702, 0.282733, 0.29907, 0.496128, 0.645953, 0.547966, 0.305743, 0.355815, 0.731546, 0.933677, 1.05362, 8.01878, 0.193064, 0.233899, 0.424834, 0.752536, 0.229386, 0.318895, 0.529877, 1.29668, 0.341896, 0.580697, 0.519686, 0.189061, 0.355912, 1.01742, 0.598775, 0.194091, 0.31108, 0.609808, 0.532074, 0.422819, 1.14844, 0.196324, 0.410935, 0.557477, 0.749531, 0.207757, 0.359112, 0.501791, 0.305005, 0.660483, 1.54707, 0.838382, 0.685479, 0.303382, 0.679131, 1.10882, 0.388115, 0.452696, 0.409691, 0.826608, 1.0327, 0.737954, 0.65902, 0.598382, 0.655214, 1.13834, 1.3112};

// sys unc on el cs
std::array <double, 84> psystup_el_cs = {0.154295, 0.191414, 0.362809, 0.540923, 0.522177, 0.150224, 0.18141, 0.373475, 0.44063, 0.281002, 0.215485, 0.249879, 0.327023, 0.993055, 0.473444, 0.319877, 0.302917, 0.680006, 1.04907, 0.718578, 0.556834, 0.153884, 0.211554, 0.325141, 0.47482, 1.06294, 0.250852, 0.394599, 0.471205, 0.636345, 0.578839, 0.339602, 0.383172, 0.725312, 0.882582, 1.21453, 8.02166, 0.201669, 0.267404, 0.477176, 1.12334, 0.246705, 0.353836, 0.511175, 1.04908, 0.334567, 0.614135, 0.516546, 0.19208, 0.521113, 1.02876, 0.581841, 0.235654, 0.413865, 0.58624, 0.483952, 0.415642, 1.02568, 0.192338, 0.467486, 0.460686, 0.68916, 0.224666, 0.38204, 0.553185, 0.368383, 0.580351, 1.50633, 0.719894, 0.929591, 0.302475, 0.823748, 1.11267, 0.410347, 0.707959, 0.416347, 0.791454, 1.0399, 0.746679, 0.624194, 0.595999, 0.629583, 1.00861, 1.4951};
std::array <double, 84> psystdown_el_cs = {0.161887, 0.19875, 0.366249, 0.54414, 0.525293, 0.158665, 0.188352, 0.377104, 0.44367, 0.286191, 0.221962, 0.255088, 0.332355, 0.99485, 0.476667, 0.324602, 0.307886, 0.682786, 1.0511, 0.721829, 0.558917, 0.161826, 0.218304, 0.329587, 0.476709, 1.06646, 0.256893, 0.397064, 0.474211, 0.638404, 0.580851, 0.345704, 0.385644, 0.727238, 0.886169, 1.21832, 8.0222, 0.207716, 0.271371, 0.480688, 1.12461, 0.25377, 0.358627, 0.513713, 1.0499, 0.339432, 0.616853, 0.519375, 0.202692, 0.525846, 1.03051, 0.585382, 0.240967, 0.41756, 0.590592, 0.487342, 0.418706, 1.02691, 0.203409, 0.4718, 0.466704, 0.693208, 0.231233, 0.387406, 0.555583, 0.371799, 0.582851, 1.50892, 0.722541, 0.930919, 0.307333, 0.825173, 1.11473, 0.415139, 0.711888, 0.420162, 0.794276, 1.04178, 0.74911, 0.628642, 0.598487, 0.63184, 1.01015, 1.4958};

std::array <double, 84> psystup;

if (use_muon_control_sample) psystup = psystup_mu_cs;
if (use_electron_control_sample) psystup = psystup_el_cs;

std::array <double, 84> psystdown;

if (use_muon_control_sample) psystdown = psystdown_mu_cs;
if (use_electron_control_sample) psystdown = psystdown_el_cs;

  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {

    if (use_electron_control_sample)
    {myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] = myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] * electron_purity;
    myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal] = myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal] * electron_purity;}

    myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal] = myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] + myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal];
//    h_cs_mus_sb->SetBinContent( i_cal+1 , nevents_mus_CS_SB_Normalized[i_cal] );
    if (storePlots)
    {	    
    myClosureHistgram.h_pred_mu_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_mus_pred_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_mu_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal]));
    myClosureHistgram.h_pred_el_sb->SetBinContent( i_cal+1 , myAccRecoIsoEffs.nevents_els_pred_SB_Normalized[i_cal] );
    myClosureHistgram.h_pred_el_sb->SetBinError( i_cal+1 , std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]));

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
    }
    //std::cout << "N events [" << i_cal << "] = " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal] << " , pred = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << std::endl;
 
    // data card
    //std::cout << "\n";
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal];
    //std::cout << " " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(60144642+59816364+30498962)/2262.946/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/4004.345/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/8000.0/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/7647.637518921/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/12877.0846508279992/831.76;
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(53057043+60494823+30682233)/36352.970569733/831.76;
    //std::cout << " 2: " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal];
    //std::cout << " 3: " << isoTrackEff_SB[i_cal];
    //std::cout << " 4: " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
    // 1.83333 c'est 1/(0.43930872+0.10614564)
    //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal];
    //std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal];
    //relative stat unc
    //if (myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]==0.0) std::cout << " 0";
    //else std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))/myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal];
    //if (myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]==0.0) std::cout << " 0";
    //else std::cout << " " << std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]);

    //std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal];
       //                //std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ",";

    // syst
    resultspred[i_cal]=myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal];
    //std::cout << "systStatus = " << systStatus << std::endl;
    //if (systStatus==1) std::cout << "Pred[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "Pred[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredEUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredSystUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;
    //std::cout << "PredESystUp[" << i_cal << "] = " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << ";" << std::endl;

    //std::cout << "uncUpPred[" << i_cal << "] = " << std::sqrt(((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])*((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])+1.8333*aveTFfromMC[i_cal]*1.8333*aveTFfromMC[i_cal]) << ";" << std::endl;
    //std::cout << "uncDownPred[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal] << ";" << std::endl;

    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(60144642+59816364+30498962)/2262.946/831.76 << ";" << std::endl;
    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.8*(59654914+51873969+30587326)/2262.946/831.76 << ";" << std::endl;
    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/8000.0/831.76 << ";" << std::endl;
    //std::cout << "aveTFErr[" << i_cal << "] = " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]*1.83333*(59654914+51873969+30587326)/7647.637518921/831.76 << ";" << std::endl;



    if (true)
    {
    if (i_cal == 0) std::cout  <<  "0 &               1 &               1 &         200-300 &         250-400 & ";
    if (i_cal == 1) std::cout  <<  "1 &               1 &               1 &         200-300 &         400-500 & ";
    if (i_cal == 2) std::cout  <<  "2 &               1 &               1 &         200-300 &         500-600 & ";
    if (i_cal == 3) std::cout  <<  "3 &               1 &               1 &         200-300 &         600-750 & ";
    if (i_cal == 4) std::cout  <<  "4 &               1 &               1 &         200-550 &            750+ & ";
    if (i_cal == 5) std::cout  <<  "5 &               1 &               1 &         300-400 &         250-400 & ";
    if (i_cal == 6) std::cout  <<  "6 &               1 &               1 &         300-400 &         400-500 & ";
    if (i_cal == 7) std::cout  <<  "7 &               1 &               1 &         300-400 &         500-600 & ";
    if (i_cal == 8) std::cout  <<  "8 &               1 &               1 &         300-400 &         600-750 & ";
    if (i_cal == 9) std::cout  <<  "9 &               1 &               1 &         400-550 &         250-400 & ";
    if (i_cal == 10) std::cout  <<  "10 &               1 &               1 &         400-550 &         400-500 & ";
    if (i_cal == 11) std::cout  <<  "11 &               1 &               1 &         400-550 &         500-600 & ";
    if (i_cal == 12) std::cout  <<  "12 &               1 &               1 &         400-550 &         600-750 & ";
    if (i_cal == 13) std::cout  <<  "13 &               1 &               1 &         550-750 &         250-400 & ";
    if (i_cal == 14) std::cout  <<  "14 &               1 &               1 &         550-750 &         400-500 & ";
    if (i_cal == 15) std::cout  <<  "15 &               1 &               1 &         550-750 &         500-600 & ";
    if (i_cal == 16) std::cout  <<  "16 &               1 &               1 &         550-750 &         600-750 & ";
    if (i_cal == 17) std::cout  <<  "17 &               1 &               1 &         550-750 &            750+ & ";
    if (i_cal == 18) std::cout  <<  "18 &               1 &               1 &            750+ &         250-600 & ";
    if (i_cal == 19) std::cout  <<  "19 &               1 &               1 &            750+ &         600-750 & ";
    if (i_cal == 20) std::cout  <<  "20 &               1 &               1 &            750+ &            750+ & ";
    if (i_cal == 21) std::cout  <<  "21 &               1 &               2 &         200-350 &         250-400 & ";
    if (i_cal == 22) std::cout  <<  "22 &               1 &               2 &         200-350 &         400-500 & ";
    if (i_cal == 23) std::cout  <<  "23 &               1 &               2 &         200-350 &         500-600 & ";
    if (i_cal == 24) std::cout  <<  "24 &               1 &               2 &         200-350 &         600-750 & ";
    if (i_cal == 25) std::cout  <<  "25 &               1 &               2 &         200-650 &            750+ & ";
    if (i_cal == 26) std::cout  <<  "26 &               1 &               2 &         350-450 &         250-400 & ";
    if (i_cal == 27) std::cout  <<  "27 &               1 &               2 &         350-450 &         400-500 & ";
    if (i_cal == 28) std::cout  <<  "28 &               1 &               2 &         350-450 &         500-600 & ";
    if (i_cal == 29) std::cout  <<  "29 &               1 &               2 &         350-450 &         600-750 & ";
    if (i_cal == 30) std::cout  <<  "30 &               1 &               2 &         450-650 &         250-400 & ";
    if (i_cal == 31) std::cout  <<  "31 &               1 &               2 &         450-650 &         400-500 & ";
    if (i_cal == 32) std::cout  <<  "32 &               1 &               2 &         450-650 &         500-600 & ";
    if (i_cal == 33) std::cout  <<  "33 &               1 &               2 &         450-650 &         600-750 & ";
    if (i_cal == 34) std::cout  <<  "34 &               1 &               2 &            650+ &         250-600 & ";
    if (i_cal == 35) std::cout  <<  "35 &               1 &               2 &            650+ &         600-750 & ";
    if (i_cal == 36) std::cout  <<  "36 &               1 &               2 &            650+ &            750+ & ";
    if (i_cal == 37) std::cout  <<  "37 &               1 &              3+ &        300-1000 &         250-350 & ";
    if (i_cal == 38) std::cout  <<  "38 &               1 &              3+ &        300-1000 &         350-450 & ";
    if (i_cal == 39) std::cout  <<  "39 &               1 &              3+ &        300-1000 &         450-550 & ";
    if (i_cal == 40) std::cout  <<  "40 &               1 &              3+ &        300-1000 &            550+ & ";
    if (i_cal == 41) std::cout  <<  "41 &               1 &              3+ &       1000-1500 &         250-350 & ";
    if (i_cal == 42) std::cout  <<  "42 &               1 &              3+ &       1000-1500 &         350-450 & ";
    if (i_cal == 43) std::cout  <<  "43 &               1 &              3+ &       1000-1500 &         450-550 & ";
    if (i_cal == 44) std::cout  <<  "44 &               1 &              3+ &       1000-1500 &            550+ & ";
    if (i_cal == 45) std::cout  <<  "45 &               1 &              3+ &           1500+ &         250-350 & ";
    if (i_cal == 46) std::cout  <<  "46 &               1 &              3+ &           1500+ &         350-550 & ";
    if (i_cal == 47) std::cout  <<  "47 &               1 &              3+ &           1500+ &            550+ & ";
    if (i_cal == 48) std::cout  <<  "48 &               2 &               1 &         200-300 &         250-350 & ";
    if (i_cal == 49) std::cout  <<  "49 &               2 &               1 &         200-300 &         350-450 & ";
    if (i_cal == 50) std::cout  <<  "50 &               2 &               1 &         200-300 &         450-600 & ";
    if (i_cal == 51) std::cout  <<  "51 &               2 &               1 &         200-450 &            600+ & ";
    if (i_cal == 52) std::cout  <<  "52 &               2 &               1 &         300-450 &         250-350 & ";
    if (i_cal == 53) std::cout  <<  "53 &               2 &               1 &         300-450 &         350-450 & ";
    if (i_cal == 54) std::cout  <<  "54 &               2 &               1 &         300-450 &         450-600 & ";
    if (i_cal == 55) std::cout  <<  "55 &               2 &               1 &            450+ &         250-450 & ";
    if (i_cal == 56) std::cout  <<  "56 &               2 &               1 &            450+ &         450-600 & ";
    if (i_cal == 57) std::cout  <<  "57 &               2 &               1 &            450+ &            600+ & ";
    if (i_cal == 58) std::cout  <<  "58 &               2 &               2 &         200-300 &         250-350 & ";
    if (i_cal == 59) std::cout  <<  "59 &               2 &               2 &         200-300 &         350-450 & ";
    if (i_cal == 60) std::cout  <<  "60 &               2 &               2 &         200-300 &         450-600 & ";
    if (i_cal == 61) std::cout  <<  "61 &               2 &               2 &         200-400 &            600+ & ";
    if (i_cal == 62) std::cout  <<  "62 &               2 &               2 &         300-400 &         250-350 & ";
    if (i_cal == 63) std::cout  <<  "63 &               2 &               2 &         300-400 &         350-450 & ";
    if (i_cal == 64) std::cout  <<  "64 &               2 &               2 &         300-400 &         450-600 & ";
    if (i_cal == 65) std::cout  <<  "65 &               2 &               2 &         400-500 &         250-450 & ";
    if (i_cal == 66) std::cout  <<  "66 &               2 &               2 &         400-500 &         450-600 & ";
    if (i_cal == 67) std::cout  <<  "67 &               2 &               2 &            400+ &            600+ & ";
    if (i_cal == 68) std::cout  <<  "68 &               2 &               2 &            500+ &         250-450 & ";
    if (i_cal == 69) std::cout  <<  "69 &               2 &               2 &            500+ &         450-600 & ";
    if (i_cal == 70) std::cout  <<  "70 &               2 &              3+ &         300-900 &         250-350 & ";
    if (i_cal == 71) std::cout  <<  "71 &               2 &              3+ &         300-900 &         350-500 & ";
    if (i_cal == 72) std::cout  <<  "72 &               2 &              3+ &        300-1300 &            500+ & ";
    if (i_cal == 73) std::cout  <<  "73 &               2 &              3+ &        900-1300 &         250-350 & ";
    if (i_cal == 74) std::cout  <<  "74 &               2 &              3+ &        900-1300 &         350-500 & ";
    if (i_cal == 75) std::cout  <<  "75 &               2 &              3+ &           1300+ &         250-350 & ";
    if (i_cal == 76) std::cout  <<  "76 &               2 &              3+ &           1300+ &         350-500 & ";
    if (i_cal == 77) std::cout  <<  "77 &               2 &              3+ &           1300+ &            500+ & ";
    if (i_cal == 78) std::cout  <<  "78 &              3+ &               1 &            200+ &         250-350 & ";
    if (i_cal == 79) std::cout  <<  "79 &              3+ &               1 &            200+ &            350+ & ";
    if (i_cal == 80) std::cout  <<  "80 &              3+ &               2 &            200+ &         250-400 & ";
    if (i_cal == 81) std::cout  <<  "81 &              3+ &               2 &            200+ &            400+ & ";
    if (i_cal == 82) std::cout  <<  "82 &              3+ &              3+ &            200+ &         250-350 & ";
    if (i_cal == 83) std::cout  <<  "83 &              3+ &              3+ &            200+ &            350+ & ";

    std::cout.precision(3);
    std::cout << std::fixed;

    std::cout << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal] << " $^{+" << std::sqrt(((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])*((std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal])+1.8333*aveTFfromMC[i_cal]*1.8333*aveTFfromMC[i_cal]) << "}_{-" << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))*isoTrackEff_SB[i_cal] << "}$ $^{+" << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]*psystup[i_cal] << "}_{-" << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]*psystdown[i_cal] << "}$" << " \\\\ " << 
std::endl;
    std::cout << "\\hline" << std::endl;
    }
  }

  std::cout  << std::endl;

  std::cout.precision(5);
  std::cout << std::fixed;

  std::cout << "\n data card # Predicted central numbers " << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ ) {std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal];}

  std::cout << "\n data card # Control sample events " << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ ) { std::cout << " " << myAccRecoIsoEffs.nevents_cs_mus_sb[i_cal]; }

  std::cout << "\n data card  # Average weight " << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ ) 
  {if (myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[i_cal] == 0) std::cout << " 0";
  else  std::cout << " " << myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]*isoTrackEff_SB[i_cal]/myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[i_cal]; }

  std::cout << "\n data card stat unc" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ ) 
  {if (myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal]==0.0) std::cout << " 0";
  else std::cout << " " << (std::sqrt(myAccRecoIsoEffs.nevents_mus_pred_SB_MC[i_cal])+std::sqrt(myAccRecoIsoEffs.nevents_els_pred_SB_MC[i_cal]))/myAccRecoIsoEffs.nevents_lept_pred_SB_Normalized[i_cal];}

  std::cout << "\n Muon CS # in Search Bins: " << std::endl;
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
  std::cout << "bin " << searchbinc << " raw events " << myAccRecoIsoEffs.nevents_mus_CS_SB_MC[searchbinc] << " weighted events " << myAccRecoIsoEffs.nevents_mus_CS_SB_Normalized[searchbinc] << std::endl;
  }

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
  if (storePlots)
  {
  (myClosureHistgram.oFile)->Write();
  std::cout << "write is done" << std::endl;
  (myClosureHistgram.oFile)->Close();
  std::cout << "close is done" << std::endl;
  }
  return ;
}

void LoopLLSyst( TTJetsSampleWeight& myTTJetsSampleWeight )
{
  //std::cout <<__LINE__<<std::endl;
  //https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#2016_analyses_80X_Data_and_MC_IC
  std::cout << "----Computing syst: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs;
  double resultspred[NSEARCH_BINS]={0};
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspred );
  myAccRecoIsoEffs.init();

  const bool domureco=false;
  const bool domuiso=false;
  const bool dodimu=false;
  const bool dodie=true;
  const bool doereco=false;
  const bool doeiso=false;
  const bool doacc=false;

  if (domureco)
  {
  std::cout << "Computing syst mu reco stat up: " << std::endl;
  double resultspredStatUp[NSEARCH_BINS]={0};
  double MuRecoEff[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_up[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_dn[PT_BINS][AC_BINS] = {{0}};
  TFile *fin = TFile::Open("v2_Effs2dPlots.root");
  TH2D * murecoeff;
  murecoeff = (TH2D*)fin->Get("mus_recoeffs")->Clone();

/*  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //std::cout << "murecoeff->GetBinContent(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinContent(i+1,j+2) << std::endl;
      //std::cout << "murecoeff->GetBinError(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinError(i+1,j+2) << std::endl;
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      MuRecoEff_Stat_Unc_up[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]>1.0) MuRecoEff_Stat_Unc_up[i][j]=1.0;
      MuRecoEff_Stat_Unc_dn[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]<0.0) MuRecoEff_Stat_Unc_up[i][j]=0.0;
      //std::cout << "MuRecoEff[" << i << "][" << j << "] = " << MuRecoEff[i][j] << std::endl;
      //std::cout << "ttbar_mus_recoeff = " << ttbar_mus_recoeff[i][j] << std::endl;
    }
  }
*/
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      ttbar_mus_recoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (ttbar_mus_recoeff[i][j] > 1.0) ttbar_mus_recoeff[i][j] = 1.0;
      //ttbar_mus_recoeff[i][j]=MuRecoEff_Stat_Unc_up[i][j];
    }
  }

  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu reco stat dn: " << std::endl;
  double resultspredStatDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_recoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (ttbar_mus_recoeff[i][j] < 0.0) {ttbar_mus_recoeff[i][j] = 0.0;
      std::cout << "mu reco eff [" << i << "][" << j << "] set to 0" << std::endl;}
      //ttbar_mus_recoeff[i][j]=MuRecoEff_Stat_Unc_dn[i][j];
    }
  }

  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatDn );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu reco syst up: " << std::endl;
  double resultspredSystUp[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_recoeff[i][j]=MuRecoEff[i][j]+0.03;
      if (ttbar_mus_recoeff[i][j]>1.0) ttbar_mus_recoeff[i][j]=1.0;
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystUp );
  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu reco syst dn: " << std::endl;
  double resultspredSystDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_recoeff[i][j]=MuRecoEff[i][j]-0.03;
      if (ttbar_mus_recoeff[i][j]<0.0) {ttbar_mus_recoeff[i][j]=0.0;
      std::cout << "mu reco eff [" << i << "][" << j << "] set to 0" << std::endl;}
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystDn );
  myAccRecoIsoEffs.init();

  std::cout << "mu reco syst up:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //std::cout << "resultspred[" << i_cal << "] = " << resultspred[i_cal] << std::endl;
    //std::cout << "resultspredStatUp[i_cal] = " << resultspredStatUp[i_cal] << std::endl;
    //std::cout << "resultspredStatDn[i_cal] = " << resultspredStatDn[i_cal] << std::endl;
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatUp[i_cal];

    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal];
    
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatUp[i_cal])*(resultspred[i_cal]-resultspredStatUp[i_cal]) + (resultspred[i_cal]-resultspredSystUp[i_cal])*(resultspred[i_cal]-resultspredSystUp[i_cal])) / resultspred[i_cal] << " ";

    else std::cout << "nan ";
  }
  std::cout << std::endl;

  std::cout << "mu reco syst down:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatDn[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatDn[i_cal])*(resultspred[i_cal]-resultspredStatDn[i_cal]) + (resultspred[i_cal]-resultspredSystDn[i_cal])*(resultspred[i_cal]-resultspredSystDn[i_cal])) / resultspred[i_cal] << " ";
    else std::cout << "nan ";
  }
  std::cout << std::endl;

  // restore default:
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_recoeff[i][j]=MuRecoEff[i][j];
    }
  }

  }

  if (domuiso)
  {
  std::cout << "Computing syst mu iso stat up: " << std::endl;
  double resultspredStatUp[NSEARCH_BINS]={0};
  double MuRecoEff[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_up[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_dn[PT_BINS][AC_BINS] = {{0}};
  TFile *fin = TFile::Open("v2_Effs2dPlots.root");
  TH2D * murecoeff;
  murecoeff = (TH2D*)fin->Get("mus_isoeffs")->Clone();

/*  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //std::cout << "murecoeff->GetBinContent(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinContent(i+1,j+2) << std::endl;
      //std::cout << "murecoeff->GetBinError(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinError(i+1,j+2) << std::endl;
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      MuRecoEff_Stat_Unc_up[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]>1.0) MuRecoEff_Stat_Unc_up[i][j]=1.0;
      MuRecoEff_Stat_Unc_dn[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]<0.0) MuRecoEff_Stat_Unc_up[i][j]=0.0;
      //std::cout << "MuRecoEff[" << i << "][" << j << "] = " << MuRecoEff[i][j] << std::endl;
      //std::cout << "ttbar_mus_isoeff = " << ttbar_mus_isoeff[i][j] << std::endl;
    }
  }
*/
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      ttbar_mus_isoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (ttbar_mus_isoeff[i][j] > 1.0) ttbar_mus_isoeff[i][j] = 1.0;
      //ttbar_mus_isoeff[i][j]=MuRecoEff_Stat_Unc_up[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu iso stat dn: " << std::endl;
  double resultspredStatDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_isoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (ttbar_mus_isoeff[i][j] < 0.0) {ttbar_mus_isoeff[i][j] = 0.0;
      std::cout << "mu iso eff [" << i << "][" << j << "] set to 0" << std::endl;}
      //ttbar_mus_isoeff[i][j]=MuRecoEff_Stat_Unc_dn[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatDn );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu iso syst up: " << std::endl;
  double resultspredSystUp[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_isoeff[i][j]=MuRecoEff[i][j]+0.01;
      if (ttbar_mus_isoeff[i][j]>1.0) ttbar_mus_isoeff[i][j]=1.0;
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst mu iso syst dn: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs5;
  double resultspredSystDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_isoeff[i][j]=MuRecoEff[i][j]-0.01;
      if (ttbar_mus_isoeff[i][j]<0.0) {ttbar_mus_isoeff[i][j]=0.0;
      std::cout << "mu iso eff [" << i << "][" << j << "] set to 0" << std::endl;}
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystDn );

  myAccRecoIsoEffs.init();

  std::cout << "mu iso syst up:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //std::cout << "resultspred[" << i_cal << "] = " << resultspred[i_cal] << std::endl;
    //std::cout << "resultspredStatUp[i_cal] = " << resultspredStatUp[i_cal] << std::endl;
    //std::cout << "resultspredStatDn[i_cal] = " << resultspredStatDn[i_cal] << std::endl;
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatUp[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatUp[i_cal])*(resultspred[i_cal]-resultspredStatUp[i_cal]) + (resultspred[i_cal]-resultspredSystUp[i_cal])*(resultspred[i_cal]-resultspredSystUp[i_cal])) / resultspred[i_cal] << " ";

    else std::cout << "nan ";
  }
  std::cout << std::endl;

  std::cout << "mu iso syst down:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatDn[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatDn[i_cal])*(resultspred[i_cal]-resultspredStatDn[i_cal]) + (resultspred[i_cal]-resultspredSystDn[i_cal])*(resultspred[i_cal]-resultspredSystDn[i_cal])) / resultspred[i_cal] << " ";

    else std::cout << "nan ";
  }
  std::cout << std::endl;


  // restore default:
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_mus_isoeff[i][j]=MuRecoEff[i][j];
    }
  }


  }

  if (dodimu)
  {
    std::cout << "Computing dimuon syst:" << std::endl;
    double resultspreddimuSyst[NSEARCH_BINS]={0};
    const double old_ttbar_corrfactor_di_mus=ttbar_corrfactor_di_mus;
    ttbar_corrfactor_di_mus=1.0-1.5*(1.0-ttbar_corrfactor_di_mus);
    if (ttbar_corrfactor_di_mus>1.0) ttbar_corrfactor_di_mus=1.0;
    if (ttbar_corrfactor_di_mus<0.0) ttbar_corrfactor_di_mus=0.0;
    LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspreddimuSyst );

    myAccRecoIsoEffs.init();

    std::cout << "Computing dimuon stat:" << std::endl;
    double resultspreddimuStat[NSEARCH_BINS]={0};
    ttbar_corrfactor_di_mus=old_ttbar_corrfactor_di_mus+0.0062;
    if (ttbar_corrfactor_di_mus>1.0) ttbar_corrfactor_di_mus=1.0;
    if (ttbar_corrfactor_di_mus<0.0) ttbar_corrfactor_di_mus=0.0;
    LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspreddimuStat );

    myAccRecoIsoEffs.init();

    std::cout << "dimuon syst:" << std::endl;
    for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
    {
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspreddimuSyst[i_cal])*(resultspred[i_cal]-resultspreddimuSyst[i_cal]) + (resultspred[i_cal]-resultspreddimuStat[i_cal])*(resultspred[i_cal]-resultspreddimuStat[i_cal])) / resultspred[i_cal] << " ";

    else std::cout << "nan ";
    }
    std::cout << std::endl;

    // restore default:
    ttbar_corrfactor_di_mus=old_ttbar_corrfactor_di_mus;
  }

  if (dodie)
  {
    std::cout << "Computing diel syst:" << std::endl;
    double resultspreddieSyst[NSEARCH_BINS]={0};
    const double old_ttbar_corrfactor_di_els=ttbar_corrfactor_di_els;
    ttbar_corrfactor_di_els=1.0-1.5*(1.0-ttbar_corrfactor_di_els);
    if (ttbar_corrfactor_di_els>1.0) ttbar_corrfactor_di_els=1.0;
    if (ttbar_corrfactor_di_els<0.0) ttbar_corrfactor_di_els=0.0;
    LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspreddieSyst );

    myAccRecoIsoEffs.init();

    std::cout << "Computing diel stat:" << std::endl;
    double resultspreddieStat[NSEARCH_BINS]={0};
    ttbar_corrfactor_di_els=old_ttbar_corrfactor_di_els+0.00847;
    if (ttbar_corrfactor_di_els>1.0) ttbar_corrfactor_di_els=1.0;
    if (ttbar_corrfactor_di_els<0.0) ttbar_corrfactor_di_els=0.0;
    LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspreddieStat );

    myAccRecoIsoEffs.init();

    std::cout << "die syst:" << std::endl;
    for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
    {
      //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspreddieSyst[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspreddieStat[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspreddieSyst[i_cal])*(resultspred[i_cal]-resultspreddieSyst[i_cal]) + (resultspred[i_cal]-resultspreddieStat[i_cal])*(resultspred[i_cal]-resultspreddieStat[i_cal])) / resultspred[i_cal] << " ";
      
    else std::cout << "nan ";
    }
    std::cout << std::endl;

    // restore default:
    ttbar_corrfactor_di_els=old_ttbar_corrfactor_di_els;
  }


  if (doereco)
  {
  std::cout << "Computing syst el reco stat up: " << std::endl;
  double resultspredStatUp[NSEARCH_BINS]={0};
  double MuRecoEff[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_up[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_dn[PT_BINS][AC_BINS] = {{0}};
  TFile *fin = TFile::Open("v2_Effs2dPlots.root");
  TH2D * murecoeff;
  murecoeff = (TH2D*)fin->Get("els_recoeffs")->Clone();

/*  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //std::cout << "murecoeff->GetBinContent(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinContent(i+1,j+2) << std::endl;
      //std::cout << "murecoeff->GetBinError(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinError(i+1,j+2) << std::endl;
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      MuRecoEff_Stat_Unc_up[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]>1.0) MuRecoEff_Stat_Unc_up[i][j]=1.0;
      MuRecoEff_Stat_Unc_dn[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]<0.0) MuRecoEff_Stat_Unc_up[i][j]=0.0;
      //std::cout << "MuRecoEff[" << i << "][" << j << "] = " << MuRecoEff[i][j] << std::endl;
      //std::cout << "ttbar_mus_isoeff = " << ttbar_mus_isoeff[i][j] << std::endl;
    }
  }
*/
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      ttbar_els_recoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (ttbar_els_recoeff[i][j] > 1.0) ttbar_els_recoeff[i][j] = 1.0;
      //ttbar_els_recoeff[i][j]=MuRecoEff_Stat_Unc_up[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el reco stat dn: " << std::endl;
  double resultspredStatDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_els_recoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (ttbar_els_recoeff[i][j] < 0.0) {ttbar_els_recoeff[i][j] = 0.0;
      std::cout << "el reco eff [" << i << "][" << j << "] set to 0" << std::endl;}
      //ttbar_els_recoeff[i][j]=MuRecoEff_Stat_Unc_dn[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatDn );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el reco syst up: " << std::endl;
  double resultspredSystUp[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //ttbar_els_recoeff[i][j]=MuRecoEff[i][j]+0.05;
      if (i==0) ttbar_els_recoeff[i][j]=MuRecoEff[i][j]+0.055;
      else ttbar_els_recoeff[i][j]=MuRecoEff[i][j]+0.025;
      if (ttbar_els_recoeff[i][j]>1.0) ttbar_els_recoeff[i][j]=1.0;
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el reco syst dn: " << std::endl;
  double resultspredSystDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //ttbar_els_recoeff[i][j]=MuRecoEff[i][j]-0.05;
      if (i==0) ttbar_els_recoeff[i][j]=MuRecoEff[i][j]-0.055;
      else ttbar_els_recoeff[i][j]=MuRecoEff[i][j]-0.025;
      if (ttbar_els_recoeff[i][j]<0.0) {ttbar_els_recoeff[i][j]=0.0;
      std::cout << "el reco eff [" << i << "][" << j << "] set to 0" << std::endl;}
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystDn );

  myAccRecoIsoEffs.init();

  std::cout << "el reco syst up:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //std::cout << "resultspred[" << i_cal << "] = " << resultspred[i_cal] << std::endl;
    //std::cout << "resultspredStatUp[i_cal] = " << resultspredStatUp[i_cal] << std::endl;
    //std::cout << "resultspredStatDn[i_cal] = " << resultspredStatDn[i_cal] << std::endl;
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatUp[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatUp[i_cal])*(resultspred[i_cal]-resultspredStatUp[i_cal]) + (resultspred[i_cal]-resultspredSystUp[i_cal])*(resultspred[i_cal]-resultspredSystUp[i_cal])) / resultspred[i_cal] << " ";

    else std::cout << "nan ";
  }
  std::cout << std::endl;

  std::cout << "el reco syst down:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatDn[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatDn[i_cal])*(resultspred[i_cal]-resultspredStatDn[i_cal]) + (resultspred[i_cal]-resultspredSystDn[i_cal])*(resultspred[i_cal]-resultspredSystDn[i_cal])) / resultspred[i_cal] << " ";
    else std::cout << "nan ";
  }
  std::cout << std::endl;


  // restore default:
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_els_recoeff[i][j]=MuRecoEff[i][j];
    }
  }

  }

  if (doeiso)
  {
  std::cout << "Computing syst el iso stat up: " << std::endl;
  double resultspredStatUp[NSEARCH_BINS]={0};
  double MuRecoEff[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_up[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Stat_Unc_dn[PT_BINS][AC_BINS] = {{0}};
  TFile *fin = TFile::Open("v2_Effs2dPlots.root");
  TH2D * murecoeff;
  murecoeff = (TH2D*)fin->Get("els_isoeffs")->Clone();

/*  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //std::cout << "murecoeff->GetBinContent(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinContent(i+1,j+2) << std::endl;
      //std::cout << "murecoeff->GetBinError(" << i+1 << "," << j+2 << ") = " << murecoeff->GetBinError(i+1,j+2) << std::endl;
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      MuRecoEff_Stat_Unc_up[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]>1.0) MuRecoEff_Stat_Unc_up[i][j]=1.0;
      MuRecoEff_Stat_Unc_dn[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (MuRecoEff_Stat_Unc_up[i][j]<0.0) MuRecoEff_Stat_Unc_up[i][j]=0.0;
      //std::cout << "MuRecoEff[" << i << "][" << j << "] = " << MuRecoEff[i][j] << std::endl;
      //std::cout << "ttbar_mus_isoeff = " << ttbar_mus_isoeff[i][j] << std::endl;
    }
  }
*/
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      MuRecoEff[i][j] = murecoeff->GetBinContent(i+1,j+2);
      ttbar_els_isoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)+murecoeff->GetBinError(i+1,j+2);
      if (ttbar_els_isoeff[i][j] > 1.0) ttbar_mus_isoeff[i][j] = 1.0; 
      //ttbar_els_isoeff[i][j]=MuRecoEff_Stat_Unc_up[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el iso stat dn: " << std::endl;
  double resultspredStatDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_els_isoeff[i][j] = murecoeff->GetBinContent(i+1,j+2)-murecoeff->GetBinError(i+1,j+2);
      if (ttbar_els_isoeff[i][j] < 0.0) {ttbar_els_isoeff[i][j] = 0.0;
      std::cout << "el iso eff [" << i << "][" << j << "] set to 0" << std::endl;}
      //ttbar_els_isoeff[i][j]=MuRecoEff_Stat_Unc_dn[i][j];
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredStatDn );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el iso syst up: " << std::endl;
  double resultspredSystUp[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //ttbar_els_isoeff[i][j]=MuRecoEff[i][j]+0.04;
      if (i==0) ttbar_els_isoeff[i][j]=MuRecoEff[i][j]+0.035;
      else ttbar_els_isoeff[i][j]=MuRecoEff[i][j]+0.02;
      if (ttbar_els_isoeff[i][j]>1.0) ttbar_els_isoeff[i][j]=1.0;
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystUp );

  myAccRecoIsoEffs.init();

  std::cout << "Computing syst el iso syst dn: " << std::endl;
  double resultspredSystDn[NSEARCH_BINS]={0};
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      //ttbar_els_isoeff[i][j]=MuRecoEff[i][j]-0.04;
      if (i==0) ttbar_els_isoeff[i][j]=MuRecoEff[i][j]-0.035;
      else ttbar_els_isoeff[i][j]=MuRecoEff[i][j]-0.02;
      if (ttbar_els_isoeff[i][j]<0.0) {ttbar_els_isoeff[i][j]=0.0;
      std::cout << "el iso eff [" << i << "][" << j << "] set to 0" << std::endl;}
    }
  }
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, resultspredSystDn );

  myAccRecoIsoEffs.init();

  std::cout << "el iso syst up:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //std::cout << "resultspred[" << i_cal << "] = " << resultspred[i_cal] << std::endl;
    //std::cout << "resultspredStatUp[i_cal] = " << resultspredStatUp[i_cal] << std::endl;
    //std::cout << "resultspredStatDn[i_cal] = " << resultspredStatDn[i_cal] << std::endl;
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatUp[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatUp[i_cal])*(resultspred[i_cal]-resultspredStatUp[i_cal]) + (resultspred[i_cal]-resultspredSystUp[i_cal])*(resultspred[i_cal]-resultspredSystUp[i_cal])) / resultspred[i_cal] << " ";
    
    else std::cout << "nan ";
  }
  std::cout << std::endl;

  std::cout << "el iso syst down:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatDn[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << sqrt((resultspred[i_cal]-resultspredStatDn[i_cal])*(resultspred[i_cal]-resultspredStatDn[i_cal]) + (resultspred[i_cal]-resultspredSystDn[i_cal])*(resultspred[i_cal]-resultspredSystDn[i_cal])) / resultspred[i_cal] << " ";
    
    else std::cout << "nan ";
  }
  std::cout << std::endl;


  // restore default:
  for(int i=0;i<PT_BINS;++i)
  {
    for(int j=0;j<AC_BINS;++j)
    {
      ttbar_els_isoeff[i][j]=MuRecoEff[i][j];
    }
  }
  }
/*
  if (doacc)
  {
  std::cout << "Computing syst acc stat up: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs14;
  double resultspredStatUp[NSEARCH_BINS]={0};
  double ttbar_mus_acc6[59] = {0.833423,0.788406,0.842337,0.863853,0.841267,0.80138,0.832275,0.832,0.923959,0.80575,0.811259,0.802173,0.829273,0.794774,0.828258,0.867018,0.843241,0.794913,0.715186,0.746076,0.7909,0.725056,0.838593,0.832711,0.832305,0.805965,0.801153,0.779009,0.897866,0.915114,0.84947,0.972856,0.816392,0.838204,0.904901,0.991136,0.918113,0.845345,0.811071,0.875981,0.907228,0.924226,0.934703,0.833212,0.803378,0.783213,0.874606,0.780441,0.847539,0.894155,0.930694,0.805498,0.794395,0.890993,0.932324,0.932169,0.902498,0.907606,1};
  double ttbar_els_acc6[59] = {0.839878,0.793378,0.818505,0.785364,0.829303,0.800936,0.822807,0.822,0.880977,0.796191,0.800274,0.858861,0.829924,0.798601,0.805243,0.805725,0.835147,0.818903,0.71145,0.868367,0.826053,0.660517,0.781286,0.825598,0.825683,0.877196,0.787516,0.737691,0.893551,0.890074,0.937998,0.991681,0.822719,0.835079,0.90463,0.978372,0.885118,0.860804,0.778661,0.883716,0.906283,0.904205,0.879994,0.81258,0.863363,0.896459,0.821712,0.77822,0.831752,0.889058,0.940685,0.831656,0.84005,0.945315,0.897005,0.869458,0.908107,0.935211,1};

  double err_acc[59];
  double err_acc_el[59];
err_acc[0] = 0.0025547787264;
err_acc[1] = 0.008044880736;
err_acc[2] = 0.017677938233;
err_acc[3] = 0.032967020858;
err_acc[4] = 0.007943290363;
err_acc[5] = 0.011614180749;
err_acc[6] = 0.021985405468;
err_acc[7] = 0.041319350952;
err_acc[8] = 0.027238687863;
err_acc[9] = 0.014787044296;
err_acc[10] = 0.01955513133;
err_acc[11] = 0.01782642436;
err_acc[12] = 0.0027855075918;
err_acc[13] = 0.0091998764441;
err_acc[14] = 0.026730719058;
err_acc[15] = 0.080920024083;
err_acc[16] = 0.012990520715;
err_acc[17] = 0.01952869548;
err_acc[18] = 0.044332931326;
err_acc[19] = 0.068816506099;
err_acc[20] = 0.025345713374;
err_acc[21] = 0.042801489792;
err_acc[22] = 0.038468552577;
err_acc[23] = 0.0066505224419;
err_acc[24] = 0.023785038319;
err_acc[25] = 0.067072614177;
err_acc[26] = 0.028757498162;
err_acc[27] = 0.029120258581;
err_acc[28] = 0.0031411003831;
err_acc[29] = 0.010078729336;
err_acc[30] = 0.03595792246;
err_acc[31] = 0.079532454142;
err_acc[32] = 0.011331035558;
err_acc[33] = 0.018389150254;
err_acc[34] = 0.039969758508;
err_acc[35] = 0.032856460033;
err_acc[36] = 0.029617978961;
err_acc[37] = 0.025635597184;
err_acc[38] = 0.037359174312;
err_acc[39] = 0.028355352218;
err_acc[40] = 0.003425895945;
err_acc[41] = 0.011261810654;
err_acc[42] = 0.03553305388;
err_acc[43] = 0.01390096196;
err_acc[44] = 0.027876603582;
err_acc[45] = 0.063412438992;
err_acc[46] = 0.068479028757;
err_acc[47] = 0.035309643975;
err_acc[48] = 0.0411564164;
err_acc[49] = 0.009412867816;
err_acc[50] = 0.026641888211;
err_acc[51] = 0.031392126842;
err_acc[52] = 0.046586510551;
err_acc[53] = 0.017827584712;
err_acc[54] = 0.027648442108;
err_acc[55] = 0.016775709933;
err_acc[56] = 0.044782537467;
err_acc[57] = 0.039130720792;
err_acc[58] = 0.045228491032;
err_acc_el[0] = 0.002485511178;
err_acc_el[1] = 0.0078358401871;
err_acc_el[2] = 0.018233420426;
err_acc_el[3] = 0.03177578297;
err_acc_el[4] = 0.0078941390535;
err_acc_el[5] = 0.011369841686;
err_acc_el[6] = 0.025335140632;
err_acc_el[7] = 0.033970932291;
err_acc_el[8] = 0.024355554371;
err_acc_el[9] = 0.016643762484;
err_acc_el[10] = 0.01875294207;
err_acc_el[11] = 0.017075742628;
err_acc_el[12] = 0.002788411153;
err_acc_el[13] = 0.0095717147149;
err_acc_el[14] = 0.028095874717;
err_acc_el[15] = 0.06861922506;
err_acc_el[16] = 0.013180822037;
err_acc_el[17] = 0.018907106966;
err_acc_el[18] = 0.051884729241;
err_acc_el[19] = 0.071389610488;
err_acc_el[20] = 0.02887668036;
err_acc_el[21] = 0.045808632739;
err_acc_el[22] = 0.038139567971;
err_acc_el[23] = 0.006697210574;
err_acc_el[24] = 0.023125787036;
err_acc_el[25] = 0.058222980551;
err_acc_el[26] = 0.030640423393;
err_acc_el[27] = 0.033487027127;
err_acc_el[28] = 0.0031628838338;
err_acc_el[29] = 0.011072717661;
err_acc_el[30] = 0.022711278805;
err_acc_el[31] = 0.064188256475;
err_acc_el[32] = 0.011120003472;
err_acc_el[33] = 0.017447980246;
err_acc_el[34] = 0.034255979511;
err_acc_el[35] = 0.052636967484;
err_acc_el[36] = 0.035987969874;
err_acc_el[37] = 0.022609518057;
err_acc_el[38] = 0.034106727801;
err_acc_el[39] = 0.030389761895;
err_acc_el[40] = 0.003376580251;
err_acc_el[41] = 0.01349109357;
err_acc_el[42] = 0.041154376948;
err_acc_el[43] = 0.014411567787;
err_acc_el[44] = 0.023919243717;
err_acc_el[45] = 0.050455032914;
err_acc_el[46] = 0.054782191295;
err_acc_el[47] = 0.040871026104;
err_acc_el[48] = 0.042801489792;
err_acc_el[49] = 0.009270851512;
err_acc_el[50] = 0.027803654334;
err_acc_el[51] = 0.031668088252;
err_acc_el[52] = 0.046573380602;
err_acc_el[53] = 0.016760293495;
err_acc_el[54] = 0.026729429441;
err_acc_el[55] = 0.021725484712;
err_acc_el[56] = 0.057163206403;
err_acc_el[57] = 0.03673688581;
err_acc_el[58] = 0.1356854731;

  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc6[binc]+err_acc[binc];
    ttbar_els_acc[binc]=ttbar_els_acc6[binc]+err_acc_el[binc];
    if (ttbar_mus_acc[binc]>1.0) ttbar_mus_acc[binc]=1.0;
    if (ttbar_els_acc[binc]>1.0) ttbar_els_acc[binc]=1.0;
  }
  LoopLLPred( myAccRecoIsoEffs14, myTTJetsSampleWeight, resultspredStatUp );

  std::cout << "Computing syst acc stat dn: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs15;
  double resultspredStatDn[NSEARCH_BINS]={0};
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc6[binc]-err_acc[binc];
    ttbar_els_acc[binc]=ttbar_els_acc6[binc]-err_acc_el[binc];
    if (ttbar_mus_acc[binc]<0.0) ttbar_mus_acc[binc]=0.0;
    if (ttbar_els_acc[binc]<0.0) ttbar_els_acc[binc]=0.0;
  }
  LoopLLPred( myAccRecoIsoEffs15, myTTJetsSampleWeight, resultspredStatDn );

  std::cout << "Computing syst acc pdf up: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs16;
  double resultspredSystUp[NSEARCH_BINS]={0};
  double ttbar_mus_acc2[59] = {0.835942,0.832299,0.797418,0.868886,0.853178,0.806223,0.850961,0.439276,0.947606,0.84147,0.797776,0.778086,0.831727,0.786137,0.845204,0.866681,0.844801,0.788654,0.711573,0.718395,0.788081,0.78812,0.89033,0.834579,0.826702,0.712294,0.804286,0.772203,0.903474,0.919929,0.840848,1,0.804077,0.834124,0.926799,1,0.907157,0.866454,0.836939,0.859473,0.908054,0.928593,0.930358,0.837006,0.767007,0.779571,0.884552,0.776942,0.844341,0.895383,0.673334,0.808444,0.800083,1.23523,0.933004,0.92923,0.903206,0.905268,1};
  double ttbar_els_acc2[59] = {0.842308,0.789583,0.819035,0.847989,0.82673,0.807729,0.781456,0.897938,0.855211,0.771078,0.762294,0.912856,0.829896,0.802186,0.802844,0.858593,0.83824,0.810224,0.707515,0.865934,0.806935,0.684365,0.78397,0.828735,0.826183,0.881464,0.784563,0.733223,0.897106,0.843994,0.934874,1,0.818133,0.826258,0.901842,1,0.898618,0.859227,0.777391,0.911587,0.904924,0.904533,0.885713,0.815396,0.864425,0.889098,0.834277,0.762281,0.824304,0.891904,0.94314,0.832135,0.848994,0.947651,0.893544,0.838926,0.920293,0.935698,1};
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc2[binc];
    ttbar_els_acc[binc]=ttbar_els_acc2[binc];
  }
  LoopLLPred( myAccRecoIsoEffs16, myTTJetsSampleWeight, resultspredSystUp );

  std::cout << "Computing syst acc pdf dn: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs17;
  double resultspredSystDn[NSEARCH_BINS]={0};
  double ttbar_mus_acc3[59] = {0.842158,0.707058,0.868398,0.868921,0.856551,0.80476,0.848929,0.448837,0.945858,0.818196,0.803146,0.756212,0.829734,0.80942,0.840484,0.874177,0.843975,0.786238,0.706613,0.718837,0.825595,0.787384,0.866403,0.830943,0.847074,0.930775,0.800331,0.767626,0.900004,0.91957,0.839868,1,0.84262,0.826749,0.922017,1,0.931568,0.85115,0.821541,0.848096,0.90873,0.9266,0.940077,0.837949,0.838783,0.769485,0.881219,0.784547,0.835651,0.895017,1.65711,0.810012,0.80732,0.678529,0.932089,0.939914,0.900571,0.907857,1};
  double ttbar_els_acc3[59] = {0.845022,0.789575,0.818076,0.71113,0.825975,0.807848,0.85276,0.894062,0.84008,0.771955,0.762126,0.908001,0.830387,0.800688,0.805924,0.81561,0.831091,0.806336,0.68826,0.862915,0.807593,0.671919,0.77522,0.826196,0.827217,0.877537,0.787408,0.733254,0.89288,0.917609,0.930889,1,0.81867,0.825223,0.89031,1,0.870045,0.85744,0.778499,0.91687,0.907382,0.903446,0.86509,0.816321,0.86142,0.89535,0.79672,0.768529,0.859507,0.88588,0.939136,0.830355,0.861618,0.944825,0.890869,0.903647,0.91362,0.934055,1};
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc3[binc];
    ttbar_els_acc[binc]=ttbar_els_acc3[binc];
  }
  LoopLLPred( myAccRecoIsoEffs17, myTTJetsSampleWeight, resultspredSystDn );

  std::cout << "Computing syst acc scale up: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs18;
  double resultspredSystUp2[NSEARCH_BINS]={0};
  double ttbar_mus_acc4[59] = {0.837651,0.790679,0.83969,0.870514,0.853486,0.804352,0.851566,0.430044,0.947874,0.831368,0.795856,0.779392,0.830218,0.796942,0.84461,0.8686,0.84614,0.787073,0.712541,0.724837,0.805657,0.782578,0.874983,0.831889,0.832903,0.790081,0.80297,0.772253,0.90127,0.919773,0.840714,1,0.821116,0.827146,0.92686,1,0.915553,0.858658,0.826881,0.855567,0.907353,0.926684,0.935511,0.836189,0.801333,0.783826,0.89065,0.774666,0.843358,0.893828,0.93142,0.810167,0.803818,0.890293,0.931763,0.933613,0.900331,0.897961,1};
  double ttbar_els_acc4[59] = {0.843227,0.789566,0.815277,0.804687,0.82565,0.807416,0.823379,0.884267,0.846802,0.772425,0.761498,0.908631,0.829896,0.802408,0.804154,0.843299,0.834023,0.809443,0.698945,0.863544,0.807888,0.67643,0.789905,0.827156,0.825616,0.885908,0.786376,0.730801,0.894148,0.888048,0.932228,1,0.818072,0.826643,0.89923,1,0.885802,0.857075,0.773192,0.916377,0.905575,0.903509,0.874765,0.81434,0.862877,0.895474,0.814872,0.766086,0.839022,0.888598,0.942568,0.83131,0.84892,0.944912,0.889125,0.867756,0.911444,0.937204,1};
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc4[binc];
    ttbar_els_acc[binc]=ttbar_els_acc4[binc];
  }
  LoopLLPred( myAccRecoIsoEffs18, myTTJetsSampleWeight, resultspredSystUp2 );

  std::cout << "Computing syst acc scale dn: " << std::endl;
  AccRecoIsoEffs myAccRecoIsoEffs19;
  double resultspredSystDn2[NSEARCH_BINS]={0};
  double ttbar_mus_acc5[59] = {0.839897,0.791441,0.837379,0.867621,0.855909,0.806546,0.848504,0.455962,0.945911,0.831714,0.804106,0.757851,0.831205,0.796551,0.841516,0.871744,0.842921,0.787879,0.70635,0.713136,0.803235,0.791712,0.884969,0.833666,0.838971,0.810113,0.802199,0.768245,0.90233,0.919809,0.839996,1,0.822364,0.833805,0.922565,1,0.920508,0.860255,0.83289,0.853232,0.909263,0.928451,0.934503,0.838329,0.797851,0.767084,0.876102,0.78564,0.837882,0.896335,0.932821,0.808221,0.803317,0.889662,0.933194,0.934768,0.90377,0.913489,1};
  double ttbar_els_acc5[59] = {0.843927,0.789582,0.821337,0.81532,0.826926,0.807982,0.822209,0.906132,0.849997,0.770664,0.762946,0.912124,0.830336,0.800632,0.804152,0.838445,0.835561,0.807472,0.698311,0.865534,0.806778,0.680833,0.770778,0.827814,0.82751,0.873795,0.785606,0.735266,0.895934,0.888802,0.933777,1,0.818556,0.825061,0.894416,1,0.887235,0.859513,0.782035,0.912185,0.906534,0.904579,0.878267,0.817158,0.863137,0.889022,0.821077,0.764775,0.840988,0.889514,0.940238,0.831354,0.859698,0.947539,0.894934,0.869582,0.922554,0.933024,1};
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc5[binc];
    ttbar_els_acc[binc]=ttbar_els_acc5[binc];
  }
  LoopLLPred( myAccRecoIsoEffs19, myTTJetsSampleWeight, resultspredSystDn2 );

  std::cout << "e acc syst dn:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //std::cout << "resultspred[" << i_cal << "] = " << resultspred[i_cal] << std::endl;
    //std::cout << "resultspredStatUp[i_cal] = " << resultspredStatUp[i_cal] << std::endl;
    //std::cout << "resultspredStatDn[i_cal] = " << resultspredStatDn[i_cal] << std::endl;
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatUp[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredSystUp[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << " " << std::sqrt(std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredStatUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp2[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystUp2[i_cal])/resultspred[i_cal]);
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::sqrt(std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystUp[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystUp2[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystUp2[i_cal])/resultspred[i_cal]);
    //else std::cout << " nan";
  }
  std::cout << std::endl;

  std::cout << "e acc syst up:" << std::endl;
  for( int i_cal = 0 ; i_cal < NSEARCH_BINS ; i_cal++ )
  {
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredStatDn[i_cal];
    //if (resultspred[i_cal]>0.0) std::cout << " " << resultspred[i_cal]-resultspredSystDn[i_cal];
    if (resultspred[i_cal]>0.0) std::cout << " " << std::sqrt(std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredStatDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn2[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystDn2[i_cal])/resultspred[i_cal]);
    //if (resultspred[i_cal]>0.0) std::cout << " " << std::sqrt(std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystDn[i_cal])/resultspred[i_cal]+std::abs(resultspred[i_cal]-resultspredSystDn2[i_cal])/resultspred[i_cal]*std::abs(resultspred[i_cal]-resultspredSystDn2[i_cal])/resultspred[i_cal]);
    //else std::cout << " nan";
  }
  std::cout << std::endl;

  // restore default:
  for (int binc=0;binc<59;++binc)
  {
    ttbar_mus_acc[binc]=ttbar_mus_acc6[binc];
    ttbar_els_acc[binc]=ttbar_els_acc6[binc];
  }

  }
  

*/
  std::cout << "syst is done" << std::endl;
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
/*
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_HT-600to800", 2.666535, 14210872, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_HT-800to1200", 1.098082, 9982765, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_HT-1200to2500", 0.198748, 2932983, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_HT-2500toInf", 0.002368413, 1519815, LUMI, inputFileList_Cal );

  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromT_", 831.76*0.5*TTbar_SingleLept_BR, 61901450, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_SingleLeptFromTbar", 831.76*0.5*TTbar_SingleLept_BR, 59860282, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "TTJets_DiLept", 831.76*TTbar_DiLept_BR, 30444678, LUMI, inputFileList_Cal );

  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "tW_top_5f_inclusiveDecays" , 35.6, 6774350, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "tW_antitop_5f_inclusiveDecays" , 35.6, 6933094, LUMI, inputFileList_Cal );

  // 1.21 is the kf
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-200To400" , 359.7*1.21, 38867206, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-400To600" , 48.91*1.21, 7759701, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-600To800" , 12.05*1.21, 17494743, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-800To1200" , 5.501*1.21, 7745467, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-1200To2500" , 1.329*1.21, 6801534, LUMI, inputFileList_Cal );
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "WJetsToLNu_HT-2500ToInf" , 0.03216*1.21, 2637821, LUMI, inputFileList_Cal );
*/

  //data
  myTTJetsSampleWeight.TTJetsSampleInfo_push_back( "MET" , 1, 1, 1.0, inputFileList_Cal );  

  //LoopLLCal( myAccRecoIsoEffs, myTTJetsSampleWeight );
  //LoopLLExp( myAccRecoIsoEffs, myTTJetsSampleWeight );
  double results[NSEARCH_BINS]={0};
  LoopLLPred( myAccRecoIsoEffs, myTTJetsSampleWeight, results );
  //std::cout << " what the hell happened!!" << std::endl;
  //LoopLLSyst( myTTJetsSampleWeight );

  std::cout << "done" << std::endl;
  //std::cout << "main: printOverview" << std::endl;
  //myAccRecoIsoEffs.printOverview();
  //myAccRecoIsoEffs.NormalizeFlowNumber();
  // myAccRecoIsoEffs.printNormalizeFlowNumber();

  for(auto a: myTTJetsSampleWeight.TTJetsSampleInfos)
  {
	a.chain->Reset();
  }
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

  std::cout << "err_acc_mu[" << searchbinc << "] = " << get_ratio_Error(nmus_acc_sb[searchbinc], nmus_sb[searchbinc]-nmus_acc_sb[searchbinc], sqrt(nmus_acc_weight_squar[searchbinc]), sqrt(nmus_not_acc_weight_squar[searchbinc])) << ";" << std::endl;
    //std::cout << "err_acc[" << searchbinc << "] = " << get_stat_Error(nmus_acc_MC_sb[searchbinc],nmus_MC_sb[searchbinc]) << ";" << std::endl;
  }


  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {

   std::cout << "err_acc_el[" << searchbinc << "] = " << get_ratio_Error(nels_acc[searchbinc], nels[searchbinc]-nels_acc[searchbinc], sqrt(nels_acc_weight_squar[searchbinc]), sqrt(nels_not_acc_weight_squar[searchbinc])) << ";" << std::endl;
    //std::cout << "err_acc_el[" << searchbinc << "] = " << get_stat_Error(nels_acc_MC[searchbinc],nels_MC[searchbinc]) << ";" << std::endl;
  }


  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
/*      mus_recoeff[i_cal][j_cal]     = nmus_reco[i_cal][j_cal]/nmus_acc_bin[i_cal][j_cal];
      //here the error calculation is wrong?

      mus_recoeff_err[i_cal][j_cal] = get_stat_Error(nmus_reco_MC[i_cal][j_cal],nmus_acc_bin_MC[i_cal][j_cal]);
      els_recoeff[i_cal][j_cal]     = nels_reco[i_cal][j_cal]/nels_acc_bin[i_cal][j_cal];
      els_recoeff_err[i_cal][j_cal] = get_stat_Error(nels_reco_MC[i_cal][j_cal],nels_acc_bin_MC[i_cal][j_cal]);

      mus_isoeff[i_cal][j_cal]     = nmus_iso[i_cal][j_cal]/nmus_reco[i_cal][j_cal];
      mus_isoeff_err[i_cal][j_cal] = get_stat_Error(nmus_iso_MC[i_cal][j_cal],nmus_reco_MC[i_cal][j_cal]);
      els_isoeff[i_cal][j_cal]     = nels_iso[i_cal][j_cal]/nels_reco[i_cal][j_cal];
      els_isoeff_err[i_cal][j_cal] = get_stat_Error(nels_iso_MC[i_cal][j_cal],nels_reco_MC[i_cal][j_cal]);
*/
      mus_isoeff_allreco[i_cal][j_cal]     = nmus_iso_allreco[i_cal][j_cal]/nmus_reco_allreco[i_cal][j_cal];
      mus_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nmus_iso_MC_allreco[i_cal][j_cal],nmus_reco_MC_allreco[i_cal][j_cal]);
      els_isoeff_allreco[i_cal][j_cal]     = nels_iso_allreco[i_cal][j_cal]/nels_reco_allreco[i_cal][j_cal];
      els_isoeff_err_allreco[i_cal][j_cal] = get_stat_Error(nels_iso_MC_allreco[i_cal][j_cal],nels_reco_MC_allreco[i_cal][j_cal]);
    }
  }
  
  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < ETA_BINS ; j_cal++)
    {
      mus_recoeff[i_cal][j_cal]     = nmus_reco[i_cal][j_cal]/nmus_acc_bin[i_cal][j_cal];
      mus_recoeff_err[i_cal][j_cal] = get_stat_Error(nmus_reco_MC[i_cal][j_cal],nmus_acc_bin_MC[i_cal][j_cal]);
      els_recoeff[i_cal][j_cal]     = nels_reco[i_cal][j_cal]/nels_acc_bin[i_cal][j_cal];
      els_recoeff_err[i_cal][j_cal] = get_stat_Error(nels_reco_MC[i_cal][j_cal],nels_acc_bin_MC[i_cal][j_cal]);
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
TFile *Effs2dPlots = new TFile("Effs2dPlots.root", "recreate");
double ptbins[PT_BINS+1]={10.0,20.0,30.0,40.0,50.0,70.0,100.0,120.0};
double acbins[AC_BINS+2]={0.0,0.005,0.02,0.05,0.15,1.0,10.0};
double etabins[ETA_BINS+1]={0.0, 0.2, 0.3, 0.9, 1.2, 1.6, 2.1, 5.0};
double njetbins[NJETS_BINS+1]={3.5,4.5,5.5,6.5,7.5,8.5,9.5};
double htbins[NHT_BINS+1]={500.0,650.0,900.0,1400.0};

//TH2D *mus_recoeffs2d  = new TH2D("mus_recoeffs","Muon RecoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
TH2D *mus_recoeffs2d  = new TH2D("mus_recoeffs","Muon RecoEffs",PT_BINS,ptbins,ETA_BINS,etabins);
TH2D *mus_isoeffs2d  = new TH2D("mus_isoeffs","Muon IsoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);

//TH2D *els_recoeffs2d  = new TH2D("els_recoeffs","Electron RecoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);
TH2D *els_recoeffs2d  = new TH2D("els_recoeffs","Electron RecoEffs",PT_BINS,ptbins,ETA_BINS,etabins);
TH2D *els_isoeffs2d  = new TH2D("els_isoeffs","Electron IsoEffs",PT_BINS,ptbins,AC_BINS+1,acbins);

TH2D *mus_acc2d  = new TH2D("mus_acc","mus_acc",NJETS_BINS,njetbins,NHT_BINS,htbins);

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < ETA_BINS ; j_cal++)
    {    
      mus_recoeffs2d->SetBinContent( i_cal+1 , j_cal+1, mus_recoeff[i_cal][j_cal] );
      els_recoeffs2d->SetBinContent( i_cal+1 , j_cal+1, els_recoeff[i_cal][j_cal] );

      mus_recoeffs2d->SetBinError( i_cal+1 , j_cal+1, mus_recoeff_err[i_cal][j_cal] );
      els_recoeffs2d->SetBinError( i_cal+1 , j_cal+1, els_recoeff_err[i_cal][j_cal] );
    }     
  }

  for(i_cal = 0 ; i_cal < PT_BINS ; i_cal++)
  {
    for(j_cal = 0 ; j_cal < AC_BINS ; j_cal++)
    {
      mus_isoeffs2d->SetBinContent( i_cal+1 , j_cal+2, mus_isoeff_allreco[i_cal][j_cal] );
      els_isoeffs2d->SetBinContent( i_cal+1 , j_cal+2, els_isoeff_allreco[i_cal][j_cal] );

      mus_isoeffs2d->SetBinError( i_cal+1 , j_cal+2, mus_isoeff_err_allreco[i_cal][j_cal] );
      els_isoeffs2d->SetBinError( i_cal+1 , j_cal+2, els_isoeff_err_allreco[i_cal][j_cal] );
    }       
  }
  mus_recoeffs2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_recoeffs2d->GetYaxis()->SetTitle("Muon Eta");
  mus_isoeffs2d->GetXaxis()->SetTitle("Muon Pt [GeV]");
  mus_isoeffs2d->GetYaxis()->SetTitle("Activity");
  els_recoeffs2d->GetXaxis()->SetTitle("Electron Pt [GeV]");
  els_recoeffs2d->GetYaxis()->SetTitle("Electron Eta");
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
	for (int etabins=0;etabins<ETA_BINS;++etabins)
	{

	  if (use_muon_control_sample)
	  {
          mus_EventWeight_iso[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0 - ttbar_mus_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal];
          mus_EventWeight_reco[i_cal][j_cal][k_cal][etabins][searchbinc] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] );
          mus_EventWeight_acc[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_mus_acc[searchbinc])/ttbar_mus_acc[searchbinc] );

          els_EventWeight_acc[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * (1.0/ttbar_mus_recoeff[j_cal][k_cal]) * ( (1.0 - ttbar_els_acc[searchbinc])/ttbar_mus_acc[searchbinc] );
          els_EventWeight_reco[i_cal][j_cal][k_cal][etabins][searchbinc] = (1.0/ttbar_mus_isoeff[j_cal][k_cal]) * ( (1 - ttbar_els_recoeff[j_cal][k_cal])/ttbar_mus_recoeff[j_cal][k_cal] )* (ttbar_els_acc[searchbinc]/ttbar_mus_acc[searchbinc]);
          els_EventWeight_iso[i_cal][j_cal][k_cal][etabins][searchbinc]  = ( (1.0 - ttbar_els_isoeff[j_cal][k_cal])/ttbar_mus_isoeff[j_cal][k_cal] ) * (ttbar_els_recoeff[j_cal][k_cal]/ttbar_mus_recoeff[j_cal][k_cal])* (ttbar_els_acc[searchbinc]/ttbar_mus_acc[searchbinc]);
	  }

	  if (use_electron_control_sample)
	  {
	  mus_EventWeight_iso[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0 - ttbar_mus_isoeff[j_cal][k_cal])/ttbar_els_isoeff[j_cal][k_cal] * ttbar_mus_recoeff[j_cal][k_cal]/ttbar_els_recoeff[j_cal][k_cal] * ttbar_mus_acc[searchbinc]/ttbar_els_acc[searchbinc];
	  mus_EventWeight_reco[i_cal][j_cal][k_cal][etabins][searchbinc] = (1.0/ttbar_els_isoeff[j_cal][k_cal]) * (1.0 - ttbar_mus_recoeff[j_cal][k_cal])/ttbar_els_recoeff[j_cal][k_cal] * ttbar_mus_acc[searchbinc]/ttbar_els_acc[searchbinc];
	  mus_EventWeight_acc[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0/ttbar_els_isoeff[j_cal][k_cal]) * (1.0/ttbar_els_recoeff[j_cal][k_cal]) * (1.0 - ttbar_mus_acc[searchbinc])/ttbar_els_acc[searchbinc];

	  els_EventWeight_acc[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0/ttbar_els_isoeff[j_cal][k_cal]) * (1.0/ttbar_els_recoeff[j_cal][k_cal]) * (1.0 - ttbar_els_acc[searchbinc])/ttbar_els_acc[searchbinc];
	  els_EventWeight_reco[i_cal][j_cal][k_cal][etabins][searchbinc] = (1.0/ttbar_els_isoeff[j_cal][k_cal]) * (1.0 - ttbar_els_recoeff[j_cal][k_cal])/ttbar_els_recoeff[j_cal][k_cal];
	  els_EventWeight_iso[i_cal][j_cal][k_cal][etabins][searchbinc]  = (1.0 - ttbar_els_isoeff[j_cal][k_cal])/ttbar_els_isoeff[j_cal][k_cal];
	  }
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
  //std::cout.precision(3);

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

  EffsHeader << "  double ttbar_mtwcorrfactor[" << PT_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    if( i_cal == 0 ) { EffsHeader << "{"; }
    EffsHeader << mtwcorrfactor[i_cal];
    if( i_cal != PT_BINS-1 ) { EffsHeader << ","; }
    if( i_cal == PT_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }

  EffsHeader << "  double ttbar_mus_acc[" << NSEARCH_BINS << "] = "; 
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {

std::cout << "mu bin " << searchbinc << " acc " << nmus_acc_sb[searchbinc] << " all " << nmus_sb[searchbinc] << " eff " << nmus_acc_sb[searchbinc]/nmus_sb[searchbinc] << " raw acc " << nmus_acc_MC_sb[searchbinc] << " raw all " << nmus_MC_sb[searchbinc] << std::endl;
//std::cout << "eff" << nmus_acc_sb[searchbinc]/nmus_sb[searchbinc] << std::endl;

    if( searchbinc == 0 ) { EffsHeader << "{"; }
    EffsHeader << nmus_acc_sb[searchbinc]/nmus_sb[searchbinc];
    if( searchbinc != NSEARCH_BINS-1 ) { EffsHeader << ","; }
    if( searchbinc == NSEARCH_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }

//for merge mu acc bins!!!
/*
std::cout << "bin 4 = " << (nmus_acc_sb[4] + nmus_acc_sb[17]) / (nmus_sb[4] + nmus_sb[17]) << " err = " << get_stat_Error(nmus_acc_MC_sb[4] + nmus_acc_MC_sb[17] ,nmus_MC_sb[4] + nmus_MC_sb[17]) << std::endl;
std::cout << "bin 13 = " << (nmus_acc_sb[13] + nmus_acc_sb[14]) / (nmus_sb[13] + nmus_sb[14]) << " err = " << get_stat_Error(nmus_acc_MC_sb[13] + nmus_acc_MC_sb[14] ,nmus_MC_sb[13] + nmus_MC_sb[14]) << std::endl;
std::cout << "bin 17 = " << (nmus_acc_sb[4] + nmus_acc_sb[17]) / (nmus_sb[4] + nmus_sb[17]) << " err = " << get_stat_Error(nmus_acc_MC_sb[4] + nmus_acc_MC_sb[17] ,nmus_MC_sb[4] + nmus_MC_sb[17]) << std::endl;
std::cout << "bin 18 = " << (nmus_acc_sb[18] + nmus_acc_sb[19]) / (nmus_sb[18] + nmus_sb[19]) << " err = " << get_stat_Error(nmus_acc_MC_sb[18] + nmus_acc_MC_sb[19] ,nmus_MC_sb[18] + nmus_MC_sb[19]) << std::endl;
std::cout << "bin 19 = " << (nmus_acc_sb[18] + nmus_acc_sb[19]) / (nmus_sb[18] + nmus_sb[19]) << " err = " << get_stat_Error(nmus_acc_MC_sb[18] + nmus_acc_MC_sb[19] ,nmus_MC_sb[18] + nmus_MC_sb[19]) << std::endl;
std::cout << "bin 20 = " << (nmus_acc_sb[18] + nmus_acc_sb[19] + nmus_acc_sb[20]) / (nmus_sb[18] + nmus_sb[19] + nmus_sb[20]) << " err = " << get_stat_Error(nmus_acc_MC_sb[18] + nmus_acc_MC_sb[19]+ nmus_acc_MC_sb[20] ,nmus_MC_sb[18] + nmus_MC_sb[19] + nmus_MC_sb[20]) << std::endl;
std::cout << "bin 25 = " << (nmus_acc_sb[25] + nmus_acc_sb[36]) / (nmus_sb[25] + nmus_sb[36]) << " err = " << get_stat_Error(nmus_acc_MC_sb[25] + nmus_acc_MC_sb[36] ,nmus_MC_sb[25] + nmus_MC_sb[36]) << std::endl;
std::cout << "bin 29 = " << (nmus_acc_sb[29] + nmus_acc_sb[28]) / (nmus_sb[29] + nmus_sb[28]) << " err = " << get_stat_Error(nmus_acc_MC_sb[29] + nmus_acc_MC_sb[28] ,nmus_MC_sb[29] + nmus_MC_sb[28]) << std::endl;
std::cout << "bin 33 = " << (nmus_acc_sb[33] + nmus_acc_sb[32]) / (nmus_sb[33] + nmus_sb[32]) << " err = " << get_stat_Error(nmus_acc_MC_sb[33] + nmus_acc_MC_sb[32] ,nmus_MC_sb[33] + nmus_MC_sb[32]) << std::endl;
std::cout << "bin 34 = " << (nmus_acc_sb[34] + nmus_acc_sb[35]) / (nmus_sb[34] + nmus_sb[35]) << " err = " << get_stat_Error(nmus_acc_MC_sb[34] + nmus_acc_MC_sb[35] ,nmus_MC_sb[34] + nmus_MC_sb[35]) << std::endl;
std::cout << "bin 35 = " << (nmus_acc_sb[34] + nmus_acc_sb[35]) / (nmus_sb[34] + nmus_sb[35]) << " err = " << get_stat_Error(nmus_acc_MC_sb[34] + nmus_acc_MC_sb[35] ,nmus_MC_sb[34] + nmus_MC_sb[35]) << std::endl;
std::cout << "bin 36 = " << (nmus_acc_sb[25] + nmus_acc_sb[36]) / (nmus_sb[25] + nmus_sb[36]) << " err = " << get_stat_Error(nmus_acc_MC_sb[25] + nmus_acc_MC_sb[36] ,nmus_MC_sb[25] + nmus_MC_sb[36]) << std::endl;
std::cout << "bin 40 = " << (nmus_acc_sb[40] + nmus_acc_sb[39]) / (nmus_sb[40] + nmus_sb[39]) << " err = " << get_stat_Error(nmus_acc_MC_sb[40] + nmus_acc_MC_sb[39] ,nmus_MC_sb[40] + nmus_MC_sb[39]) << std::endl;
std::cout << "bin 47 = " << (nmus_acc_sb[47] + nmus_acc_sb[46]) / (nmus_sb[47] + nmus_sb[46]) << " err = " << get_stat_Error(nmus_acc_MC_sb[47] + nmus_acc_MC_sb[46] ,nmus_MC_sb[47] + nmus_MC_sb[46]) << std::endl;
std::cout << "bin 51 = " << (nmus_acc_sb[51] + nmus_acc_sb[57]) / (nmus_sb[51] + nmus_sb[57]) << " err = " << get_stat_Error(nmus_acc_MC_sb[51] + nmus_acc_MC_sb[57] ,nmus_MC_sb[51] + nmus_MC_sb[57]) << std::endl;
std::cout << "bin 57 = " << (nmus_acc_sb[51] + nmus_acc_sb[57]) / (nmus_sb[51] + nmus_sb[57]) << " err = " << get_stat_Error(nmus_acc_MC_sb[51] + nmus_acc_MC_sb[57] ,nmus_MC_sb[51] + nmus_MC_sb[57]) << std::endl;
std::cout << "bin 61 = " << (nmus_acc_sb[61] + nmus_acc_sb[67]) / (nmus_sb[61] + nmus_sb[67]) << " err = " << get_stat_Error(nmus_acc_MC_sb[61] + nmus_acc_MC_sb[67] ,nmus_MC_sb[61] + nmus_MC_sb[67]) << std::endl;
std::cout << "bin 67 = " << (nmus_acc_sb[61] + nmus_acc_sb[67]) / (nmus_sb[61] + nmus_sb[67]) << " err = " << get_stat_Error(nmus_acc_MC_sb[61] + nmus_acc_MC_sb[67] ,nmus_MC_sb[61] + nmus_MC_sb[67]) << std::endl;
std::cout << "bin 69 = " << (nmus_acc_sb[69] + nmus_acc_sb[68]) / (nmus_sb[69] + nmus_sb[68]) << " err = " << get_stat_Error(nmus_acc_MC_sb[69] + nmus_acc_MC_sb[68] ,nmus_MC_sb[69] + nmus_MC_sb[68]) << std::endl;
std::cout << "bin 72 = " << (nmus_acc_sb[72] + nmus_acc_sb[77]) / (nmus_sb[72] + nmus_sb[77]) << " err = " << get_stat_Error(nmus_acc_MC_sb[72] + nmus_acc_MC_sb[77] ,nmus_MC_sb[72] + nmus_MC_sb[77]) << std::endl;
std::cout << "bin 77 = " << (nmus_acc_sb[72] + nmus_acc_sb[77]) / (nmus_sb[72] + nmus_sb[77]) << " err = " << get_stat_Error(nmus_acc_MC_sb[72] + nmus_acc_MC_sb[77] ,nmus_MC_sb[72] + nmus_MC_sb[77]) << std::endl;
std::cout << "bin 79 = " << (nmus_acc_sb[79] + nmus_acc_sb[78]) / (nmus_sb[79] + nmus_sb[78]) << " err = " << get_stat_Error(nmus_acc_MC_sb[79] + nmus_acc_MC_sb[78] ,nmus_MC_sb[79] + nmus_MC_sb[78]) << std::endl;
std::cout << "bin 81 = " << (nmus_acc_sb[81] + nmus_acc_sb[80]) / (nmus_sb[81] + nmus_sb[80]) << " err = " << get_stat_Error(nmus_acc_MC_sb[81] + nmus_acc_MC_sb[80] ,nmus_MC_sb[81] + nmus_MC_sb[80]) << std::endl;
std::cout << "bin 82 = " << (nmus_acc_sb[82] + nmus_acc_sb[83]) / (nmus_sb[82] + nmus_sb[83]) << " err = " << get_stat_Error(nmus_acc_MC_sb[82] + nmus_acc_MC_sb[83] ,nmus_MC_sb[82] + nmus_MC_sb[83]) << std::endl;
std::cout << "bin 83 = " << (nmus_acc_sb[82] + nmus_acc_sb[83]) / (nmus_sb[82] + nmus_sb[83]) << " err = " << get_stat_Error(nmus_acc_MC_sb[82] + nmus_acc_MC_sb[83] ,nmus_MC_sb[82] + nmus_MC_sb[83]) << std::endl;
*/


  EffsHeader << "  double ttbar_mus_recoeff[" << PT_BINS << "][" << ETA_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < ETA_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << mus_recoeff[i_cal][j_cal];            
      if( j_cal != ETA_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == ETA_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == ETA_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }

  EffsHeader << "  double ttbar_mus_isoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
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

  EffsHeader << "  double ttbar_els_acc[" << NSEARCH_BINS << "] = ";                            
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {

std::cout << "el bin " << searchbinc << " acc " << nels_acc[searchbinc] << " all " << nels[searchbinc] << " eff " << nels_acc[searchbinc]/nels[searchbinc] << " raw acc " << nels_acc_MC[searchbinc] << " raw all " << nels_MC[searchbinc] << std::endl;

//els_acc[searchbinc] = nels_acc[searchbinc]/nels[searchbinc];
    if( searchbinc == 0 ) { EffsHeader << "{"; }
    //EffsHeader << els_acc[searchbinc];
    EffsHeader << nels_acc[searchbinc]/nels[searchbinc];
    if( searchbinc != NSEARCH_BINS-1 ) { EffsHeader << ","; }
    if( searchbinc == NSEARCH_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }
/*
std::cout << "bin 4 = " << (nels_acc[4] + nels_acc[17]) / (nels[4] + nels[17]) << " err = " << get_stat_Error(nels_acc_MC[4] + nels_acc_MC[17] ,nels_MC[4] + nels_MC[17]) << std::endl;
std::cout << "bin 13 = " << (nels_acc[13] + nels_acc[14]) / (nels[13] + nels[14]) << " err = " << get_stat_Error(nels_acc_MC[13] + nels_acc_MC[14] ,nels_MC[13] + nels_MC[14]) << std::endl;
std::cout << "bin 17 = " << (nels_acc[4] + nels_acc[17]) / (nels[4] + nels[17]) << " err = " << get_stat_Error(nels_acc_MC[4] + nels_acc_MC[17] ,nels_MC[4] + nels_MC[17]) << std::endl;
std::cout << "bin 18 = " << (nels_acc[18] + nels_acc[19] + nels_acc[20]) / (nels[18] + nels[19] + nels[20]) << " err = " << get_stat_Error(nels_acc_MC[18] + nels_acc_MC[19]+ nels_acc_MC[20] ,nels_MC[18] + nels_MC[19] + nels_MC[20]) << std::endl;
std::cout << "bin 19 = " << (nels_acc[18] + nels_acc[19] + nels_acc[20]) / (nels[18] + nels[19] + nels[20]) << " err = " << get_stat_Error(nels_acc_MC[18] + nels_acc_MC[19]+ nels_acc_MC[20] ,nels_MC[18] + nels_MC[19] + nels_MC[20]) << std::endl;
std::cout << "bin 20 = " << (nels_acc[18] + nels_acc[19] + nels_acc[20]) / (nels[18] + nels[19] + nels[20]) << " err = " << get_stat_Error(nels_acc_MC[18] + nels_acc_MC[19]+ nels_acc_MC[20] ,nels_MC[18] + nels_MC[19] + nels_MC[20]) << std::endl;
std::cout << "bin 25 = " << (nels_acc[25] + nels_acc[36]) / (nels[25] + nels[36]) << " err = " << get_stat_Error(nels_acc_MC[25] + nels_acc_MC[36] ,nels_MC[25] + nels_MC[36]) << std::endl;
std::cout << "bin 29 = " << (nels_acc[29] + nels_acc[28]) / (nels[29] + nels[28]) << " err = " << get_stat_Error(nels_acc_MC[29] + nels_acc_MC[28] ,nels_MC[29] + nels_MC[28]) << std::endl;
std::cout << "bin 33 = " << (nels_acc[33] + nels_acc[32]) / (nels[33] + nels[32]) << " err = " << get_stat_Error(nels_acc_MC[33] + nels_acc_MC[32] ,nels_MC[33] + nels_MC[32]) << std::endl;
std::cout << "bin 34 = " << (nels_acc[34] + nels_acc[35]) / (nels[34] + nels[35]) << " err = " << get_stat_Error(nels_acc_MC[34] + nels_acc_MC[35] ,nels_MC[34] + nels_MC[35]) << std::endl;
std::cout << "bin 35 = " << (nels_acc[34] + nels_acc[35]) / (nels[34] + nels[35]) << " err = " << get_stat_Error(nels_acc_MC[34] + nels_acc_MC[35] ,nels_MC[34] + nels_MC[35]) << std::endl;
std::cout << "bin 36 = " << (nels_acc[25] + nels_acc[36]) / (nels[25] + nels[36]) << " err = " << get_stat_Error(nels_acc_MC[25] + nels_acc_MC[36] ,nels_MC[25] + nels_MC[36]) << std::endl;
std::cout << "bin 40 = " << (nels_acc[40] + nels_acc[39]) / (nels[40] + nels[39]) << " err = " << get_stat_Error(nels_acc_MC[40] + nels_acc_MC[39] ,nels_MC[40] + nels_MC[39]) << std::endl;
std::cout << "bin 47 = " << (nels_acc[47] + nels_acc[46]) / (nels[47] + nels[46]) << " err = " << get_stat_Error(nels_acc_MC[47] + nels_acc_MC[46] ,nels_MC[47] + nels_MC[46]) << std::endl;
std::cout << "bin 51 = " << (nels_acc[51] + nels_acc[57]) / (nels[51] + nels[57]) << " err = " << get_stat_Error(nels_acc_MC[51] + nels_acc_MC[57] ,nels_MC[51] + nels_MC[57]) << std::endl;
std::cout << "bin 57 = " << (nels_acc[51] + nels_acc[57]) / (nels[51] + nels[57]) << " err = " << get_stat_Error(nels_acc_MC[51] + nels_acc_MC[57] ,nels_MC[51] + nels_MC[57]) << std::endl;
std::cout << "bin 61 = " << (nels_acc[61] + nels_acc[67]) / (nels[61] + nels[67]) << " err = " << get_stat_Error(nels_acc_MC[61] + nels_acc_MC[67] ,nels_MC[61] + nels_MC[67]) << std::endl;
std::cout << "bin 67 = " << (nels_acc[61] + nels_acc[67]) / (nels[61] + nels[67]) << " err = " << get_stat_Error(nels_acc_MC[61] + nels_acc_MC[67] ,nels_MC[61] + nels_MC[67]) << std::endl;
std::cout << "bin 69 = " << (nels_acc[69] + nels_acc[68]) / (nels[69] + nels[68]) << " err = " << get_stat_Error(nels_acc_MC[69] + nels_acc_MC[68] ,nels_MC[69] + nels_MC[68]) << std::endl;
std::cout << "bin 72 = " << (nels_acc[72] + nels_acc[77]) / (nels[72] + nels[77]) << " err = " << get_stat_Error(nels_acc_MC[72] + nels_acc_MC[77] ,nels_MC[72] + nels_MC[77]) << std::endl;
std::cout << "bin 77 = " << (nels_acc[72] + nels_acc[77]) / (nels[72] + nels[77]) << " err = " << get_stat_Error(nels_acc_MC[72] + nels_acc_MC[77] ,nels_MC[72] + nels_MC[77]) << std::endl;
std::cout << "bin 79 = " << (nels_acc[79] + nels_acc[78]) / (nels[79] + nels[78]) << " err = " << get_stat_Error(nels_acc_MC[79] + nels_acc_MC[78] ,nels_MC[79] + nels_MC[78]) << std::endl;
std::cout << "bin 81 = " << (nels_acc[81] + nels_acc[80]) / (nels[81] + nels[80]) << " err = " << get_stat_Error(nels_acc_MC[81] + nels_acc_MC[80] ,nels_MC[81] + nels_MC[80]) << std::endl;
std::cout << "bin 82 = " << (nels_acc[82] + nels_acc[83]) / (nels[82] + nels[83]) << " err = " << get_stat_Error(nels_acc_MC[82] + nels_acc_MC[83] ,nels_MC[82] + nels_MC[83]) << std::endl;
std::cout << "bin 83 = " << (nels_acc[82] + nels_acc[83]) / (nels[82] + nels[83]) << " err = " << get_stat_Error(nels_acc_MC[82] + nels_acc_MC[83] ,nels_MC[82] + nels_MC[83]) << std::endl;
*/
  EffsHeader << "  double ttbar_els_recoeff[" << PT_BINS << "][" << ETA_BINS << "] = ";
  for( i_cal = 0 ; i_cal < PT_BINS ; i_cal++ )
  {
    for( j_cal = 0 ; j_cal < ETA_BINS ; j_cal++ )
    {
      if( i_cal == 0 && j_cal == 0 ) { EffsHeader << "{{"; }
      if( i_cal != 0 && j_cal == 0 ) { EffsHeader << "{"; }

      EffsHeader << els_recoeff[i_cal][j_cal];
      if( j_cal != ETA_BINS-1 ) { EffsHeader << ","; }

      if( i_cal != PT_BINS-1 && j_cal == ETA_BINS-1 ) { EffsHeader << "},"; }
      if( i_cal == PT_BINS-1 && j_cal == ETA_BINS-1 ) { EffsHeader << "}};" << std::endl; }
    }
  }

  EffsHeader << "  double ttbar_els_isoeff[" << PT_BINS << "][" << AC_BINS << "] = ";
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


  EffsHeader << "  double ttbar_corrfactor_di_mus = " << corrfactor_di_mus << ";" << std::endl;
  EffsHeader << "  double ttbar_corrfactor_di_els = " << corrfactor_di_els << ";" << std::endl;

  EffsHeader << "  double isoTrackEff_SB[" << NSEARCH_BINS << "] = {";
  for( int searchbinc = 0 ; searchbinc < NSEARCH_BINS ; ++searchbinc )
  {
    EffsHeader << isotrkeff[searchbinc];
    if( searchbinc != NSEARCH_BINS-1 ) { EffsHeader << ","; }
    if( searchbinc == NSEARCH_BINS-1 ) { EffsHeader << "};" << std::endl; }
  }
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

/*std::cout << "pid=" << pid << std::endl;

std::cout << "genpfActivity=" << genpfActivity << std::endl;

for ( int i=0; i < IdFlag.size(); i++)
{
std::cout << "IdFlag=" << IdFlag.at(i) << std::endl;
}

for ( int i=0; i < recoleptactivityVec.size(); i++)
{
std::cout << "recoleptactivityVec=" << recoleptactivityVec.at(i) << std::endl;
}

for ( int i=0; i < MiniIso.size(); i++)
{
std::cout << "MiniIso=" << MiniIso.at(i) << std::endl;
}
*/

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
