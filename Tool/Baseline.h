#include "SusyAnaTools/Tools/NTupleReader.h"

const std::string spec = "lostlept";

void passBaselineFunc(NTupleReader &tr)
{
  bool debug = false;
  bool doIsoTrksVeto = true;
  bool doMuonVeto = true;
  bool doElectronVeto = true;  

  bool incZEROtop = false;

  bool passBaseline = true;
  bool passBaselineNoTag = true;

  std::string jetVecLabel = "jetsLVec";
  std::string CSVVecLabel = "recoJetsBtag_0";
  std::string METLabel    = "met";
  std::string METPhiLabel = "metphi";
 
  if( spec.compare("noIsoTrksVeto") == 0)
  {
     doIsoTrksVeto = false;
  }
  if( spec.compare("incZEROtop") == 0)
  {
     incZEROtop = true;
  }
  if( spec.compare("hadtau") == 0)
  {
     doMuonVeto = false;
  }
  if( spec.compare("lostlept") == 0)
  {
     doMuonVeto = false;
     doElectronVeto = false;
     doIsoTrksVeto = false;
  }
  if(spec.compare("Zinv") == 0) 
  {
    jetVecLabel = "cleanJetpt30ArrVec";//"jetsLVec";//"prodJetsNoMu_jetsLVec";
    CSVVecLabel = "cleanJetpt30ArrBTag";//"recoJetsBtag_0";
    METLabel    = "cleanMetPt";
    METPhiLabel = "cleanMetPhi";
    doMuonVeto  = false;
  }

  // Form TLorentzVector of MET
  TLorentzVector metLVec; metLVec.SetPtEtaPhiM(tr.getVar<double>(METLabel), 0, tr.getVar<double>(METPhiLabel), 0);

  // Calculate number of leptons
  int nMuons = AnaFunctions::countMuons(tr.getVec<TLorentzVector>("muonsLVec"), tr.getVec<double>("muonsMiniIso"), tr.getVec<double>("muonsMtw"), AnaConsts::muonsMiniIsoArr);
  int nElectrons = AnaFunctions::countElectrons(tr.getVec<TLorentzVector>("elesLVec"), tr.getVec<double>("elesMiniIso"), tr.getVec<double>("elesMtw"), tr.getVec<unsigned int>("elesisEB"), AnaConsts::elesMiniIsoArr);
  int nIsoTrks = AnaFunctions::countIsoTrks(tr.getVec<TLorentzVector>("loose_isoTrksLVec"), tr.getVec<double>("loose_isoTrks_iso"), tr.getVec<double>("loose_isoTrks_mtw"), tr.getVec<int>("loose_isoTrks_pdgId"));

  // Calculate number of jets and b-tagged jets
  int cntCSVS = AnaFunctions::countCSVS(tr.getVec<TLorentzVector>(jetVecLabel), tr.getVec<double>(CSVVecLabel), AnaConsts::cutCSVS, AnaConsts::bTagArr);
  int cntNJetsPt50Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt50Eta24Arr);
  int cntNJetsPt30Eta24 = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt30Eta24Arr);
        int cntNJetsPt30      = AnaFunctions::countJets(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt30Arr);

  // Calculate deltaPhi
  std::vector<double> * dPhiVec = new std::vector<double>();
  (*dPhiVec) = AnaFunctions::calcDPhi(tr.getVec<TLorentzVector>(jetVecLabel), metLVec.Phi(), 3, AnaConsts::dphiArr);

  // Prepare jets and b-tag working points for top tagger
  std::vector<TLorentzVector> *jetsLVec_forTagger = new std::vector<TLorentzVector>(); std::vector<double> *recoJetsBtag_forTagger = new std::vector<double>();
  AnaFunctions::prepareJetsForTagger(tr.getVec<TLorentzVector>(jetVecLabel), tr.getVec<double>(CSVVecLabel), (*jetsLVec_forTagger), (*recoJetsBtag_forTagger));
  if( debug ) std::cout<<"\njetsLVec_forTagger->size : "<<jetsLVec_forTagger->size()<<"  recoJetsBtag_forTagger->size : "<<recoJetsBtag_forTagger->size()<<"  passBaseline : "<<passBaseline<<std::endl;

  // Pass lepton veto?
  bool passLeptVeto = true, passMuonVeto = true, passEleVeto = true, passIsoTrkVeto = true;
  if( doMuonVeto && nMuons != AnaConsts::nMuonsSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passMuonVeto = false; }
  if( doElectronVeto && nElectrons != AnaConsts::nElectronsSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passEleVeto = false; }
  // Isolated track veto is disabled for now
  if( doIsoTrksVeto && nIsoTrks != AnaConsts::nIsoTrksSel ){ passBaseline = false; passBaselineNoTag = false; passLeptVeto = false; passIsoTrkVeto = false; }
  if( debug ) std::cout<<"nMuons : "<<nMuons<<"  nElectrons : "<<nElectrons<<"  nIsoTrks : "<<nIsoTrks<<"  passBaseline : "<<passBaseline<<std::endl;

  // Pass number of jets?
  bool passnJets = true;
  if( cntNJetsPt50Eta24 < AnaConsts::nJetsSelPt50Eta24 ){ passBaseline = false; passBaselineNoTag = false; passnJets = false; }
  if( cntNJetsPt30Eta24 < AnaConsts::nJetsSelPt30Eta24 ){ passBaseline = false; passBaselineNoTag = false; passnJets = false; }
  if( debug ) std::cout<<"cntNJetsPt50Eta24 : "<<cntNJetsPt50Eta24<<"  cntNJetsPt30Eta24 : "<<cntNJetsPt30Eta24<<"  cntNJetsPt30 : "<<cntNJetsPt30<<"  passBaseline : "<<passBaseline<<std::endl;

  // Pass deltaPhi?
  bool passdPhis = true;
  if( dPhiVec->at(0) < AnaConsts::dPhi0_CUT || dPhiVec->at(1) < AnaConsts::dPhi1_CUT || dPhiVec->at(2) < AnaConsts::dPhi2_CUT ){ passBaseline = false; passBaselineNoTag = false; passdPhis = false; }
  if( debug ) std::cout<<"dPhi0 : "<<dPhiVec->at(0)<<"  dPhi1 : "<<dPhiVec->at(1)<<"  dPhi2 : "<<dPhiVec->at(2)<<"  passBaseline : "<<passBaseline<<std::endl;

  // Pass number of b-tagged jets?
  bool passBJets = true;
  if( !( (AnaConsts::low_nJetsSelBtagged == -1 || cntCSVS >= AnaConsts::low_nJetsSelBtagged) && (AnaConsts::high_nJetsSelBtagged == -1 || cntCSVS < AnaConsts::high_nJetsSelBtagged ) ) ){ passBaseline = false; passBJets = false; }
  if( debug ) std::cout<<"cntCSVS : "<<cntCSVS<<"  passBaseline : "<<passBaseline<<std::endl;

  // Pass the baseline MET requirement?
  bool passMET = true;
  if( metLVec.Pt() < AnaConsts::defaultMETcut ){ passBaseline = false; passBaselineNoTag = false; passMET = false; }
  if( debug ) std::cout<<"met : "<<tr.getVar<double>("met")<<"  defaultMETcut : "<<AnaConsts::defaultMETcut<<"  passBaseline : "<<passBaseline<<std::endl;

  // Calculate top tagger related variables. 
  // Note that to save speed, only do the calculation after previous base line requirements.
  int nTopCandSortedCnt = -1;
  double bestTopJetMass = -1;
  if( passnJets && cntNJetsPt30 >= AnaConsts::nJetsSel ){
      type3Ptr->processEvent((*jetsLVec_forTagger), (*recoJetsBtag_forTagger), metLVec);
      nTopCandSortedCnt = type3Ptr->nTopCandSortedCnt;
      bestTopJetMass = type3Ptr->bestTopJetLVec.M();
  }

  // Pass top tagger requirement?
  bool passTagger = type3Ptr->passNewTaggerReq() && (incZEROtop || nTopCandSortedCnt >= AnaConsts::low_nTopCandSortedSel);

  if( !passTagger ) passBaseline = false;

  //bool passNewCuts = type3Ptr->passNewCuts();

  // Register all the calculated variables
  tr.registerDerivedVar("nMuons_CUT" + spec, nMuons);
  tr.registerDerivedVar("nElectrons_CUT" + spec, nElectrons);
  tr.registerDerivedVar("nIsoTrks_CUT" + spec, nIsoTrks);

  //tr.registerDerivedVar("cntNJetsPt50Eta24" + spec, cntNJetsPt50Eta24);
  tr.registerDerivedVar("cntNJetsPt30Eta24" + spec, cntNJetsPt30Eta24);

  tr.registerDerivedVec("dPhiVec" + spec, dPhiVec);

  tr.registerDerivedVar("cntCSVS" + spec, cntCSVS);

  tr.registerDerivedVec("jetsLVec_forTagger" + spec, jetsLVec_forTagger);
  tr.registerDerivedVec("recoJetsBtag_forTagger" + spec, recoJetsBtag_forTagger);

  tr.registerDerivedVar("cntNJetsPt30" + spec, cntNJetsPt30);

  //tr.registerDerivedVar("passLeptVeto" + spec, passLeptVeto);
  //tr.registerDerivedVar("passMuonVeto" + spec, passMuonVeto);
  //tr.registerDerivedVar("passEleVeto" + spec, passEleVeto);
  //tr.registerDerivedVar("passIsoTrkVeto" + spec, passIsoTrkVeto);
  //tr.registerDerivedVar("passnJets" + spec, passnJets);
  //tr.registerDerivedVar("passdPhis" + spec, passdPhis);
  //tr.registerDerivedVar("passBJets" + spec, passBJets);
  //tr.registerDerivedVar("passMET" + spec, passMET);
  //tr.registerDerivedVar("passTagger" + spec, passTagger);
  tr.registerDerivedVar("passBaseline" + spec, passBaseline);
  //tr.registerDerivedVar("passBaselineNoTag" + spec, passBaselineNoTag);
  //tr.registerDerivedVar("passNewCuts" + spec, passNewCuts);

  tr.registerDerivedVar("nTopCandSortedCnt" + spec, nTopCandSortedCnt);

  tr.registerDerivedVar("best_lept_brJet_MT" + spec, type3Ptr->best_lept_brJet_MT);
  tr.registerDerivedVar("best_had_brJet_MT" + spec, type3Ptr->best_had_brJet_MT);
  tr.registerDerivedVar("best_had_brJet_mTcomb" + spec, type3Ptr->best_had_brJet_mTcomb);
  tr.registerDerivedVar("best_had_brJet_MT2" + spec, type3Ptr->best_had_brJet_MT2);

  double HT = AnaFunctions::calcHT(tr.getVec<TLorentzVector>(jetVecLabel), AnaConsts::pt50Eta24Arr);
  tr.registerDerivedVar("HT" + spec, HT);
  tr.registerDerivedVar("bestTopJetMass" + spec, bestTopJetMass);

  if( debug ) std::cout<<"passBaseline : "<<passBaseline<<"  passBaseline : "<<passBaseline<<std::endl;
} 
