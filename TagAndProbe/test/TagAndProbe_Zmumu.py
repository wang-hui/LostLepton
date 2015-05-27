import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

######### EXAMPLE CFG 
###  A simple test of runnning T&P on Zmumu to determine muon isolation and identification efficiencies
###  More a showcase of the tool than an actual physics example

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Source
process.source = cms.Source("PoolSource",
          fileNames = cms.untracked.vstring ("root://cmsxrootd.fnal.gov//store/mc/Phys14DR/DYJetsToLL_M-50_13TeV-madgraph-pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/06C61714-7E6C-E411-9205-002590DB92A8.root")
          #fileNames = cms.untracked.vstring ("file:06C61714-7E6C-E411-9205-002590DB92A8.root")
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )    

## Tags. In a real analysis we should require that the tag muon fires the trigger, 
##       that's easy with PAT muons but not RECO/AOD ones, so we won't do it here
##       (the J/Psi example shows it)
process.load("SusyAnaTools.SkimsAUX.prodMuons_cfi")
process.load("SusyAnaTools.SkimsAUX.prodJets_cfi")

process.goodVertices = cms.EDFilter(
  "VertexSelector",
  filter = cms.bool(False),
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
)

## Tags. In a real analysis we should require that the tag muon fires the trigger, that's easy with PAT muons but not RECO/AOD ones, so we won't do it here. (the J/Psi example shows it)
process.tagMuons = process.prodMuons.clone()
process.tagMuons.MinMuPt = cms.double(20.0)
process.tagMuons.MaxMuEta = cms.double(2.0)
process.tagMuons.DoMuonIsolation = cms.bool(False)
process.tagMuons.DoMuonMiniIsolation = cms.bool(True)

## Probes
process.probeMuons = process.prodMuons.clone()
process.probeMuons.MinMuPt = cms.double(5.0)
process.probeMuons.MaxMuEta = cms.double(2.4)
process.probeMuons.DoMuonIsolation = cms.bool(False)

process.load("LostLepton.TagAndProbe.tagAndProbe_cfi")
process.TagAndProbeVariables = process.tagAndProbe.clone()
process.TagAndProbeVariables.tagMuons = cms.InputTag('tagMuons')


## Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("70.0 < mass < 130.0"),
)



## Match muons to MC
#process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
#    pdgId = cms.vint32(13),
#    src = cms.InputTag("muons"),
#    distMin = cms.double(0.3),
#    matched = cms.InputTag("genParticles")
#)

## Make the tree
process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # pairs
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # variables to use
    variables = cms.PSet(
        ## methods of reco::Candidate
        eta = cms.string("eta"),
        pt  = cms.string("pt"),
        #muonsMiniIso = cms.string("muonsMiniIso"),
        #miniIso = cms.string("miniIso"),
        ## a method of the reco::Muon object (thanks to the 3.4.X StringParser)
        #nsegm = cms.string("numberOfMatches"), 
        ## this one is an external variable
        #drj = cms.InputTag("drToNearestJet"),
        #nGV = cms.InputTag("goodVertices.size()"),
        #miniIso = cms.InputTag("probeMuons","muonsMiniIso"),
        #miniIso = cms.InputTag("probeMuons","miniIso"),
        miniIso = cms.InputTag("TagAndProbeVariables","miniIso"),
        activity = cms.InputTag("TagAndProbeVariables","activity"),
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ## one defined by an external collection of passing probes
        #passingCal = cms.InputTag("probesPassingCal"), 
        ## two defined by simple string cuts
        passingGlb = cms.string("isGlobalMuon"),
        #passingIso = cms.string("muonsMiniIso<0.2"),
        #passingIso = cms.InputTag("probePassingIso"),
    ),
    # mc-truth info
    isMC = cms.bool(False),
    motherPdgId = cms.vint32(22,23),
    #makeMCUnbiasTree = cms.bool(True),
    #checkMotherInUnbiasEff = cms.bool(True),
    #tagMatches = cms.InputTag("muMcMatch"),
    #probeMatches  = cms.InputTag("muMcMatch"),
    #allProbes     = cms.InputTag("probeMuons"),
)
##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
process.tagAndProbe = cms.Path( 
    process.goodVertices *
    (process.tagMuons + process.probeMuons) *   # 'A*B' means 'B needs output of A' 'A+B' means 'if you want you can re-arrange the order'
    (process.tpPairs) * process.TagAndProbeVariables *
    process.muonEffs
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("testTagProbeFitTreeProducer_ZMuMu.root"))

#process.display = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('display_ttbar_fastsim_310pre10.root'),
#    outputCommands = cms.untracked.vstring(
#        'keep *'
#    )
#)
#process.outpath = cms.EndPath(process.display)

