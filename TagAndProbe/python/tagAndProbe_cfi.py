import FWCore.ParameterSet.Config as cms

tagAndProbe = cms.EDProducer('TagAndProbe',
  MuonSource   = cms.InputTag('probeMuons'),
  muonsMiniIsoSource = cms.InputTag("probeMuons","muonsMiniIso"),
  jetSrc = cms.InputTag('slimmedJets'),
)
