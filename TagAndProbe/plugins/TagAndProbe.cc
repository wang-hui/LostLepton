// -*- C++ -*-
//
// Package:    LostLepton/TagAndProbe
// Class:      TagAndProbe
// 
/**\class TagAndProbe TagAndProbe.cc LostLepton/TagAndProbe/plugins/TagAndProbe.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Florent Sylvain Lacroix
//         Created:  Fri, 22 May 2015 13:37:22 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "LostLepton/Tool/Activity.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//
// class declaration
//

class TagAndProbe : public edm::EDProducer {
   public:
      explicit TagAndProbe(const edm::ParameterSet&);
      ~TagAndProbe();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      edm::InputTag muonSrc_;
      edm::InputTag muonsMiniIsoSrc_;
      edm::InputTag jetSrc_;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TagAndProbe::TagAndProbe(const edm::ParameterSet& iConfig)
{
  muonSrc_       = iConfig.getParameter<edm::InputTag>("MuonSource");
  muonsMiniIsoSrc_ = iConfig.getParameter<edm::InputTag>("muonsMiniIsoSource");
  jetSrc_ = iConfig.getParameter<edm::InputTag>("jetSrc");

  produces<edm::ValueMap<float> >("miniIso");
  produces<edm::ValueMap<float> >("activity");
  
}


TagAndProbe::~TagAndProbe()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
TagAndProbe::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle<std::vector<pat::Muon> > muons;
  iEvent.getByLabel(muonSrc_, muons);

  edm::Handle<std::vector<double> > muonsMiniIso;
  iEvent.getByLabel(muonsMiniIsoSrc_, muonsMiniIso);

  edm::Handle<std::vector<pat::Jet> > jets;
  iEvent.getByLabel(jetSrc_, jets);

  std::vector<TLorentzVector> jetsLVec;
  std::vector<double> recoJetschargedHadronEnergyFraction;
  std::vector<double> recoJetschargedEmEnergyFraction;
  std::vector<pat::Jet> extJets = (*jets);

  for(unsigned int ij=0; ij < extJets.size(); ij++)
  {
    const pat::Jet& jet = extJets[ij];
    TLorentzVector perJetLVec;
    perJetLVec.SetPtEtaPhiE( jet.pt(), jet.eta(), jet.phi(), jet.energy() );
    jetsLVec.push_back(perJetLVec);
    double chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
    recoJetschargedHadronEnergyFraction.push_back( chargedHadronEnergyFraction );
    double chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
    recoJetschargedEmEnergyFraction.push_back( chargedEmEnergyFraction );
  }

  // prepare vector for output
  std::vector<float> values;
  std::vector<float> values2;
  //for (std::vector<pat::Muon>::const_iterator m = muons->begin(); m != muons->end(); ++m)
  for (unsigned int mc = 0; mc<muons->size(); ++mc)
  {
    //std::cout << "muonsMiniIso[" << mc << "] = " << (*muonsMiniIso)[mc] << std::endl;
    values.push_back((*muonsMiniIso)[mc]);

    Activity myActivity;
    //std::cout << "muonc eta = " << ((*muons)[mc]).eta() << std::endl;
    //std::cout << "muonc phi = " << ((*muons)[mc]).phi() << std::endl;
    myActivity.getVariables(((*muons)[mc]).eta(), ((*muons)[mc]).phi(), jetsLVec, recoJetschargedHadronEnergyFraction, recoJetschargedEmEnergyFraction);
    const float activity = myActivity.getMuActivity();
    myActivity.reset();
    values2.push_back(activity);

  }

  // convert into ValueMap and store
  std::auto_ptr<ValueMap<float> > valMap(new ValueMap<float>());
  ValueMap<float>::Filler filler(*valMap);
  filler.insert(muons, values.begin(), values.end());
  filler.fill();
  iEvent.put(valMap, "miniIso");

  std::auto_ptr<ValueMap<float> > valMap2(new ValueMap<float>());
  ValueMap<float>::Filler filler2(*valMap2);
  filler2.insert(muons, values2.begin(), values2.end());
  filler2.fill();
  iEvent.put(valMap2, "activity");

 
}

// ------------ method called once each job just before starting event loop  ------------
void 
TagAndProbe::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TagAndProbe::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TagAndProbe::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TagAndProbe::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TagAndProbe::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TagAndProbe::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TagAndProbe::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagAndProbe);
