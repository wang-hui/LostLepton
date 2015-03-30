#include <algorithm>
#include <vector>
#include "TLorentzVector.h"

class Activity
{
 public:
  void getVariables(
                    double eta_lepton,
                    double phi_lepton,
                    const std::vector<TLorentzVector>& jetLVec,
                    const std::vector<double>& recoJetschargedHadronEnergyFraction,
                    const std::vector<double>& recoJetschargedEmEnergyFraction
                   );
  double getActivity();
  void reset();

 private:
  double eta_lepton_in = 0;
  double phi_lepton_in = 0;
  std::vector<TLorentzVector> jetLVec_in;
  std::vector<double> recoJetschargedHadronEnergyFraction_in;
  std::vector<double> recoJetschargedEmEnergyFraction_in;

  double DeltaPhi(double phi1, double phi2);
  double DeltaR(double eta1, double phi1, double eta2, double phi2);
};

void Activity::getVariables(
                            double eta_lepton,
                            double phi_lepton, 
                            const std::vector<TLorentzVector>& jetLVec,
                            const std::vector<double>& recoJetschargedHadronEnergyFraction,
                            const std::vector<double>& recoJetschargedEmEnergyFraction
                           )
{
  eta_lepton_in = eta_lepton;
  phi_lepton_in = phi_lepton;
  jetLVec_in.resize( jetLVec.size() );
  std::copy ( jetLVec.begin() , jetLVec.end() , jetLVec_in.begin() );
  recoJetschargedHadronEnergyFraction_in.resize( recoJetschargedHadronEnergyFraction.size() );
  std::copy ( recoJetschargedHadronEnergyFraction.begin() , recoJetschargedHadronEnergyFraction.end() , recoJetschargedHadronEnergyFraction_in.begin() );
  recoJetschargedEmEnergyFraction_in.resize( recoJetschargedEmEnergyFraction.size() );
  std::copy ( recoJetschargedEmEnergyFraction.begin() , recoJetschargedEmEnergyFraction.end() , recoJetschargedEmEnergyFraction_in.begin() );
}

double Activity::getActivity()
{
  double activity = 0;
  for( unsigned int i = 0 ; i < jetLVec_in.size() ; i++ )
  {
    if( DeltaR( eta_lepton_in , phi_lepton_in , (jetLVec_in.at(i)).Eta() , (jetLVec_in.at(i)).Phi() ) < 1.0 )
    {
      activity+= (jetLVec_in.at(i)).Pt() * (recoJetschargedEmEnergyFraction_in[i] + recoJetschargedHadronEnergyFraction_in[i]);
    }
    else
      continue;
  }

  return activity;
}

void Activity::reset()
{
  eta_lepton_in = 0;
  phi_lepton_in = 0;
  jetLVec_in.clear();
  recoJetschargedHadronEnergyFraction_in.clear();
  recoJetschargedEmEnergyFraction_in.clear();
}

double Activity::DeltaPhi(double phi1, double phi2)
{
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double Activity::DeltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}



