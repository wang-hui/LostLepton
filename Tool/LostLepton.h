#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "Math/QuantFuncMathCore.h"
#include "TMath.h"

#include "NTupleReader.h"

double DeltaPhi(double phi1, double phi2) 
{
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}


double DeltaR(double eta1, double phi1, double eta2, double phi2) 
{
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}
/*
bool isGoodMuonID(double id_GlobalMuonPromptTight,
                  //double isGlobalMuon,
                  double isPFMuon,
                  //double tk_chi2,
                  //double tk_numvalhits,
                  double numberOfMatchedStations,
                  double id_d0,
                  double id_tk_dz,
                  double tk_numvalPixelhits,
                  double tk_LayersWithMeasurement
                  )
{

  if((int)id_GlobalMuonPromptTight==0)
    return false;
  //if((int)isGlobalMuon==0)
    //return false;
  if((int)isPFMuon==0)
    return false;
  
  //if(tk_chi2>=10)
    //return false;
  //if(tk_numvalhits<=0)
    //return false;
  if(numberOfMatchedStations<=1)
    return false;
  if(id_d0>=0.2 || id_d0<=-0.2)
    return false;
  if(id_tk_dz>=0.5 || id_tk_dz<=-0.5)
    return false;
  if(tk_numvalPixelhits<=0)
    return false;
  if(tk_LayersWithMeasurement<=5)
    return false;
  
  return true;
}
*/
/*
bool isGoodElectronID(double isEB,
                      double isEE,
                      double dEtaIn,
                      double dPhiIn,
                      double full5x5_sigmaIetaIeta,
                      double hadOverEm
                      )
{ 
  if(( (int)isEB==1 && (int)isEE==1 )||( (int)isEB==0 && (int)isEE==0 ) )
    return false;

  if((int)isEB==1)
  {
    if(dEtaIn > 0.021)
      return false;
    if(dPhiIn > 0.25)
      return false;
    if(full5x5_sigmaIetaIeta > 0.012)
      return false;
    if(hadOverEm > 0.24)
      return false;
  }
  
  if((int)isEE==1)
  {
    if(dEtaIn > 0.028)
      return false;
    if(dPhiIn > 0.23)
      return false;
    if(full5x5_sigmaIetaIeta > 0.035)
      return false;
    if(hadOverEm > 0.19)
      return false;
  }

  return true;
}
*/
/*
bool isMatchedDeltaR(double eta1, double phi1, double eta2, double phi2)
{
  double deltaR;
  deltaR = DeltaR(eta1,phi1,eta2,phi2);

  if(deltaR<0.01)
    return true;
  else
    return false;
}
*/

/*
bool isPassMuonIso(double mus_pfIsoR04_sumNeutralHadronEt,
                   double mus_pfIsoR04_sumPhotonEt,
                   double mus_pfIsoR04_sumPUPt,
                   double mus_pfIsoR04_sumChargedHadronPt,
                   double reco_mus_pt
                   )
{
  double isoNeutral = mus_pfIsoR04_sumNeutralHadronEt + mus_pfIsoR04_sumPhotonEt - 0.5*mus_pfIsoR04_sumPUPt;
  isoNeutral = ( isoNeutral > 0) ? isoNeutral : 0;

  double  muRelIso = (mus_pfIsoR04_sumChargedHadronPt + isoNeutral) / reco_mus_pt;
  if( muRelIso < 0.2 )
  {
    return true;
  }
  else
    return false;
}
*/



double get_stat_Error(double a,
                      double an
                     )
{
  double n;
  n = an - a;

  double err;
  err = 1000;

  double alpha;
  alpha = 1-0.6827;

  if( a>=0 && n>=0 )
  {
    err = std::sqrt(n/(a+n)/(a+n)*n/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,a+1,1)+a/(a+n)/(a+n)*a/(a+n)/(a+n)*ROOT::Math::gamma_quantile_c(alpha/2,n+1,1));
    return err;
  }
  else
  {
    return -1;
  }
}

double get_sys_Error(double r,
                     double p
                    )
{
  double err;
  err = 1000;

  err = r - (1-(1-r)*(1+p));

  if(err < 0)
  {
    return -1;
  }
  return err;
}

