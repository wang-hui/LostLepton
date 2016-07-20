#include "SystUncs.h"

void SystUncs::GetEffAndStatUnc()
{
  TFile *fin = TFile::Open("RootForPlotting/Effs2dPlots.root");
  TH2D * murecoeff;
  murecoeff = (TH2D*)fin->Get("mus_recoeffs")->Clone();

  for(int i=1;i<PT_BINS+1;i++)
  {
    for(int j=1;j<AC_BINS+1;j++)
    {
      MuRecoEff[i][j] = murecoeff->GetBinContent(i,j+1);
      MuRecoEff_Stat_Unc[i][j] = murecoeff->GetBinError(i,j+1);
    }
  }
}

void SystUncs::CombineTotalUnc()
{
  for(int i=1;i<PT_BINS+1;i++)
  {
    for(int j=1;j<AC_BINS+1;j++)
    {
      MuRecoEff_Unc_up[i][j] = std::sqrt(MuRecoEff_Stat_Unc[i][j]*MuRecoEff_Stat_Unc[i][j] + MuRecoEff[i][j]*MuRecoEff[i][j]*MuRecoSystPercent_up*MuRecoSystPercent_up);
      MuRecoEff_Unc_dn[i][j] = std::sqrt(MuRecoEff_Stat_Unc[i][j]*MuRecoEff_Stat_Unc[i][j] + MuRecoEff[i][j]*MuRecoEff[i][j]*MuRecoSystPercent_dn*MuRecoSystPercent_dn);
    }
  }
}
