#include "TFile.h"
#include "TH2D.h"
#include "ConstantsSnippet.h"

class SystUncs
{
 public:
  double MuRecoEff[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Unc_up[PT_BINS][AC_BINS] = {{0}}, MuRecoEff_Unc_dn[PT_BINS][AC_BINS] = {{0}};
  double MuRecoSystPercent_up = 0.01, MuRecoSystPercent_dn = 0.01;
  double pred_MuReco_up[NSEARCH_BINS] = {0}, pred_MuReco_dn[NSEARCH_BINS] = {0};

  void GetEffAndStatUnc();
	void CombineTotalUnc();
 private:
  double MuRecoEff_Stat_Unc[PT_BINS][AC_BINS] = {{0}};
};
