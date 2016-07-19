#include "ConstantsSnippet.h"
//############determine the pt bin number############
int Set_ptbin_number(
                     double lep_pt
                    )
{
  int ptbin_num;

  //if(gen_pt < 10)
  //{
    //ptbin_num = 0;
  //}
  if(lep_pt >= 10 && lep_pt < 20)
  {
    ptbin_num = 0;
  }
  else if(lep_pt >= 20 && lep_pt < 30)
  {
    ptbin_num = 1;
  }
  else if(lep_pt >= 30 && lep_pt < 40)
  {
    ptbin_num = 2;
  }
  else if(lep_pt >= 40 && lep_pt < 50)
  {
    ptbin_num = 3;
  }
  else if(lep_pt >= 50 && lep_pt < 70)
  {
    ptbin_num = 4;
  }
  else if(lep_pt >= 70 && lep_pt < 100)
  {
    ptbin_num = 5;
  }
  else if(lep_pt >= 100 )
  {
    ptbin_num = 6;
  }
  else
  {
    ptbin_num = -100;
    std::cout << "Invalid Pt bin! Please check!" << std::endl;
  }

  return ptbin_num;
}

//############determine the activity bin number############

int Set_acbin_number(double activity)
{
  int acbin_num;

  if(activity < 0.02)
  {
    acbin_num = 0;
  }
  else if(activity >= 0.02 && activity < 0.05)
  {
    acbin_num = 1;
  }
  else if(activity >= 0.05 && activity < 0.15)
  {
    acbin_num = 2;
  }
  else if(activity >= 0.15 && activity < 1)
  {
    acbin_num = 3;
  }
  else if(activity >= 1)
  {
    acbin_num = 4;
  }
  else
  {
    acbin_num = -100;
    std::cout << "Invalid Activity bin! Please check!" << std::endl;
  }

  return acbin_num;
}

//############determine the nJets bin number############
int Set_HTbin_number(int ht)
{
  int htbin_num;
  if (ht<650) htbin_num=0;
  else if (ht<900) htbin_num=1;
  else htbin_num=2;
  return htbin_num;
}

int Set_MT2bin_number(int mt2)
{
  int mt2bin_num;
  if (mt2<230) mt2bin_num=0;
  else if (mt2<300) mt2bin_num=1;
  else mt2bin_num=2;
  return mt2bin_num;
}

int Set_njetsbin_number(int njets)
{
  int njetsbin_num;

  if(njets == 4)
  {
    njetsbin_num = 0;
  }
  else if(njets == 5)
  {
    njetsbin_num = 1;
  }
  else if(njets == 6)
  {
    njetsbin_num = 2;
  }
  else if(njets == 7)
  {
    njetsbin_num = 3;
  }
  else if(njets ==8)
  {
    njetsbin_num = 4;
  }
  else if(njets >= 9)
  {
    njetsbin_num = 5;
  }
  else
  {
    njetsbin_num = -100;
    std::cout << "Invalid NJet bin! Please check!" << std::endl;
  }

  return njetsbin_num;
}

