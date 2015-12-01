#define PT_BINS 7
#define AC_BINS 5
#define NJETS_BINS 6

#define LL_BINS 1
#define NSEARCH_BINS 45
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

int Set_acbin_number(
                     double activity
                    )
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

int Set_njetsbin_number(
                        int njets
                       )
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

