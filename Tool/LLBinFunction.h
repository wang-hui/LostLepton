#define PT_BINS 7
#define AC_BINS 8
#define NJETS_BINS 6

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

  return ptbin_num;
}

//############determine the activity bin number############

int Set_acbin_number(
                     double activity
                    )
{
  int acbin_num;

  if(activity < 5)
  {
    acbin_num = 0;
  }
  else if(activity >= 5 && activity < 10)
  {
    acbin_num = 1;
  }
  else if(activity >= 10 && activity < 20)
  {
    acbin_num = 2;
  }
  else if(activity >= 20 && activity < 40)
  {
    acbin_num = 3;
  }
  else if(activity >= 40 && activity < 60)
  {
    acbin_num = 4;
  }
  else if(activity >= 60 && activity < 80)
  {
    acbin_num = 5;
  }
  else if(activity >= 80 && activity < 100)
  {
    acbin_num = 6;
  }
  else if(activity >= 100 )
  {
    acbin_num = 7;
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

  return njetsbin_num;
}

