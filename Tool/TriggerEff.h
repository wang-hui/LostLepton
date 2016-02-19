#include<iostream>

//real met efficiency
namespace RealMET
{
  double GetTriggerEffWeight( double met ) 
  {
    if (met<100) return 0.0130;
    else if (met<150) return 0.3239;
    else if (met<175) return 0.7855;
    else if (met<200) return 0.9429;
    else if (met<275) return 0.9736;
    else if (met<400) return 0.9952; 
    else return 1.0;
  }

  double GetTriggerEffStatUncHi( double met ) 
  {
    if (met<100) return 0.0011;
    else if (met<150) return 0.0131;
    else if (met<175) return 0.0226;
    else if (met<200) return 0.0149;
    else if (met<275) return 0.0078;
    else if (met<400) return 0.0040; 
    else return 0.0;
  }

  double GetTriggerEffStatfUncLo( double met ) 
  {
    if (met<100) return 0.0010;
    else if (met<150) return 0.0128;
    else if (met<175) return 0.0244;
    else if (met<200) return 0.0190;
    else if (met<275) return 0.0104;
    else if (met<400) return 0.0109; 
    else return 0.0109;
  }

  double GetTriggerEffSystUncHi( double met ) 
  {
    if (met<100) return 0.0272;
    else if (met<150) return 0.0872;
    else if (met<175) return 0.1505;
    else if (met<200) return 0.0423;
    else if (met<275) return 0.0112;
    else if (met<400) return 0.0001; 
    else return 0.0000;
  }

  double GetTriggerEffSystUncLo( double met ) 
  {
    if (met<100) return 0.0120;
    else if (met<150) return 0.0872;
    else if (met<175) return 0.1505;
    else if (met<200) return 0.0792;
    else if (met<275) return 0.0112;
    else if (met<400) return 0.0001; 
    else return 0.0018;
  }
}

//fake met efficiency
namespace FakeMET
{
  double GetTriggerEffWeight( double met ) 
  {
    if (met<100) return 0.0149;
    else if (met<150) return 0.2482;
    else if (met<175) return 0.5727;
    else if (met<200) return 0.7152;
    else if (met<275) return 0.8098;
    else if (met<400) return 0.8998;
    else return 0.9304; 
  }

  double GetTriggerEffStatUncHi( double met ) 
  {
    if (met<100) return 0.0001;
    else if (met<150) return 0.0017;
    else if (met<175) return 0.0058;
    else if (met<200) return 0.0074;
    else if (met<275) return 0.0061;
    else if (met<400) return 0.0079; 
    else return 0.0138;
  }

  double GetTriggerEffStatfUncLo( double met ) 
  {
    if (met<100) return 0.0001;
    else if (met<150) return 0.0017;
    else if (met<175) return 0.0058;
    else if (met<200) return 0.0075;
    else if (met<275) return 0.0062;
    else if (met<400) return 0.0085; 
    else return 0.0165;
  }
}

double QCDGetTriggerEff( std::string tag, double met )
{
  double METEff = 1;
  if( tag.find("HTMHT") != std::string::npos 
    ||tag.find("SingleMuon") != std::string::npos
    ){ METEff = 1; return METEff; }
  else if( tag.find("QCD") != std::string::npos ){ METEff = FakeMET::GetTriggerEffWeight( met ); return METEff; }
  else if( 
           tag.find("TTJets") != std::string::npos 
        || tag.find("ST_tW_") != std::string::npos
        || tag.find("WJetsToLNu_HT") != std::string::npos
        || tag.find("ZJetsToNuNu_HT") != std::string::npos
        || tag.find("TTZ") != std::string::npos
         ){ METEff = RealMET::GetTriggerEffWeight( met ); return METEff; }
  else{ std::cout << "QCDTag not in the list! What the fuck is going on ??!! Please check TriggerEff.h" << std::endl; return 1; }
}
