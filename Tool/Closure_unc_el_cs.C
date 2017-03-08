{
  const unsigned int nsb=84;
  double Pred[nsb];
  double Ratio[nsb+1];
  double ratio_error[nsb+1];

 for (unsigned int sbc=0;sbc<nsb;++sbc)
 {
   Pred[sbc]=0.0;
   ratio_error[sbc]=0.0;
   Ratio[sbc]=0.0;
 }


//more systclosure2.txt | awk -F\( '{print "Ratio["$2}' | awk -F]= '{print $1"]=("$2}'
//more ratio_error.txt | awk -F\( '{print "ratio_error["$2}' | awk -F, '{print $1"]=("$2}'
//more ratio.txt | awk -F\( '{print "Ratio["$2}' | awk -F, '{print $1"]=("$2}'

Ratio[1]=(0.9884203);
Ratio[2]=(1.033649);
Ratio[3]=(0.7435854);
Ratio[4]=(0.8176571);
Ratio[5]=(1.214167);
Ratio[6]=(0.9591138);
Ratio[7]=(0.97335);
Ratio[8]=(0.7797272);
Ratio[9]=(1.070087);
Ratio[10]=(0.7978055);
Ratio[11]=(0.8994714);
Ratio[12]=(0.8406575);
Ratio[13]=(0.8355354);
Ratio[14]=(0.7124541);
Ratio[15]=(0.6975654);
Ratio[16]=(0.8902523);
Ratio[17]=(0.8275742);
Ratio[18]=(0.5935574);
Ratio[19]=(1.313069);
Ratio[20]=(0.811009);
Ratio[21]=(0.7161122);
Ratio[22]=(1.038124);
Ratio[23]=(1.069919);
Ratio[24]=(0.9909859);
Ratio[25]=(1.171522);
Ratio[26]=(1.759594);
Ratio[27]=(0.8670179);
Ratio[28]=(0.6970008);
Ratio[29]=(0.7760005);
Ratio[30]=(0.8591061);
Ratio[31]=(0.6710428);
Ratio[32]=(0.8741049);
Ratio[33]=(0.8778337);
Ratio[34]=(0.8755505);
Ratio[35]=(0.7582772);
Ratio[36]=(1.66266);
Ratio[37]=(0.7106275);
Ratio[38]=(0.8812125);
Ratio[39]=(1.184184);
Ratio[40]=(1.125379);
Ratio[41]=(1.767031);
Ratio[42]=(1.093713);
Ratio[43]=(1.178614);
Ratio[44]=(0.783576);
Ratio[45]=(1.041209);
Ratio[46]=(1.079829);
Ratio[47]=(0.7392985);
Ratio[48]=(0.7141941);
Ratio[49]=(0.9680466);
Ratio[50]=(0.56811);
Ratio[51]=(0.6808082);
Ratio[52]=(0.6606721);
Ratio[53]=(0.8393321);
Ratio[54]=(0.6693025);
Ratio[55]=(0.7815802);
Ratio[56]=(0.7612);
Ratio[57]=(0.7465641);
Ratio[58]=(0.3799404);
Ratio[59]=(0.9232347);
Ratio[60]=(0.7350386);
Ratio[61]=(0.7736119);
Ratio[62]=(0.8611409);
Ratio[63]=(0.8986869);
Ratio[64]=(1.047536);
Ratio[65]=(0.6231507);
Ratio[66]=(0.7078662);
Ratio[67]=(0.8022052);
Ratio[68]=(0.8473031);
Ratio[69]=(0.6009861);
Ratio[70]=(1.545587);
Ratio[71]=(0.8935586);
Ratio[72]=(1.402141);
Ratio[73]=(0.1548402);
Ratio[74]=(1.034932);
Ratio[75]=(0.3905511);
Ratio[76]=(1.075536);
Ratio[77]=(1.090827);
Ratio[78]=(0.9915192);
Ratio[79]=(1.312174);
Ratio[80]=(1.185088);
Ratio[81]=(0.5830235);
Ratio[82]=(1.175213);
Ratio[83]=(0.6524969);
Ratio[84]=(1.891661);

ratio_error[1]=(0.02073298);
ratio_error[2]=(0.09616983);
ratio_error[3]=(0.1373299);
ratio_error[4]=(0.2305737);
ratio_error[5]=(0.359505);
ratio_error[6]=(0.02536258);
ratio_error[7]=(0.07751819);
ratio_error[8]=(0.1523804);
ratio_error[9]=(0.2551308);
ratio_error[10]=(0.05527153);
ratio_error[11]=(0.0721951);
ratio_error[12]=(0.1416578);
ratio_error[13]=(0.1973698);
ratio_error[14]=(0.2280569);
ratio_error[15]=(0.1882027);
ratio_error[16]=(0.1550803);
ratio_error[17]=(0.1493763);
ratio_error[18]=(0.2387108);
ratio_error[19]=(0.6306409);
ratio_error[20]=(0.4443311);
ratio_error[21]=(0.241366);
ratio_error[22]=(0.03140375);
ratio_error[23]=(0.1021737);
ratio_error[24]=(0.1889654);
ratio_error[25]=(0.362047);
ratio_error[26]=(0.7271961);
ratio_error[27]=(0.07406401);
ratio_error[28]=(0.1397887);
ratio_error[29]=(0.2524594);
ratio_error[30]=(0.4151634);
ratio_error[31]=(0.1849884);
ratio_error[32]=(0.2222374);
ratio_error[33]=(0.2807533);
ratio_error[34]=(0.3345359);
ratio_error[35]=(0.305356);
ratio_error[36]=(0.8644201);
ratio_error[37]=(0.4758851);
ratio_error[38]=(0.07288472);
ratio_error[39]=(0.184783);
ratio_error[40]=(0.3798058);
ratio_error[41]=(1.035037);
ratio_error[42]=(0.1318903);
ratio_error[43]=(0.2242489);
ratio_error[44]=(0.2215891);
ratio_error[45]=(0.3981997);
ratio_error[46]=(0.1002662);
ratio_error[47]=(0.1844165);
ratio_error[48]=(0.2213066);
ratio_error[49]=(0.08006642);
ratio_error[50]=(0.1310076);
ratio_error[51]=(0.2119708);
ratio_error[52]=(0.2952076);
ratio_error[53]=(0.1018506);
ratio_error[54]=(0.1417345);
ratio_error[55]=(0.2351651);
ratio_error[56]=(0.2181376);
ratio_error[57]=(0.2717231);
ratio_error[58]=(0.2300206);
ratio_error[59]=(0.07730054);
ratio_error[60]=(0.1695555);
ratio_error[61]=(0.2538326);
ratio_error[62]=(0.3863418);
ratio_error[63]=(0.1267116);
ratio_error[64]=(0.231197);
ratio_error[65]=(0.2744034);
ratio_error[66]=(0.215731);
ratio_error[67]=(0.4293512);
ratio_error[68]=(0.2378044);
ratio_error[69]=(0.3080019);
ratio_error[70]=(0.7700442);
ratio_error[71]=(0.2033218);
ratio_error[72]=(0.6544813);
ratio_error[73]=(0.07959774);
ratio_error[74]=(0.2458851);
ratio_error[75]=(0.1792213);
ratio_error[76]=(0.2194803);
ratio_error[77]=(0.2555155);
ratio_error[78]=(0.2298276);
ratio_error[79]=(0.5752239);
ratio_error[80]=(0.3053936);
ratio_error[81]=(0.3125893);
ratio_error[82]=(0.3645174);
ratio_error[83]=(0.5924307);
ratio_error[84]=(0.7528444);

 //for (unsigned int sbc=0;sbc<37;++sbc)
 //{
 //  Ratio[sbc]*=pred[sbc];
 //}

 //std::cout.precision(2);
 std::cout.precision(5);
 std::cout << std::fixed;

 for (unsigned int sbc=0;sbc<nsb;++sbc)
 {
   //if (Ratio[sbc+1]>ratio_error[sbc+1]) std::cout << " " << Ratio[sbc+1]/Pred[sbc];
   //else std::cout << " " << ratio_error[sbc+1]/Pred[sbc];
   //if (Ratio[sbc+1]>ratio_error[sbc+1]) std::cout << " " << Ratio[sbc+1];
   //else std::cout << " " << ratio_error[sbc+1];
   //if (std::abs(1.0-ratio_error[sbc+1])>Ratio[sbc+1]) std::cout << " " << std::abs(1.0-ratio_error[sbc+1]);
   //else std::cout << " " << Ratio[sbc+1];
   if (std::abs(1.0-Ratio[sbc+1])>ratio_error[sbc+1]) std::cout << " " << std::abs(1.0-Ratio[sbc+1]);
   else std::cout << " " << ratio_error[sbc+1];
 }
 std::cout << std::endl;

 //for (unsigned int sbc=1;sbc<46;++sbc)
 //{
 //  if (ratio_error[sbc]>Ratio[sbc])  std::cout << " " << -ratio_error[sbc];
 //  else std::cout << " " << -Ratio[sbc];
 //}
 //std::cout << std::endl;

 //for (unsigned int sbc=1;sbc<46;++sbc)
 //{
 //  if (ratio_error[sbc]>Ratio[sbc])  std::cout << " " << ratio_error[sbc]/pred[sbc];
 //  else std::cout << " " << Ratio[sbc]/pred[sbc];
 //}
 //std::cout << std::endl;



}
