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

Ratio[1]=(0.968537);
Ratio[2]=(0.935385);
Ratio[3]=(0.9129613);
Ratio[4]=(0.8986328);
Ratio[5]=(1.530435);
Ratio[6]=(0.9401805);
Ratio[7]=(0.9890564);
Ratio[8]=(0.9196478);
Ratio[9]=(0.9258184);
Ratio[10]=(0.8583236);
Ratio[11]=(0.8979938);
Ratio[12]=(0.8400889);
Ratio[13]=(1.147406);
Ratio[14]=(0.5807772);
Ratio[15]=(0.880593);
Ratio[16]=(0.9359876);
Ratio[17]=(1.003738);
Ratio[18]=(1.03541);
Ratio[19]=(0.9208812);
Ratio[20]=(0.3771421);
Ratio[21]=(0.7325917);
Ratio[22]=(0.9897927);
Ratio[23]=(0.9117913);
Ratio[24]=(0.9428808);
Ratio[25]=(1.021299);
Ratio[26]=(1.181983);
Ratio[27]=(0.8094993);
Ratio[28]=(0.8362463);
Ratio[29]=(0.7402043);
Ratio[30]=(0.9638116);
Ratio[31]=(0.7289972);
Ratio[32]=(0.865868);
Ratio[33]=(0.7912123);
Ratio[34]=(0.8340461);
Ratio[35]=(0.5752173);
Ratio[36]=(1.233712);
Ratio[37]=(0.6600421);
Ratio[38]=(0.8869914);
Ratio[39]=(0.9654633);
Ratio[40]=(1.037705);
Ratio[41]=(1.121401);
Ratio[42]=(0.9184685);
Ratio[43]=(0.9346037);
Ratio[44]=(0.9113354);
Ratio[45]=(1.859621);
Ratio[46]=(1.11988);
Ratio[47]=(0.9021496);
Ratio[48]=(0.7119925);
Ratio[49]=(1.016488);
Ratio[50]=(0.8012248);
Ratio[51]=(0.8895122);
Ratio[52]=(0.9091942);
Ratio[53]=(0.9041906);
Ratio[54]=(1.061936);
Ratio[55]=(0.7067897);
Ratio[56]=(1.128412);
Ratio[57]=(0.8070478);
Ratio[58]=(1.474097);
Ratio[59]=(0.9134268);
Ratio[60]=(0.8795629);
Ratio[61]=(1.202501);
Ratio[62]=(1.190503);
Ratio[63]=(0.927054);
Ratio[64]=(0.9363595);
Ratio[65]=(0.8355527);
Ratio[66]=(0.878577);
Ratio[67]=(1.060915);
Ratio[68]=(0.5574358);
Ratio[69]=(0.4129339);
Ratio[70]=(0.947288);
Ratio[71]=(0.908693);
Ratio[72]=(1.145919);
Ratio[73]=(0.160938);
Ratio[74]=(0.7925978);
Ratio[75]=(0.733383);
Ratio[76]=(1.184583);
Ratio[77]=(0.6678912);
Ratio[78]=(0.8191631);
Ratio[79]=(1.302093);
Ratio[80]=(0.6203585);
Ratio[81]=(1.066478);
Ratio[82]=(0.5892703);
Ratio[83]=(1.180116);
Ratio[84]=(0.475027);


ratio_error[1]=(0.01776367);
ratio_error[2]=(0.07063173);
ratio_error[3]=(0.131238);
ratio_error[4]=(0.2589536);
ratio_error[5]=(0.4288047);
ratio_error[6]=(0.02227062);
ratio_error[7]=(0.07629891);
ratio_error[8]=(0.1421196);
ratio_error[9]=(0.213899);
ratio_error[10]=(0.05302853);
ratio_error[11]=(0.06305043);
ratio_error[12]=(0.1223674);
ratio_error[13]=(0.2390102);
ratio_error[14]=(0.2361902);
ratio_error[15]=(0.2312047);
ratio_error[16]=(0.1600144);
ratio_error[17]=(0.1754791);
ratio_error[18]=(0.3958552);
ratio_error[19]=(0.523786);
ratio_error[20]=(0.2245677);
ratio_error[21]=(0.252259);
ratio_error[22]=(0.02635747);
ratio_error[23]=(0.08163054);
ratio_error[24]=(0.1777908);
ratio_error[25]=(0.267725);
ratio_error[26]=(0.4590217);
ratio_error[27]=(0.06340259);
ratio_error[28]=(0.1495665);
ratio_error[29]=(0.2959931);
ratio_error[30]=(0.4268352);
ratio_error[31]=(0.1825015);
ratio_error[32]=(0.1719988);
ratio_error[33]=(0.2448922);
ratio_error[34]=(0.3489558);
ratio_error[35]=(0.2547223);
ratio_error[36]=(0.6114975);
ratio_error[37]=(0.4313757);
ratio_error[38]=(0.06030109);
ratio_error[39]=(0.1368159);
ratio_error[40]=(0.3164328);
ratio_error[41]=(0.6133089);
ratio_error[42]=(0.1109097);
ratio_error[43]=(0.1756688);
ratio_error[44]=(0.2612419);
ratio_error[45]=(0.6285859);
ratio_error[46]=(0.0983893);
ratio_error[47]=(0.1681367);
ratio_error[48]=(0.1927795);
ratio_error[49]=(0.07423779);
ratio_error[50]=(0.1400197);
ratio_error[51]=(0.2822171);
ratio_error[52]=(0.3693024);
ratio_error[53]=(0.08751507);
ratio_error[54]=(0.1900716);
ratio_error[55]=(0.2001648);
ratio_error[56]=(0.325644);
ratio_error[57]=(0.2848757);
ratio_error[58]=(0.8055518);
ratio_error[59]=(0.06957191);
ratio_error[60]=(0.1464985);
ratio_error[61]=(0.4038857);
ratio_error[62]=(0.4864875);
ratio_error[63]=(0.1056009);
ratio_error[64]=(0.1976075);
ratio_error[65]=(0.2996778);
ratio_error[66]=(0.210383);
ratio_error[67]=(0.5353317);
ratio_error[68]=(0.1868792);
ratio_error[69]=(0.1792645);
ratio_error[70]=(0.4477796);
ratio_error[71]=(0.2064104);
ratio_error[72]=(0.4614447);
ratio_error[73]=(0.0822959);
ratio_error[74]=(0.1876632);
ratio_error[75]=(0.2704946);
ratio_error[76]=(0.2067722);
ratio_error[77]=(0.1696837);
ratio_error[78]=(0.178012);
ratio_error[79]=(0.5602336);
ratio_error[80]=(0.183838);
ratio_error[81]=(0.4159032);
ratio_error[82]=(0.2435197);
ratio_error[83]=(0.7923935);
ratio_error[84]=(0.2226101);


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
