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

Ratio[1]=(0.9718345);
Ratio[2]=(0.986182);
Ratio[3]=(0.8009616);
Ratio[4]=(0.7874288);
Ratio[5]=(1.416859);
Ratio[6]=(0.9463098);
Ratio[7]=(0.9719813);
Ratio[8]=(0.8365575);
Ratio[9]=(0.8815675);
Ratio[10]=(0.8106605);
Ratio[11]=(0.9191522);
Ratio[12]=(0.9333003);
Ratio[13]=(0.7287472);
Ratio[14]=(0.6248317);
Ratio[15]=(0.7573403);
Ratio[16]=(0.8646906);
Ratio[17]=(0.8642015);
Ratio[18]=(0.6944388);
Ratio[19]=(1.597669);
Ratio[20]=(1.300004);
Ratio[21]=(0.5948491);
Ratio[22]=(1.0218);
Ratio[23]=(0.9984454);
Ratio[24]=(1.085972);
Ratio[25]=(1.210296);
Ratio[26]=(1.55508);
Ratio[27]=(0.8610924);
Ratio[28]=(0.6905412);
Ratio[29]=(0.7542769);
Ratio[30]=(0.690005);
Ratio[31]=(0.6635276);
Ratio[32]=(1.023728);
Ratio[33]=(0.6773096);
Ratio[34]=(0.8726796);
Ratio[35]=(0.6186126);
Ratio[36]=(0.720963);
Ratio[37]=(1.382609);
Ratio[38]=(0.8853608);
Ratio[39]=(1.051431);
Ratio[40]=(1.266609);
Ratio[41]=(1.467028);
Ratio[42]=(1.00727);
Ratio[43]=(0.9775799);
Ratio[44]=(0.7359429);
Ratio[45]=(1.147687);
Ratio[46]=(1.039578);
Ratio[47]=(0.7124492);
Ratio[48]=(0.7682716);
Ratio[49]=(0.8734521);
Ratio[50]=(0.6663753);
Ratio[51]=(0.6217177);
Ratio[52]=(0.708963);
Ratio[53]=(0.9078072);
Ratio[54]=(0.730188);
Ratio[55]=(0.793012);
Ratio[56]=(0.8114608);
Ratio[57]=(0.6084986);
Ratio[58]=(0.5895243);
Ratio[59]=(0.8858714);
Ratio[60]=(0.7441432);
Ratio[61]=(0.7161821);
Ratio[62]=(1.030901);
Ratio[63]=(0.9619501);
Ratio[64]=(0.9671668);
Ratio[65]=(0.6765407);
Ratio[66]=(0.8866561);
Ratio[67]=(0.4762409);
Ratio[68]=(1.541928);
Ratio[69]=(0.954134);
Ratio[70]=(1.065689);
Ratio[71]=(0.9548199);
Ratio[72]=(1.176354);
Ratio[73]=(1.219151);
Ratio[74]=(1.111632);
Ratio[75]=(0.5416517);
Ratio[76]=(1.055069);
Ratio[77]=(1.152274);
Ratio[78]=(0.9740988);
Ratio[79]=(1.461218);
Ratio[80]=(0.45029);
Ratio[81]=(0.8296573);
Ratio[82]=(0.9007439);
Ratio[83]=(0.4129778);
Ratio[84]=(0.645457);


ratio_error[1]=(0.01807439);
ratio_error[2]=(0.08326459);
ratio_error[3]=(0.1349349);
ratio_error[4]=(0.1918357);
ratio_error[5]=(0.3675581);
ratio_error[6]=(0.02286612);
ratio_error[7]=(0.07432087);
ratio_error[8]=(0.1481403);
ratio_error[9]=(0.2538467);
ratio_error[10]=(0.05244661);
ratio_error[11]=(0.0684373);
ratio_error[12]=(0.1374817);
ratio_error[13]=(0.162112);
ratio_error[14]=(0.1922385);
ratio_error[15]=(0.1984148);
ratio_error[16]=(0.1507403);
ratio_error[17]=(0.1512814);
ratio_error[18]=(0.2523047);
ratio_error[19]=(0.733111);
ratio_error[20]=(0.657729);
ratio_error[21]=(0.2083145);
ratio_error[22]=(0.02656109);
ratio_error[23]=(0.09291625);
ratio_error[24]=(0.1932896);
ratio_error[25]=(0.3343525);
ratio_error[26]=(0.6129504);
ratio_error[27]=(0.0659157);
ratio_error[28]=(0.1175864);
ratio_error[29]=(0.2319808);
ratio_error[30]=(0.2997899);
ratio_error[31]=(0.1636283);
ratio_error[32]=(0.2080988);
ratio_error[33]=(0.1938647);
ratio_error[34]=(0.3253316);
ratio_error[35]=(0.2540754);
ratio_error[36]=(0.4328768);
ratio_error[37]=(0.8662748);
ratio_error[38]=(0.06351929);
ratio_error[39]=(0.1538868);
ratio_error[40]=(0.4136842);
ratio_error[41]=(0.8314365);
ratio_error[42]=(0.1211436);
ratio_error[43]=(0.2118333);
ratio_error[44]=(0.2039387);
ratio_error[45]=(0.4182139);
ratio_error[46]=(0.1247009);
ratio_error[47]=(0.1715393);
ratio_error[48]=(0.263227);
ratio_error[49]=(0.06087764);
ratio_error[50]=(0.1187999);
ratio_error[51]=(0.1641709);
ratio_error[52]=(0.2707877);
ratio_error[53]=(0.08408691);
ratio_error[54]=(0.119343);
ratio_error[55]=(0.1981037);
ratio_error[56]=(0.1936598);
ratio_error[57]=(0.2218239);
ratio_error[58]=(0.3111376);
ratio_error[59]=(0.06470924);
ratio_error[60]=(0.1402417);
ratio_error[61]=(0.1998314);
ratio_error[62]=(0.4397516);
ratio_error[63]=(0.106462);
ratio_error[64]=(0.197678);
ratio_error[65]=(0.2691838);
ratio_error[66]=(0.1530444);
ratio_error[67]=(0.2780859);
ratio_error[68]=(0.4250503);
ratio_error[69]=(0.6238258);
ratio_error[70]=(0.4814053);
ratio_error[71]=(0.1873596);
ratio_error[72]=(0.4294216);
ratio_error[73]=(0.5620064);
ratio_error[74]=(0.2428761);
ratio_error[75]=(0.1990581);
ratio_error[76]=(0.1807634);
ratio_error[77]=(0.3145234);
ratio_error[78]=(0.2002069);
ratio_error[79]=(0.4978858);
ratio_error[80]=(0.261787);
ratio_error[81]=(0.3113966);
ratio_error[82]=(0.5461902);
ratio_error[83]=(0.2311598);
ratio_error[84]=(0.4065562);


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
