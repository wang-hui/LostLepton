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

Ratio[1]=(0.9544345);
Ratio[2]=(0.8830448);
Ratio[3]=(0.8905826);
Ratio[4]=(0.9235389);
Ratio[5]=(1.983199);
Ratio[6]=(0.9454236);
Ratio[7]=(0.9721154);
Ratio[8]=(0.9386503);
Ratio[9]=(0.9362441);
Ratio[10]=(0.8711202);
Ratio[11]=(0.9150323);
Ratio[12]=(0.884246);
Ratio[13]=(1.119265);
Ratio[14]=(0.6526005);
Ratio[15]=(0.9757022);
Ratio[16]=(0.9076948);
Ratio[17]=(0.9813681);
Ratio[18]=(1.31217);
Ratio[19]=(1.02893);
Ratio[20]=(0.6580543);
Ratio[21]=(0.6567935);
Ratio[22]=(0.9884901);
Ratio[23]=(0.8677674);
Ratio[24]=(0.9752451);
Ratio[25]=(1.170926);
Ratio[26]=(1.016413);
Ratio[27]=(0.8387792);
Ratio[28]=(0.882192);
Ratio[29]=(0.8134968);
Ratio[30]=(1.058452);
Ratio[31]=(0.8515541);
Ratio[32]=(0.808888);
Ratio[33]=(0.7592695);
Ratio[34]=(0.8309969);
Ratio[35]=(0.5509781);
Ratio[36]=(0.7452685);
Ratio[37]=(0.9290462);
Ratio[38]=(0.8957436);
Ratio[39]=(0.9147847);
Ratio[40]=(1.124083);
Ratio[41]=(0.8750033);
Ratio[42]=(0.991218);
Ratio[43]=(0.861733);
Ratio[44]=(0.8750403);
Ratio[45]=(1.571506);
Ratio[46]=(1.148648);
Ratio[47]=(0.8755001);
Ratio[48]=(0.8293156);
Ratio[49]=(0.9426852);
Ratio[50]=(0.89099);
Ratio[51]=(0.8691972);
Ratio[52]=(0.9540701);
Ratio[53]=(0.8327674);
Ratio[54]=(1.031063);
Ratio[55]=(0.8718982);
Ratio[56]=(0.8633632);
Ratio[57]=(0.9646341);
Ratio[58]=(1.47918);
Ratio[59]=(0.8619262);
Ratio[60]=(0.8362744);
Ratio[61]=(0.9675111);
Ratio[62]=(1.545111);
Ratio[63]=(0.8975248);
Ratio[64]=(1.031416);
Ratio[65]=(0.8289773);
Ratio[66]=(1.053137);
Ratio[67]=(0.9432301);
Ratio[68]=(0.4004564);
Ratio[69]=(0.4607779);
Ratio[70]=(0.9005635);
Ratio[71]=(0.9759313);
Ratio[72]=(1.48784);
Ratio[73]=(0.9656318);
Ratio[74]=(0.7516079);
Ratio[75]=(0.7671837);
Ratio[76]=(1.120266);
Ratio[77]=(0.79928);
Ratio[78]=(0.5898687);
Ratio[79]=(0.8300681);
Ratio[80]=(0.7213006);
Ratio[81]=(1.522591);
Ratio[82]=(0.7952509);
Ratio[83]=(0.5709852);
Ratio[84]=(0.4348763);


ratio_error[1]=(0.01596791);
ratio_error[2]=(0.06288371);
ratio_error[3]=(0.1197618);
ratio_error[4]=(0.2265736);
ratio_error[5]=(0.4581946);
ratio_error[6]=(0.02087802);
ratio_error[7]=(0.0714834);
ratio_error[8]=(0.1405029);
ratio_error[9]=(0.2038543);
ratio_error[10]=(0.05164125);
ratio_error[11]=(0.06283503);
ratio_error[12]=(0.1208679);
ratio_error[13]=(0.226294);
ratio_error[14]=(0.2516294);
ratio_error[15]=(0.2445345);
ratio_error[16]=(0.1540404);
ratio_error[17]=(0.1655876);
ratio_error[18]=(0.4390905);
ratio_error[19]=(0.5501014);
ratio_error[20]=(0.3607544);
ratio_error[21]=(0.2188732);
ratio_error[22]=(0.02315425);
ratio_error[23]=(0.07305376);
ratio_error[24]=(0.1715409);
ratio_error[25]=(0.2718126);
ratio_error[26]=(0.3916002);
ratio_error[27]=(0.05902796);
ratio_error[28]=(0.1444041);
ratio_error[29]=(0.3251018);
ratio_error[30]=(0.4750551);
ratio_error[31]=(0.1894771);
ratio_error[32]=(0.1581925);
ratio_error[33]=(0.2163221);
ratio_error[34]=(0.3422492);
ratio_error[35]=(0.2403051);
ratio_error[36]=(0.4405263);
ratio_error[37]=(0.5561582);
ratio_error[38]=(0.05523135);
ratio_error[39]=(0.1234675);
ratio_error[40]=(0.3237045);
ratio_error[41]=(0.4589862);
ratio_error[42]=(0.1049753);
ratio_error[43]=(0.1634565);
ratio_error[44]=(0.2376159);
ratio_error[45]=(0.5703327);
ratio_error[46]=(0.1353466);
ratio_error[47]=(0.1541921);
ratio_error[48]=(0.1809847);
ratio_error[49]=(0.0593117);
ratio_error[50]=(0.1315849);
ratio_error[51]=(0.2326068);
ratio_error[52]=(0.3396534);
ratio_error[53]=(0.06306222);
ratio_error[54]=(0.1481705);
ratio_error[55]=(0.2140036);
ratio_error[56]=(0.2084165);
ratio_error[57]=(0.278111);
ratio_error[58]=(0.6837773);
ratio_error[59]=(0.05810432);
ratio_error[60]=(0.125184);
ratio_error[61]=(0.2945776);
ratio_error[62]=(0.6034345);
ratio_error[63]=(0.08467055);
ratio_error[64]=(0.1935356);
ratio_error[65]=(0.2642161);
ratio_error[66]=(0.1862919);
ratio_error[67]=(0.4128062);
ratio_error[68]=(0.1635145);
ratio_error[69]=(0.2781128);
ratio_error[70]=(0.3640705);
ratio_error[71]=(0.1935745);
ratio_error[72]=(0.4913127);
ratio_error[73]=(0.4965168);
ratio_error[74]=(0.177019);
ratio_error[75]=(0.2634212);
ratio_error[76]=(0.1684022);
ratio_error[77]=(0.2268759);
ratio_error[78]=(0.2069362);
ratio_error[79]=(0.3193201);
ratio_error[80]=(0.1964104);
ratio_error[81]=(0.4674566);
ratio_error[82]=(0.4569813);
ratio_error[83]=(0.1933512);
ratio_error[84]=(0.1990018);


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
