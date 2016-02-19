  const double ttbar_mtwcorrfactor[7] = {1.08948,1.1203,1.15576,1.14949,1.16511,1.10954,1.08833};
  const double ttbar_mus_acc[37] = {0.893988,0.887344,0.884429,0.932975,0.898295,0.92276,0.902324,0.852677,0.878031,0.934175,0.880519,0.886159,0.888121,0.891822,0.892749,0.901712,0.903349,0.911844,0.848102,0.887506,0.899916,0.956571,0.954662,0.951531,0.901521,0.943196,0.943899,0.928401,0.924249,0.950237,0.949835,0.951417,0.904682,0.947323,0.930223,0.900489,0.896446};
  const double ttbar_mus_recoeff[7][5] = {{0.967074,0.970334,0.963136,0.961717,0.942104},{0.972964,0.982557,0.995233,0.98292,0.956931},{0.985442,0.980715,0.986284,0.977126,0.991116},{0.994154,0.997404,0.995239,0.98344,0.975489},{0.990752,0.991442,0.985563,0.983474,0.979278},{0.989979,0.992727,0.987739,0.97664,0.967262},{0.992183,0.990902,0.984855,0.971469,0.957987}};
  const double ttbar_mus_isoeff[7][5] = {{0.97575,0.947795,0.939957,0.823867,0.4563},{0.986559,0.968583,0.954654,0.826281,0.569047},{0.981604,0.98089,0.961994,0.79578,0.591509},{0.990263,0.975144,0.965147,0.820812,0.690826},{0.991402,0.982374,0.962076,0.83693,0.734935},{0.988633,0.992423,0.978562,0.928087,0.885404},{0.997301,0.997071,0.994214,0.987146,0.953837}};
  const double ttbar_els_acc[37] = {0.893718,0.872375,0.873775,0.893484,0.919653,0.900366,0.933132,0.83841,0.893395,0.911325,0.920967,0.878107,0.867511,0.877812,0.818209,0.892195,0.89264,0.809101,0.899241,0.894666,0.832858,0.937981,0.949277,0.94484,0.914741,0.943376,0.9289,0.928381,0.912331,0.93931,0.95288,0.950333,0.883437,0.911351,0.937081,0.883915,0.896208};
  const double ttbar_els_recoeff[7][5] = {{0.810734,0.84548,0.766867,0.70774,0.610694},{0.88292,0.85877,0.902208,0.809064,0.715883},{0.923351,0.918439,0.885527,0.836874,0.677156},{0.94249,0.921015,0.938607,0.844962,0.830169},{0.926711,0.927073,0.908588,0.857045,0.740409},{0.941898,0.946637,0.94076,0.872037,0.756808},{0.951667,0.943596,0.915037,0.899843,0.765731}};
  const double ttbar_els_isoeff[7][5] = {{0.933274,0.939303,0.871556,0.788969,0.485563},{0.964741,0.940634,0.929215,0.784611,0.530577},{0.975382,0.967552,0.935811,0.781673,0.589081},{0.9863,0.964127,0.916923,0.793694,0.60057},{0.983339,0.980982,0.953163,0.875331,0.748616},{0.987022,0.985178,0.974876,0.936,0.891747},{0.995839,0.993557,0.989797,0.962004,0.93847}};
  const double ttbar_corrfactor_di_mus = 0.993711;
  const double ttbar_corrfactor_di_els = 0.964613;
  double isoTrackEff_SB[37] = {
0.58265, 0.601425, 0.646018, 0.795819, 0.609429, 0.671346, 0.617747, 0.705308, 0.610776, 0.739056, 0.649424, 0.590203, 0.617803, 0.681299, 0.670849, 0.600549, 0.602114, 0.572457, 0.525956, 0.630472, 0.779585, 0.542299, 0.57843, 0.646578, 0.616819, 0.54589, 0.598434, 0.651662, 0.743122, 0.546526, 0.634391, 0.685644, 0.594616, 0.652518, 0.628895, 0.66193, 0.698136 };
