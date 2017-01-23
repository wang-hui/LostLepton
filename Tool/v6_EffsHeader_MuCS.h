  double ttbar_mtwcorrfactor[7] = {1.04193,1.06267,1.07404,1.08663,1.11306,1.18561,1.65656};
  double ttbar_mus_acc[84] = {0.780193,0.82393,0.835125,0.769594,0.92018298817,0.779563,0.812786,0.860071,0.856293,0.838565,0.797106,0.805551,0.739294,0.8693357056,0.847759,0.862227,0.823519,0.92018298817,0.83023204195,0.83023204195,0.72810308556,0.786707,0.836489,0.886124,0.791415,0.92762955774,0.82535,0.768684,0.924996,0.91924762334,0.839694,0.891171,0.648817,0.7681374129,0.98157375349,0.98157375349,0.92762955774,0.777854,0.751782,0.855028,0.85675474408,0.874987,0.726759,0.875533,0.699366,0.868238,0.916838,0.91870532928,0.912479,0.947013,0.943535,0.91613353663,0.815765,0.8287  ,0.921988,0.8453  ,0.816404,0.91613353663,0.922984,0.925083,0.978797,0.93268490502,0.805151,0.83331,0.830446,0.819422,0.889238,0.93268490502,0.88376,0.82742708279,0.850515,0.813213,0.91549745287,0.899288,0.933333,0.917803,1,0.91549745287,0.915387,0.92552234269,0.974176,0.96126901587,0.98633526132,0.98633526132};
  double ttbar_mus_recoeff[7][5] = {{0.964408,0.96532,0.969147,0.95984,0.950648},{0.974389,0.986208,0.976967,0.979465,0.971187},{0.984965,0.982843,0.980217,0.983335,0.966799},{0.979301,0.985525,0.979976,0.97267,0.925914},{0.986028,0.980716,0.97916,0.981516,0.965567},{0.978663,0.984728,0.975522,0.978571,0.96401},{0.913796,0.931964,0.93648,0.938051,0.985049}};
  double ttbar_mus_isoeff[7][5] = {{0.97614,0.970054,0.940133,0.836283,0.460842},{0.977676,0.973717,0.958495,0.849012,0.591276},{0.989179,0.990064,0.97591,0.871752,0.617503},{0.986618,0.985757,0.978102,0.861549,0.615821},{0.994342,0.99403,0.987992,0.891221,0.735342},{0.997674,0.995598,0.994254,0.958739,0.915876},{0.998144,0.99728,0.998538,0.98572,0.976709}};
  double ttbar_els_acc[84] = {0.783808,0.822412,0.819415,0.869192,0.81225645483,0.801268,0.793529,0.847683,0.807488,0.823154,0.795892,0.795172,0.880392,0.91777690471,0.910419,0.830248,0.819692,0.81225645483,0.88581878044,0.88581878044,0.88581878044,0.809308,0.819631,0.868003,0.699888,0.944690084,0.839415,0.796494,0.797563,0.74976375022,0.747043,0.802366,0.827576,0.81182570731,0.91527731034,0.91527731034,0.944690084,0.788042,0.788136,0.781693,0.795206077,0.867721,0.843156,0.817941,0.88733,0.874958,0.908222,0.91221164719,0.886755,0.921615,0.781377,0.96731901867,0.79078,0.862422,0.868224,0.892126,0.835305,0.96731901867,0.911929,0.906479,0.907165,0.93145303237,0.81403,0.885412,0.824022,0.789837,0.689519,0.93145303237,0.815704,0.79974031804,0.812691,0.822507,0.92848121132,0.932499,0.953683,0.899394,0.920393,0.92848121132,0.916087,0.90494159848,0.787582,0.80155764506,0.75807141613,0.75807141613};
  double ttbar_els_recoeff[7][5] = {{0.77249,0.780312,0.763259,0.684666,0.552774},{0.876397,0.857118,0.865194,0.813286,0.720952},{0.900748,0.884238,0.886188,0.849244,0.771461},{0.879998,0.884294,0.895708,0.846441,0.822061},{0.903027,0.894578,0.879066,0.839691,0.854942},{0.898415,0.883988,0.876204,0.835222,0.787409},{0.88445,0.839395,0.857475,0.844684,0.776513}};
  double ttbar_els_isoeff[7][5] = {{0.901282,0.901512,0.823654,0.723731,0.360348},{0.938337,0.924106,0.875705,0.74488,0.446034},{0.970133,0.956317,0.907111,0.769051,0.513676},{0.965441,0.961667,0.931106,0.755764,0.474133},{0.988162,0.978673,0.962848,0.870834,0.715482},{0.996453,0.98689,0.970765,0.922334,0.911529},{0.992467,0.98965,0.979693,0.9458,0.913362}};
  double ttbar_corrfactor_di_mus = 0.993714;
  double ttbar_corrfactor_di_els = 0.971324;
  double isoTrackEff_SB[84] = {0.594375, 0.546478, 0.7282, 0.688987, 0.595184, 0.591746, 0.630295, 0.563533, 0.614842, 0.633231, 0.599336, 0.608467, 0.574718, 0.49790383114, 0.519775, 0.697071, 0.524727, 0.444574, 0.79088364396, 0.79088364396, 0.79088364396, 0.592583, 0.577918, 0.603168, 0.440157, 0.577365, 0.591467, 0.482023, 0.643849, 0.63944231766, 0.6993, 0.570809, 0.45026, 0.460367, 0.71155236655 , 0.71155236655, 0.71155236655, 0.568682, 0.627962, 0.71398, 0.6915667528, 0.541915, 0.525106, 0.730244, 0.710579, 0.709875, 0.653734, 0.864201, 0.546497, 0.545002, 0.836265, 0.80460967188, 0.636125, 0.567026, 0.441162, 0.529601, 0.693581, 0.70000631589, 0.526511, 0.548097, 0.50195, 0.52638378831, 0.561576, 0.585059, 0.609084, 0.561132, 0.593056, 0.56913185298, 0.906777, 0.475208, 0.490372, 0.660005, 0.5, 0.501237, 0.64193, 0.700626, 0.606326, 0.55556847701, 0.578813, 0.61165800818, 0.621471, 0.59577751302, 0.604189, 0.55815360214};
