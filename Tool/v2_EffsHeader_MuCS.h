  const double ttbar_mtwcorrfactor[7] = {1.0489,1.08228,1.09569,1.12216,1.15341,1.21229,1.66882};
  const double ttbar_mus_acc[69] = {0.826654,0.790422,0.853053,0.836018,0.855013,0.806637,0.817756,0.813614,0.932188,0.829861,0.865473,0.870932,0.826392,0.78402,0.854565,0.887736,0.848666,0.80455,0.717903,0.662087,0.792813,0.78136,0.468032,0.763785,0.829625,0.870847,0.770509,0.846343,0.780344,0.934779,0.894247,0.903966,0.903519,0.780827,0.82567,0.815524,0.867723,0.994045,0.960914,0.854931,0.868011,0.875747,0.90127,0.914258,0.902125,1,0.842876,0.82205,0.802963,0.81743,0.841172,0.749516,0.632915,0.782489,0.894722,0.913091,1,0.773945,0.889123,0.656878,0.885187,0.946405,0.893779,0.919449,0.887991,0.994867,0.908563,1,1};
  const double ttbar_mus_recoeff[7][5] = {{0.97757,0.969873,0.965631,0.953099,0.96145},{0.981485,0.985595,0.981391,0.97938,0.972094},{0.987326,0.991531,0.989286,0.980015,0.986805},{0.981692,0.986204,0.986938,0.991805,0.961839},{0.986111,0.989642,0.989055,0.984345,0.96924},{0.986647,0.985886,0.984101,0.958744,0.960267},{0.954972,0.974749,0.970804,0.953664,0.991724}};
  const double ttbar_mus_isoeff[7][5] = {{0.976805,0.960362,0.936416,0.835856,0.455684},{0.978089,0.975318,0.950507,0.825607,0.561071},{0.992942,0.984556,0.969419,0.849932,0.608887},{0.987595,0.985429,0.966941,0.826598,0.666476},{0.992009,0.992491,0.984482,0.894254,0.772629},{0.997313,0.997008,0.990712,0.951208,0.901405},{0.99641,0.997014,0.998276,0.984904,0.92814}};
  const double ttbar_els_acc[69] = {0.844153,0.801092,0.844289,0.780409,0.825387,0.805292,0.814966,0.885285,0.899731,0.797879,0.818245,0.935399,0.829605,0.799304,0.851686,0.811315,0.875146,0.833682,0.721366,0.795379,0.998018,0.826519,0.877385,0.786196,0.828854,0.839014,0.816552,0.848405,0.850378,0.958754,0.896056,0.900667,0.94036,0.991229,0.827426,0.859295,0.913409,0.982899,0.843043,0.807333,0.831561,0.418362,0.905483,0.90358,0.826544,1,0.813224,0.877972,0.93889,0.648385,0.89119,0.778059,0.791618,0.863262,0.8833,0.902062,1,0.804188,0.790657,0.762895,0.919288,0.889932,1,0.896836,0.984696,1,0.949056,1,1};
  const double ttbar_els_recoeff[7][5] = {{0.789447,0.822228,0.792728,0.705712,0.547003},{0.905342,0.892278,0.888698,0.834654,0.725883},{0.938773,0.917128,0.919071,0.862649,0.734204},{0.919304,0.93242,0.904567,0.87073,0.824924},{0.93334,0.919459,0.901739,0.848637,0.843719},{0.928724,0.923061,0.912459,0.885021,0.767513},{0.927758,0.916681,0.887866,0.876542,0.699021}};
  const double ttbar_els_isoeff[7][5] = {{0.890221,0.860531,0.854124,0.720503,0.363349},{0.938323,0.916269,0.909245,0.760213,0.476106},{0.967181,0.961907,0.921274,0.78245,0.530351},{0.967368,0.974895,0.942114,0.762851,0.513852},{0.984546,0.982689,0.962916,0.860603,0.674364},{0.991967,0.985128,0.985159,0.927422,0.88907},{0.993544,0.991588,0.982422,0.944961,0.952693}};
  const double ttbar_corrfactor_di_mus = 0.992976;
  const double ttbar_corrfactor_di_els = 0.966526;
  double isoTrackEff_SB[69] = {
0.59178, 0.603306, 0.644624, 0.538872, 0.607068, 0.60189, 0.62154, 0.523643, 0.498594, 0.568278, 0.508693, 0.74624, 0.587798, 0.635241, 0.671636, 0.508456, 0.58892, 0.629514, 0.625723, 0.928806, 0.617246, 0.685992, 0.252431, 0.797209, 0.563732, 0.672117, 0.844977, 0.745292, 0.834687, 0, 0.5556, 0.542734, 0.538438, 0.662089, 0.619676, 0.630295, 0.701177, 0.590044, 0.693856, 0.61528, 0.633424, 0.806679, 0.55336, 0.601654, 0.560064, 0.5, 0.600977, 0.604798, 0.430853, 1, 0.361529, 0.575604, 0.645187, 0.597788, 0.542484, 0.472388, 1, 0.506642, 0.442787, 0.22472, 0.638169, 0.621614, 0.910914, 0.639699, 0.405993, 0.052778, 0.609838, 0.696969, 0.5, };
