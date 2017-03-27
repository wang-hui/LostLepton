//in vim, use
//10,87s/\(\d\)\s\+\(\d\)/\1, \2/g
//to convert data card to arries.

void min_and_max(int nsb, double pred[], double unc[])
{
double min = unc[0];
double max = unc[0];
for (int i = 1; i < nsb; i++)
{
if(pred[i] != 0 && unc[i] < min) min = unc[i];
if(pred[i] != 0 && unc[i] > max) max = unc[i];
}
std::cout << " min " << min << " max " << max << std::endl;
}


int Plot_compare_el_cs()
{
const unsigned nsb = 84;
const bool print_min_and_max = false;
const bool print_plot = true;

//A: TF method
double PredA[nsb] = {558.050, 22.829, 7.546, 0.948, 1.280, 313.443, 18.883, 3.719, 0.000, 36.908, 23.434, 8.786, 1.159, 0.000, 2.037, 2.686, 5.918, 0.000, 0.000, 1.377, 0.000, 232.412, 16.532, 2.585, 0.000, 1.002, 18.641, 2.530, 1.170, 0.000, 2.030, 5.578, 1.185, 1.447, 0.000, 0.000, 0.000, 34.602, 5.754, 1.782, 1.259, 1.498, 3.148, 0.000, 0.832, 2.184, 0.588, 0.588, 17.721, 2.641, 0.548, 0.565, 12.811, 5.259, 0.480, 0.550, 0.000, 0.000, 16.082, 3.329, 0.313, 0.505, 9.894, 1.390, 1.292, 0.495, 0.767, 0.000, 0.000, 0.000, 2.404, 0.000, 0.000, 1.354, 0.000, 1.392, 1.043, 0.000, 0.000, 0.589, 0.285, 0.000, 0.000, 0.000};

double stat_unc_abs_up[nsb] = {22.491, 4.975, 3.030, 1.251, 1.688, 16.405, 4.695, 2.516, 1.460, 4.960, 5.010, 3.031, 1.528, 0.720, 1.378, 1.817, 2.525, 0.738, 1.208, 1.816, 0.917, 14.572, 4.452, 2.044, 2.282, 2.304, 4.231, 2.461, 1.544, 1.119, 1.605, 2.547, 1.153, 1.408, 0.710, 2.172, 0.667, 6.292, 3.437, 2.351, 2.895, 1.976, 2.489, 1.151, 1.913, 2.124, 1.353, 1.352, 3.136, 1.302, 0.723, 1.300, 3.922, 2.244, 1.104, 1.264, 1.646, 0.643, 2.884, 1.520, 0.719, 1.161, 3.271, 1.834, 1.705, 1.139, 1.763, 1.158, 0.969, 2.737, 1.626, 1.559, 0.335, 1.786, 0.512, 1.355, 1.376, 0.954, 1.295, 1.356, 0.656, 0.998, 0.993, 1.220};

double ErrorUpA[nsb];
double ErrorDownA[nsb];

double ErrorUpA1[nsb] = {0.040, 0.218, 0.402, 1.319, 1.319, 0.052, 0.249, 0.677, 0.000, 0.134, 0.214, 0.345, 1.319, 0.000, 0.677, 0.677, 0.427, 0.000, 0.000, 1.319, 0.000, 0.063, 0.269, 0.791, 0.000, 2.300, 0.227, 0.973, 1.319, 0.000, 0.791, 0.457, 0.973, 0.973, 0.000, 0.000, 0.000, 0.182, 0.597, 1.319, 2.300, 1.319, 0.791, 0.000, 2.300, 0.973, 2.300, 2.300, 0.177, 0.493, 1.319, 2.300, 0.306, 0.427, 2.300, 2.300, 0.000, 0.000, 0.179, 0.457, 2.300, 2.300, 0.331, 1.319, 1.319, 2.300, 2.300, 0.000, 0.000, 0.000, 0.677, 0.000, 0.000, 1.319, 0.000, 0.973, 1.319, 0.000, 0.000, 2.300, 2.300, 0.000, 0.000, 0.000};
double ErrorDownA1[nsb] = {0.039, 0.182, 0.297, 0.646, 0.646, 0.050, 0.203, 0.432, 0.000, 0.119, 0.179, 0.264, 0.646, 0.000, 0.432, 0.432, 0.311, 0.000, 0.000, 0.646, 0.000, 0.059, 0.216, 0.479, 0.000, 0.827, 0.188, 0.544, 0.646, 0.000, 0.479, 0.327, 0.544, 0.544, 0.000, 0.000, 0.000, 0.156, 0.397, 0.646, 0.827, 0.646, 0.479, 0.000, 0.827, 0.544, 0.827, 0.827, 0.152, 0.346, 0.646, 0.827, 0.240, 0.311, 0.827, 0.827, 0.000, 0.000, 0.154, 0.327, 0.827, 0.827, 0.255, 0.646, 0.646, 0.827, 0.827, 0.000, 0.000, 0.000, 0.432, 0.000, 0.000, 0.646, 0.000, 0.544, 0.646, 0.000, 0.000, 0.827, 0.827, 0.000, 0.000, 0.000};

double ErrorUpA2[nsb] = {0.017, 0.073, 0.143, 0.260, 0.244, 0.022, 0.073, 0.141, 0.218, 0.058, 0.067, 0.135, 0.207, 0.262, 0.244, 0.159, 0.169, 0.353, 0.439, 0.459, 0.299, 0.025, 0.086, 0.170, 0.250, 0.395, 0.075, 0.161, 0.297, 0.417, 0.232, 0.193, 0.294, 0.364, 0.507, 0.512, 0.638, 0.067, 0.144, 0.304, 0.620, 0.100, 0.170, 0.258, 0.385, 0.076, 0.150, 0.324, 0.069, 0.161, 0.298, 0.410, 0.097, 0.165, 0.255, 0.263, 0.332, 0.542, 0.072, 0.164, 0.307, 0.384, 0.113, 0.205, 0.341, 0.219, 0.480, 0.281, 0.610, 0.469, 0.220, 0.432, 0.479, 0.217, 0.354, 0.165, 0.170, 0.204, 0.378, 0.241, 0.335, 0.301, 0.717, 0.378};
double ErrorDownA2[nsb] = {0.017, 0.073, 0.143, 0.260, 0.244, 0.022, 0.073, 0.141, 0.218, 0.058, 0.067, 0.135, 0.207, 0.262, 0.244, 0.159, 0.169, 0.353, 0.439, 0.459, 0.299, 0.025, 0.086, 0.170, 0.250, 0.395, 0.075, 0.161, 0.297, 0.417, 0.232, 0.193, 0.294, 0.364, 0.507, 0.512, 0.638, 0.067, 0.144, 0.304, 0.620, 0.100, 0.170, 0.258, 0.385, 0.076, 0.150, 0.324, 0.069, 0.161, 0.298, 0.410, 0.097, 0.165, 0.255, 0.263, 0.332, 0.542, 0.072, 0.164, 0.307, 0.384, 0.113, 0.205, 0.341, 0.219, 0.480, 0.281, 0.610, 0.469, 0.220, 0.432, 0.479, 0.217, 0.354, 0.165, 0.170, 0.204, 0.378, 0.241, 0.335, 0.301, 0.717, 0.378};

double ErrorUpA3[nsb] = {0.012, 0.008, 0.020, 0.006, 0.006, 0.016, 0.018, 0.011, 0.012, 0.015, 0.008, 0.031, 0.000, 0.034, 0.011, 0.024, 0.004, 0.033, 0.008, 0.031, 0.029, 0.014, 0.013, 0.050, 0.024, 0.068, 0.013, 0.012, 0.017, 0.034, 0.002, 0.018, 0.012, 0.034, 0.005, 0.068, 0.030, 0.007, 0.012, 0.013, 0.044, 0.022, 0.020, 0.016, 0.009, 0.022, 0.029, 0.017, 0.059, 0.042, 0.021, 0.075, 0.025, 0.043, 0.058, 0.000, 0.039, 0.001, 0.035, 0.035, 0.044, 0.012, 0.013, 0.026, 0.065, 0.013, 0.054, 0.031, 0.055, 0.002, 0.036, 0.003, 0.014, 0.056, 0.008, 0.049, 0.020, 0.021, 0.033, 0.004, 0.004, 0.049, 0.052, 0.054};
double ErrorDownA3[nsb] = {0.018, 0.014, 0.025, 0.008, 0.005, 0.023, 0.025, 0.015, 0.013, 0.021, 0.010, 0.053, 0.006, 0.032, 0.013, 0.033, 0.003, 0.047, 0.029, 0.017, 0.033, 0.018, 0.019, 0.068, 0.026, 0.077, 0.020, 0.026, 0.026, 0.032, 0.005, 0.022, 0.013, 0.055, 0.000, 0.076, 0.030, 0.010, 0.018, 0.022, 0.042, 0.033, 0.032, 0.024, 0.027, 0.037, 0.052, 0.039, 0.086, 0.069, 0.048, 0.133, 0.038, 0.068, 0.093, 0.015, 0.076, 0.001, 0.051, 0.058, 0.067, 0.000, 0.021, 0.036, 0.093, 0.021, 0.087, 0.052, 0.088, 0.016, 0.054, 0.006, 0.021, 0.096, 0.019, 0.098, 0.009, 0.046, 0.076, 0.010, 0.053, 0.103, 0.092, 0.130};

double ErrorUpA4[nsb] = {0.000, 0.000, 0.003, 0.010, 0.004, 0.000, 0.003, 0.004, 0.004, 0.002, 0.001, 0.001, 0.006, 0.000, 0.006, 0.002, 0.002, 0.002, 0.003, 0.004, 0.002, 0.000, 0.000, 0.004, 0.007, 0.001, 0.000, 0.004, 0.004, 0.006, 0.013, 0.007, 0.005, 0.004, 0.003, 0.002, 0.002, 0.002, 0.004, 0.003, 0.010, 0.004, 0.002, 0.001, 0.006, 0.003, 0.009, 0.003, 0.002, 0.000, 0.010, 0.010, 0.003, 0.000, 0.000, 0.002, 0.003, 0.012, 0.001, 0.002, 0.001, 0.003, 0.000, 0.001, 0.004, 0.001, 0.010, 0.000, 0.005, 0.021, 0.004, 0.006, 0.009, 0.010, 0.003, 0.008, 0.009, 0.009, 0.001, 0.001, 0.008, 0.006, 0.033, 0.014};
double ErrorDownA4[nsb] = {0.039, 0.182, 0.297, 0.646, 0.646, 0.050, 0.203, 0.432, 0.004, 0.119, 0.179, 0.264, 0.646, 0.000, 0.432, 0.432, 0.311, 0.002, 0.004, 0.646, 0.002, 0.059, 0.216, 0.479, 0.008, 0.827, 0.188, 0.544, 0.646, 0.006, 0.479, 0.327, 0.544, 0.544, 0.003, 0.003, 0.000, 0.156, 0.397, 0.646, 0.827, 0.646, 0.479, 0.002, 0.827, 0.544, 0.827, 0.827, 0.152, 0.346, 0.646, 0.827, 0.240, 0.311, 0.827, 0.827, 0.003, 0.012, 0.154, 0.327, 0.827, 0.827, 0.255, 0.646, 0.646, 0.827, 0.827, 0.001, 0.005, 0.021, 0.432, 0.006, 0.009, 0.646, 0.003, 0.544, 0.646, 0.010, 0.001, 0.827, 0.827, 0.006, 0.036, 0.015};

double ErrorUpA5[nsb] = {0.001, 0.003, 0.373, 0.008, 0.042, 0.001, 0.002, 0.019, 0.002, 0.009, 0.003, 0.002, 0.029, 0.016, 0.040, 0.017, 0.001, 0.022, 0.087, 0.020, 0.072, 0.001, 0.013, 0.007, 0.025, 0.021, 0.000, 0.000, 0.010, 0.064, 0.128, 0.005, 0.009, 0.067, 0.013, 0.008, 0.023, 0.004, 0.008, 0.029, 0.041, 0.000, 0.009, 0.008, 0.006, 0.014, 0.028, 0.007, 0.002, 0.010, 0.005, 0.021, 0.001, 0.000, 0.048, 0.012, 0.012, 0.017, 0.008, 0.004, 0.005, 0.034, 0.004, 0.012, 0.006, 0.001, 0.003, 0.027, 0.012, 0.056, 0.001, 0.028, 0.038, 0.006, 0.017, 0.006, 0.004, 0.003, 0.019, 0.016, 0.011, 0.013, 0.021, 0.029};
double ErrorDownA5[nsb] = {0.003, 0.023, 0.434, 0.052, 0.049, 0.001, 0.005, 0.022, 0.010, 0.016, 0.003, 0.005, 0.011, 0.026, 0.082, 0.007, 0.008, 0.000, 0.069, 0.019, 0.064, 0.002, 0.018, 0.003, 0.011, 0.038, 0.007, 0.002, 0.008, 0.093, 0.477, 0.003, 0.007, 0.043, 0.020, 0.020, 0.002, 0.006, 0.006, 0.018, 0.058, 0.002, 0.020, 0.003, 0.017, 0.013, 0.037, 0.085, 0.002, 0.014, 0.027, 0.032, 0.001, 0.001, 0.001, 0.017, 0.010, 0.015, 0.005, 0.002, 0.016, 0.016, 0.003, 0.010, 0.023, 0.007, 0.014, 0.056, 0.010, 0.045, 0.000, 0.031, 0.058, 0.001, 0.014, 0.010, 0.005, 0.030, 0.120, 0.004, 0.010, 0.041, 0.008, 0.025};

double ErrorUpA6[nsb] = {0.003, 0.007, 0.092, 0.074, 0.000, 0.001, 0.058, 0.081, 0.155, 0.034, 0.063, 0.067, 0.001, 0.507, 0.036, 0.080, 0.002, 0.217, 0.403, 0.125, 0.755, 0.030, 0.006, 0.009, 0.065, 0.000, 0.008, 0.149, 0.059, 0.189, 0.167, 0.116, 0.012, 0.022, 0.455, 0.019, 0.386, 0.011, 0.063, 0.112, 0.546, 0.007, 0.012, 0.062, 0.044, 0.037, 0.029, 0.011, 0.029, 0.031, 0.066, 0.400, 0.025, 0.019, 0.048, 0.047, 0.153, 0.075, 0.015, 0.017, 0.107, 0.040, 0.011, 0.048, 0.175, 0.082, 0.073, 0.000, 0.060, 0.095, 0.091, 0.068, 0.010, 0.005, 0.100, 0.004, 0.000, 0.167, 0.153, 0.038, 0.069, 0.002, 0.315, 0.485};
double ErrorDownA6[nsb] = {0.000, 0.025, 0.003, 0.023, 0.099, 0.004, 0.047, 0.126, 0.105, 0.036, 0.003, 0.016, 0.032, 0.047, 0.095, 0.047, 0.078, 0.066, 0.231, 0.115, 0.025, 0.000, 0.004, 0.061, 0.093, 0.146, 0.012, 0.076, 0.072, 0.147, 0.011, 0.034, 0.072, 0.035, 0.051, 0.138, 0.310, 0.012, 0.008, 0.062, 0.067, 0.085, 0.145, 0.180, 0.063, 0.075, 0.068, 0.007, 0.018, 0.004, 0.004, 0.017, 0.000, 0.037, 0.167, 0.018, 0.063, 0.011, 0.007, 0.003, 0.115, 0.113, 0.026, 0.033, 0.171, 0.047, 0.009, 0.012, 0.005, 0.001, 0.019, 0.246, 0.130, 0.087, 0.003, 0.017, 0.017, 0.009, 0.022, 0.035, 0.122, 0.013, 0.045, 0.001};

double ErrorUpA7[nsb] = {0.003, 0.008, 0.011, 0.006, 0.005, 0.004, 0.012, 0.003, 0.017, 0.008, 0.008, 0.003, 0.003, 0.053, 0.003, 0.000, 0.006, 0.004, 0.000, 0.001, 0.004, 0.003, 0.028, 0.042, 0.055, 0.104, 0.006, 0.000, 0.074, 0.036, 0.002, 0.033, 0.004, 0.033, 0.014, 0.013, 0.002, 0.002, 0.017, 0.000, 0.000, 0.017, 0.034, 0.004, 0.018, 0.013, 0.017, 0.003, 0.010, 0.078, 0.046, 0.013, 0.020, 0.035, 0.005, 0.032, 0.005, 0.004, 0.020, 0.067, 0.050, 0.004, 0.010, 0.024, 0.006, 0.066, 0.123, 0.002, 0.030, 0.151, 0.046, 0.008, 0.000, 0.000, 0.003, 0.012, 0.007, 0.009, 0.002, 0.023, 0.023, 0.015, 0.031, 0.001};
double ErrorDownA7[nsb] = {0.002, 0.006, 0.024, 0.005, 0.010, 0.002, 0.014, 0.022, 0.011, 0.005, 0.005, 0.001, 0.003, 0.050, 0.057, 0.004, 0.012, 0.002, 0.008, 0.025, 0.000, 0.002, 0.003, 0.031, 0.015, 0.018, 0.023, 0.014, 0.003, 0.006, 0.009, 0.006, 0.006, 0.000, 0.039, 0.000, 0.000, 0.018, 0.007, 0.045, 0.041, 0.027, 0.012, 0.013, 0.007, 0.000, 0.012, 0.029, 0.004, 0.039, 0.001, 0.012, 0.005, 0.004, 0.000, 0.044, 0.002, 0.002, 0.003, 0.003, 0.051, 0.013, 0.019, 0.027, 0.000, 0.034, 0.008, 0.010, 0.048, 0.033, 0.003, 0.000, 0.006, 0.054, 0.000, 0.033, 0.013, 0.009, 0.006, 0.018, 0.004, 0.018, 0.019, 0.000};

double ErrorUpA8[nsb] = {0.005, 0.009, 0.050, 0.053, 0.015, 0.018, 0.014, 0.044, 0.030, 0.011, 0.004, 0.021, 0.036, 0.000, 0.104, 0.066, 0.035, 0.143, 0.003, 0.424, 0.046, 0.020, 0.047, 0.008, 0.003, 0.102, 0.033, 0.052, 0.053, 0.008, 0.028, 0.108, 0.005, 0.083, 0.276, 0.026, 0.002, 0.033, 0.047, 0.070, 0.014, 0.074, 0.015, 0.007, 0.072, 0.002, 0.041, 0.040, 0.018, 0.065, 0.002, 0.066, 0.007, 0.136, 0.176, 0.000, 0.247, 0.084, 0.009, 0.092, 0.122, 0.324, 0.012, 0.026, 0.263, 0.186, 0.001, 0.554, 0.030, 0.214, 0.049, 0.164, 0.000, 0.163, 0.201, 0.002, 0.027, 0.140, 0.348, 0.089, 0.115, 0.219, 0.075, 0.312};
double ErrorDownA8[nsb] = {0.000, 0.004, 0.011, 0.025, 0.107, 0.006, 0.078, 0.011, 0.021, 0.028, 0.053, 0.060, 0.057, 0.010, 0.041, 0.005, 0.006, 0.045, 0.217, 0.222, 0.053, 0.000, 0.030, 0.072, 0.069, 0.029, 0.005, 0.017, 0.114, 0.043, 0.044, 0.000, 0.020, 0.053, 0.643, 0.161, 0.090, 0.026, 0.060, 0.090, 0.174, 0.018, 0.023, 0.010, 0.005, 0.002, 0.024, 0.036, 0.056, 0.033, 0.104, 0.102, 0.021, 0.030, 0.068, 0.035, 0.127, 0.070, 0.094, 0.046, 0.000, 0.134, 0.041, 0.214, 0.150, 0.063, 0.546, 0.031, 0.002, 0.324, 0.007, 0.301, 0.103, 0.029, 0.071, 0.118, 0.232, 0.105, 0.231, 0.368, 0.234, 0.410, 0.065, 0.261};

double ErrorUpA9[nsb] = {0.143, 0.166, 0.138, 0.146, 0.178, 0.146, 0.135, 0.131, 0.128, 0.175, 0.138, 0.178, 0.167, 0.259, 0.211, 0.144, 0.132, 0.132, 0.188, 0.253, 0.128, 0.153, 0.148, 0.148, 0.127, 0.105, 0.161, 0.120, 0.166, 0.105, 0.165, 0.157, 0.221, 0.149, 0.131, 0.107, 0.234, 0.147, 0.140, 0.134, 0.158, 0.151, 0.154, 0.137, 0.144, 0.165, 0.154, 0.149, 0.196, 0.242, 0.191, 0.146, 0.137, 0.151, 0.146, 0.147, 0.087, 0.227, 0.218, 0.226, 0.187, 0.164, 0.170, 0.122, 0.158, 0.179, 0.141, 0.312, 0.150, 0.110, 0.212, 0.220, 0.542, 0.134, 0.212, 0.183, 0.263, 0.203, 0.162, 0.351, 0.292, 0.269, 0.135, 0.542};
double ErrorDownA9[nsb] = {0.146, 0.170, 0.141, 0.149, 0.183, 0.149, 0.138, 0.134, 0.131, 0.179, 0.142, 0.182, 0.171, 0.264, 0.216, 0.146, 0.134, 0.135, 0.193, 0.258, 0.131, 0.156, 0.152, 0.151, 0.130, 0.108, 0.165, 0.123, 0.170, 0.108, 0.168, 0.160, 0.225, 0.152, 0.133, 0.109, 0.239, 0.151, 0.143, 0.139, 0.162, 0.154, 0.157, 0.140, 0.147, 0.168, 0.158, 0.153, 0.200, 0.247, 0.194, 0.149, 0.140, 0.154, 0.148, 0.150, 0.089, 0.235, 0.222, 0.230, 0.190, 0.167, 0.174, 0.125, 0.164, 0.184, 0.145, 0.319, 0.153, 0.113, 0.216, 0.224, 0.553, 0.136, 0.216, 0.186, 0.268, 0.207, 0.164, 0.357, 0.297, 0.274, 0.137, 0.551};

//B: LL method

double PredB[nsb] = { 526.62236, 20.68731, 7.25702, 0.77669, 0.54157, 315.17141, 15.47174, 3.99243, 0.00000, 36.86093, 21.51228, 7.38420, 2.35528, 0.00000, 1.83546, 2.77413, 5.45928, 0.00000, 0.00000, 1.58572, 0.00000, 206.07879, 20.65984, 1.51328, 0.00000, 0.53091, 18.25767, 2.06294, 0.88196, 0.00000, 2.59554, 5.12647, 0.76886, 1.35641, 0.00000, 0.00000, 0.00000, 31.04311, 3.71694, 0.98639, 1.87903, 1.29103, 2.12926, 0.00000, 6.03990, 1.45210, 0.43629, 0.64102, 18.64841, 3.66814, 0.73710, 0.43403, 11.13384, 6.16763, 2.10920, 0.59264, 0.00000, 0.00000, 15.62546, 2.66099, 0.25198, 0.59319, 8.56282, 1.07309, 2.37665, 0.46946, 0.47615, 0.00000, 0.00000, 0.00000, 1.88514, 0.00000, 0.00000, 1.13073, 0.00000, 1.46772, 0.73784, 0.00000, 0.00000, 0.36582, 0.18131, 0.00000, 0.00000, 0.00000};

double avg_weight[nsb] = { 0.76111, 0.64307, 0.81815, 0.49388, 0.40288, 0.71423, 0.71977, 0.86432, 0.52828, 0.54922, 0.75247, 0.67931, 0.60780, 0.33887, 0.42287, 0.52772, 0.60815, 0.55598, 0.42422, 0.68920, 0.51671, 0.70016, 0.63413, 0.55442, 0.97103, 0.62896, 0.67011, 1.01647, 0.61673, 0.52730, 0.59042, 0.56515, 0.33646, 0.43478, 0.78056, 0.82275, 0.47275, 0.81923, 0.68846, 0.72558, 0.78916, 0.58734, 0.61198, 0.68876, 0.76163, 0.57087, 0.66599, 0.67437, 0.38373, 0.46969, 0.32441, 0.60692, 0.83025, 0.72363, 0.51538, 0.62168, 0.92336, 0.67919, 0.35023, 0.44408, 0.33970, 0.45095, 0.63113, 0.64244, 0.79089, 0.59042, 0.67512, 0.60941, 0.39625, 0.84927, 0.46653, 0.58492, 0.72392, 0.64514, 0.51673, 0.37487, 0.35070, 0.41767, 0.48407, 0.38646, 0.32989, 0.35159, 0.63034, 0.31140};

double ErrorUpB[nsb];
double ErrorDownB[nsb];

double ErrorUpB1[nsb] = { 0.05350, 0.21347, 0.33285, 0.73657, 0.74079, 0.06554, 0.22010, 0.55641, 0, 0.13818, 0.20029, 0.28590, 0.77184, 0, 0.48135, 0.49539, 0.33785, 0, 0, 0.75106, 0, 0.07458, 0.36880, 0.52127, 0, 1.04167, 0.21751, 0.60356, 0.73884, 0, 0.54674, 0.41029, 0.60267, 0.64461, 0, 0, 0, 0.18530, 0.44082, 0.73660, 1.04167, 0.79740, 0.53986, 0, 1.04167, 0.63484, 1.04167, 1.04167, 0.30495, 0.49749, 0.77756, 1.04167, 0.26300, 0.38300, 1.04167, 1.04167, 0, 0, 0.22610, 0.35265, 1.04167, 1.04167, 0.29576, 0.74301, 0.80136, 1.04167, 1.04167, 0, 0, 0, 0.47755, 0, 0, 0.74174, 0, 0.62369, 0.77277, 0, 0, 1.04167, 1.04167, 0, 0, 0};
double ErrorDownB1[nsb] = {0.05350, 0.21347, 0.33285, 0.73657, 0.74079, 0.06554, 0.22010, 0.55641, 0, 0.13818, 0.20029, 0.28590, 0.77184, 0, 0.48135, 0.49539, 0.33785, 0, 0, 0.75106, 0, 0.07458, 0.36880, 0.52127, 0, 1.04167, 0.21751, 0.60356, 0.73884, 0, 0.54674, 0.41029, 0.60267, 0.64461, 0, 0, 0, 0.18530, 0.44082, 0.73660, 1.04167, 0.79740, 0.53986, 0, 1.04167, 0.63484, 1.04167, 1.04167, 0.30495, 0.49749, 0.77756, 1.04167, 0.26300, 0.38300, 1.04167, 1.04167, 0, 0, 0.22610, 0.35265, 1.04167, 1.04167, 0.29576, 0.74301, 0.80136, 1.04167, 1.04167, 0, 0, 0, 0.47755, 0, 0, 0.74174, 0, 0.62369, 0.77277, 0, 0, 1.04167, 1.04167, 0, 0, 0};

double ErrorUpB2[nsb] = {0.02073, 0.09617, 0.25641, 0.23057, 0.35951, 0.04089, 0.07752, 0.22027, 0.25513, 0.20219, 0.10053, 0.15934, 0.19737, 0.28755, 0.30243, 0.15508, 0.17243, 0.40644, 0.63064, 0.44433, 0.28389, 0.03812, 0.10217, 0.18897, 0.36205, 0.75959, 0.13298, 0.30300, 0.25246, 0.41516, 0.32896, 0.22224, 0.28075, 0.33454, 0.30536, 0.86442, 0.47589, 0.11879, 0.18478, 0.37981, 1.03504, 0.13189, 0.22425, 0.22159, 0.39820, 0.10027, 0.26070, 0.28581, 0.08007, 0.43189, 0.31919, 0.33933, 0.16067, 0.33070, 0.23517, 0.23880, 0.27172, 0.62006, 0.07730, 0.26496, 0.25383, 0.38634, 0.12671, 0.23120, 0.37685, 0.29213, 0.42935, 0.23780, 0.39901, 0.77004, 0.20332, 0.65448, 0.84516, 0.24589, 0.60945, 0.21948, 0.25552, 0.22983, 0.57522, 0.30539, 0.41698, 0.36452, 0.59243, 0.89166};
double ErrorDownB2[nsb] = {0.02073, 0.09617, 0.25641, 0.23057, 0.35951, 0.04089, 0.07752, 0.22027, 0.25513, 0.20219, 0.10053, 0.15934, 0.19737, 0.28755, 0.30243, 0.15508, 0.17243, 0.40644, 0.63064, 0.44433, 0.28389, 0.03812, 0.10217, 0.18897, 0.36205, 0.75959, 0.13298, 0.30300, 0.25246, 0.41516, 0.32896, 0.22224, 0.28075, 0.33454, 0.30536, 0.86442, 0.47589, 0.11879, 0.18478, 0.37981, 1.03504, 0.13189, 0.22425, 0.22159, 0.39820, 0.10027, 0.26070, 0.28581, 0.08007, 0.43189, 0.31919, 0.33933, 0.16067, 0.33070, 0.23517, 0.23880, 0.27172, 0.62006, 0.07730, 0.26496, 0.25383, 0.38634, 0.12671, 0.23120, 0.37685, 0.29213, 0.42935, 0.23780, 0.39901, 0.77004, 0.20332, 0.65448, 0.84516, 0.24589, 0.60945, 0.21948, 0.25552, 0.22983, 0.57522, 0.30539, 0.41698, 0.36452, 0.59243, 0.89166};

double ErrorUpB3[nsb] = {0.00291, 0.00282, 0.00274, 0.00337, 0.00202, 0.00298, 0.00274, 0.00262, 0.00251, 0.00277, 0.00290, 0.00284, 0.00325, 0.00308, 0.00314, 0.00262, 0.00269, 0.00201, 0.00331, 0.00297, 0.00365, 0.00296, 0.00272, 0.00252, 0.00264, 0.00265, 0.00288, 0.00305, 0.00218, 0.00166, 0.00253, 0.00235, 0.00388, 0.00314, 0.00190, 0.00144, 0.00245, 0.00295, 0.00306, 0.00257, 0.00249, 0.00265, 0.00341, 0.00251, 0.00391, 0.00280, 0.00261, 0.00259, 0.00244, 0.00215, 0.00171, 0.00292, 0.00272, 0.00289, 0.00247, 0.00300, 0.00294, 0.00267, 0.00241, 0.00244, 0.00172, 0.00227, 0.00287, 0.00299, 0.00282, 0.00276, 0.00204, 0.00198, 0.00242, 0.00275, 0.00259, 0.00285, 0.00241, 0.00279, 0.00251, 0.00251, 0.00197, 0.00267, 0.00243, 0.00218, 0.00157, 0.00163, 0.00135, 0.00180};
double ErrorDownB3[nsb] = {0.00291, 0.00282, 0.00274, 0.00337, 0.00202, 0.00298, 0.00274, 0.00262, 0.00251, 0.00277, 0.00290, 0.00284, 0.00325, 0.00308, 0.00314, 0.00262, 0.00269, 0.00201, 0.00331, 0.00297, 0.00365, 0.00296, 0.00272, 0.00252, 0.00264, 0.00265, 0.00288, 0.00305, 0.00218, 0.00166, 0.00253, 0.00235, 0.00388, 0.00314, 0.00190, 0.00144, 0.00245, 0.00295, 0.00306, 0.00257, 0.00249, 0.00265, 0.00341, 0.00251, 0.00391, 0.00280, 0.00261, 0.00259, 0.00244, 0.00215, 0.00171, 0.00292, 0.00272, 0.00289, 0.00247, 0.00300, 0.00294, 0.00267, 0.00241, 0.00244, 0.00172, 0.00227, 0.00287, 0.00299, 0.00282, 0.00276, 0.00204, 0.00198, 0.00242, 0.00275, 0.00259, 0.00285, 0.00241, 0.00279, 0.00251, 0.00251, 0.00197, 0.00267, 0.00243, 0.00218, 0.00157, 0.00163, 0.00135, 0.00180};

double ErrorUpB4[nsb] = {0.01002, 0.01023, 0.01043, 0.00888, 0.01219, 0.00984, 0.01044, 0.01072, 0.01099, 0.01036, 0.01004, 0.01020, 0.00919, 0.00961, 0.00945, 0.01073, 0.01056, 0.01223, 0.00903, 0.00987, 0.00819, 0.00990, 0.01048, 0.01097, 0.01068, 0.01066, 0.01009, 0.00968, 0.01181, 0.01308, 0.01093, 0.01139, 0.00764, 0.00945, 0.01249, 0.01362, 0.01113, 0.00991, 0.00964, 0.01084, 0.01105, 0.01064, 0.00878, 0.01100, 0.00756, 0.01029, 0.01074, 0.01081, 0.01117, 0.01189, 0.01296, 0.01000, 0.01048, 0.01005, 0.01108, 0.00980, 0.00993, 0.01060, 0.01123, 0.01116, 0.01292, 0.01159, 0.01012, 0.00982, 0.01024, 0.01038, 0.01215, 0.01230, 0.01122, 0.01041, 0.01079, 0.01015, 0.01124, 0.01031, 0.01100, 0.01100, 0.01232, 0.01059, 0.01118, 0.01180, 0.01329, 0.01315, 0.01383, 0.01274};
double ErrorDownB4[nsb] = {0.01002, 0.01023, 0.01043, 0.00888, 0.01219, 0.00984, 0.01044, 0.01072, 0.01099, 0.01036, 0.01004, 0.01020, 0.00919, 0.00961, 0.00945, 0.01073, 0.01056, 0.01223, 0.00903, 0.00987, 0.00819, 0.00990, 0.01048, 0.01097, 0.01068, 0.01066, 0.01009, 0.00968, 0.01181, 0.01308, 0.01093, 0.01139, 0.00764, 0.00945, 0.01249, 0.01362, 0.01113, 0.00991, 0.00964, 0.01084, 0.01105, 0.01064, 0.00878, 0.01100, 0.00756, 0.01029, 0.01074, 0.01081, 0.01117, 0.01189, 0.01296, 0.01000, 0.01048, 0.01005, 0.01108, 0.00980, 0.00993, 0.01060, 0.01123, 0.01116, 0.01292, 0.01159, 0.01012, 0.00982, 0.01024, 0.01038, 0.01215, 0.01230, 0.01122, 0.01041, 0.01079, 0.01015, 0.01124, 0.01031, 0.01100, 0.01100, 0.01232, 0.01059, 0.01118, 0.01180, 0.01329, 0.01315, 0.01383, 0.01274};

double ErrorUpB5[nsb] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double ErrorDownB5[nsb] = {0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002};

double ErrorUpB6[nsb] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};
double ErrorDownB6[nsb] = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

double ErrorUpB7[nsb] = {0.06751, 0.08337, 0.21559, 0.46934, 0.34960, 0.04372, 0.08707, 0.26661, 0.33029, 0.13711, 0.13148, 0.13874, 0.21323, 0.93903, 0.33359, 0.24116, 0.20059, 0.52455, 0.82548, 0.54332, 0.45935, 0.05460, 0.11698, 0.22028, 0.27946, 0.72325, 0.15775, 0.21446, 0.37248, 0.46241, 0.45733, 0.20872, 0.22677, 0.62925, 0.81001, 0.83513, 8.00558, 0.08784, 0.14011, 0.24900, 0.41354, 0.14565, 0.23303, 0.43913, 0.96204, 0.28562, 0.53460, 0.40129, 0.08431, 0.24341, 0.96759, 0.44133, 0.10439, 0.20168, 0.51342, 0.39421, 0.28315, 0.80266, 0.08138, 0.35211, 0.34805, 0.54350, 0.11967, 0.26434, 0.38034, 0.17772, 0.36286, 1.47511, 0.58133, 0.50255, 0.17379, 0.48122, 0.70439, 0.28882, 0.32052, 0.32035, 0.73136, 1.00000, 0.45013, 0.51762, 0.40119, 0.48862, 0.80247, 1.19207};
double ErrorDownB7[nsb] = {0.06751, 0.08337, 0.21559, 0.46934, 0.34960, 0.04372, 0.08707, 0.26661, 0.33029, 0.13711, 0.13148, 0.13874, 0.21323, 0.93903, 0.33359, 0.24116, 0.20059, 0.52455, 0.82548, 0.54332, 0.45935, 0.05460, 0.11698, 0.22028, 0.27946, 0.72325, 0.15775, 0.21446, 0.37248, 0.46241, 0.45733, 0.20872, 0.22677, 0.62925, 0.81001, 0.83513, 8.00558, 0.08784, 0.14011, 0.24900, 0.41354, 0.14565, 0.23303, 0.43913, 0.96204, 0.28562, 0.53460, 0.40129, 0.08431, 0.24341, 0.96759, 0.44133, 0.10439, 0.20168, 0.51342, 0.39421, 0.28315, 0.80266, 0.08138, 0.35211, 0.34805, 0.54350, 0.11967, 0.26434, 0.38034, 0.17772, 0.36286, 1.47511, 0.58133, 0.50255, 0.17379, 0.48122, 0.70439, 0.28882, 0.32052, 0.32035, 0.73136, 1.00000, 0.45013, 0.51762, 0.40119, 0.48862, 0.80247, 1.19207};

double ErrorUpB8[nsb] = {0.01113, 0.01326, 0.01288, 0.01026, 0.01600, 0.01092, 0.01229, 0.01557, 0.01355, 0.01195, 0.01072, 0.01128, 0.01038, 0.01455, 0.01504, 0.01374, 0.01234, 0.01684, 0.00645, 0.00967, 0.01073, 0.01162, 0.01335, 0.01870, 0.00944, 0.01735, 0.01283, 0.00943, 0.01587, 0.01847, 0.01158, 0.01350, 0.00756, 0.00893, 0.01584, 0.02457, 0.01795, 0.01104, 0.01001, 0.01326, 0.01163, 0.01577, 0.00943, 0.01405, 0.00946, 0.01343, 0.01858, 0.02564, 0.01775, 0.02219, 0.01606, 0.01769, 0.01218, 0.01344, 0.01553, 0.01619, 0.01310, 0.02221, 0.01987, 0.02000, 0.02256, 0.02424, 0.01183, 0.01394, 0.01221, 0.01180, 0.01255, 0.01958, 0.01268, 0.01083, 0.01381, 0.01231, 0.02112, 0.01802, 0.01863, 0.01927, 0.02284, 0.01350, 0.01842, 0.02143, 0.02068, 0.01760, 0.01776, 0.02007};
double ErrorDownB8[nsb] = {0.01247, 0.01462, 0.01453, 0.01396, 0.01967, 0.01268, 0.01374, 0.01758, 0.01577, 0.01546, 0.01320, 0.01441, 0.01107, 0.02082, 0.01942, 0.01703, 0.01402, 0.01949, 0.01651, 0.01299, 0.01263, 0.01306, 0.01497, 0.02000, 0.01199, 0.02890, 0.01464, 0.01240, 0.01872, 0.01912, 0.01479, 0.01634, 0.00961, 0.01252, 0.02442, 0.03269, 0.02620, 0.01234, 0.01168, 0.01463, 0.01514, 0.01849, 0.01114, 0.01823, 0.01191, 0.01937, 0.02133, 0.02715, 0.02102, 0.02604, 0.02059, 0.02671, 0.01405, 0.01594, 0.01749, 0.01879, 0.01523, 0.02354, 0.02310, 0.02256, 0.02797, 0.02566, 0.01342, 0.01641, 0.01534, 0.01413, 0.01333, 0.02133, 0.01664, 0.01419, 0.01562, 0.01447, 0.02228, 0.02332, 0.02263, 0.02170, 0.02976, 0.02602, 0.02253, 0.02525, 0.02256, 0.02249, 0.02007, 0.02355};

double ErrorUpB9[nsb] = {0.05184, 0.05929, 0.05626, 0.04671, 0.05378, 0.05121, 0.05241, 0.05303, 0.05257, 0.04938, 0.05030, 0.04352, 0.05966, 0.04556, 0.05445, 0.05137, 0.05656, 0.05464, 0.03254, 0.05636, 0.04845, 0.05339, 0.05757, 0.05811, 0.03957, 0.05905, 0.05926, 0.04645, 0.05660, 0.04956, 0.04264, 0.06198, 0.04001, 0.04282, 0.05831, 0.06474, 0.06777, 0.05355, 0.04719, 0.06398, 0.05104, 0.05967, 0.06133, 0.05065, 0.04181, 0.04498, 0.07056, 0.06911, 0.06349, 0.06663, 0.04959, 0.06665, 0.04908, 0.05485, 0.07566, 0.05943, 0.05337, 0.06637, 0.06562, 0.06759, 0.06131, 0.07128, 0.05592, 0.06236, 0.04560, 0.05081, 0.05533, 0.07978, 0.04875, 0.04505, 0.05567, 0.04949, 0.07107, 0.06207, 0.06967, 0.06428, 0.05978, 0.06456, 0.06101, 0.06630, 0.05821, 0.05918, 0.05848, 0.05682};
double ErrorDownB9[nsb] = {0.05827, 0.06661, 0.06277, 0.06195, 0.06008, 0.05886, 0.05882, 0.05901, 0.05756, 0.05942, 0.05936, 0.05210, 0.06623, 0.06062, 0.06431, 0.06120, 0.06208, 0.06298, 0.05953, 0.07019, 0.05624, 0.06018, 0.06485, 0.06341, 0.04561, 0.07746, 0.06830, 0.05241, 0.06661, 0.05201, 0.05257, 0.07286, 0.04849, 0.05367, 0.07676, 0.07925, 0.08836, 0.06063, 0.05371, 0.07260, 0.05723, 0.06875, 0.07244, 0.05660, 0.04993, 0.05778, 0.08150, 0.07444, 0.07274, 0.07410, 0.05701, 0.07316, 0.05591, 0.06140, 0.08848, 0.06665, 0.06289, 0.07279, 0.07391, 0.07596, 0.07182, 0.07867, 0.06439, 0.07320, 0.05478, 0.05891, 0.06004, 0.09031, 0.06060, 0.05389, 0.06263, 0.05519, 0.07544, 0.07093, 0.08276, 0.07200, 0.07445, 0.06950, 0.06787, 0.07506, 0.06374, 0.06250, 0.06209, 0.06223};

double ErrorUpB10[nsb] = {0.02404, 0.02699, 0.02675, 0.02862, 0.04028, 0.02486, 0.02658, 0.03328, 0.03211, 0.03245, 0.02687, 0.02792, 0.02221, 0.05079, 0.04105, 0.03584, 0.03064, 0.04060, 0.05043, 0.03248, 0.02526, 0.02503, 0.02909, 0.03718, 0.02625, 0.06870, 0.02674, 0.02762, 0.03457, 0.03622, 0.03129, 0.03190, 0.02014, 0.02640, 0.06495, 0.06201, 0.05337, 0.02313, 0.02197, 0.02675, 0.02990, 0.03612, 0.02095, 0.03676, 0.02445, 0.04483, 0.03811, 0.04474, 0.04039, 0.04846, 0.04085, 0.06384, 0.02772, 0.03307, 0.02894, 0.03364, 0.02699, 0.04277, 0.04324, 0.04187, 0.05943, 0.04920, 0.02550, 0.03025, 0.03566, 0.02768, 0.02695, 0.04944, 0.03809, 0.03362, 0.02830, 0.02823, 0.04467, 0.04921, 0.04331, 0.04182, 0.06517, 0.06985, 0.04763, 0.05426, 0.03857, 0.05540, 0.04135, 0.03714};
double ErrorDownB10[nsb] = {0.03209, 0.03562, 0.03544, 0.03752, 0.05349, 0.03330, 0.03561, 0.04590, 0.04327, 0.04158, 0.03569, 0.04057, 0.02727, 0.05893, 0.04965, 0.04517, 0.03772, 0.05372, 0.05228, 0.03384, 0.03440, 0.03312, 0.03714, 0.05021, 0.03457, 0.08719, 0.03561, 0.03395, 0.04336, 0.05407, 0.03716, 0.04194, 0.02769, 0.03495, 0.07119, 0.09349, 0.07201, 0.03122, 0.03013, 0.03415, 0.04126, 0.04586, 0.02638, 0.04851, 0.03217, 0.05418, 0.04292, 0.05672, 0.05490, 0.06857, 0.05715, 0.07535, 0.03742, 0.04365, 0.03926, 0.04718, 0.03499, 0.04802, 0.05978, 0.05528, 0.07775, 0.06535, 0.03346, 0.04188, 0.04305, 0.03567, 0.03209, 0.05752, 0.04749, 0.04033, 0.04019, 0.03890, 0.06041, 0.06199, 0.05588, 0.05344, 0.07312, 0.08058, 0.05992, 0.07180, 0.05364, 0.06224, 0.05129, 0.04539};

double ErrorUpB11[nsb] = {0.07226, 0.07621, 0.07171, 0.07654, 0.07860, 0.07361, 0.07415, 0.07416, 0.07594, 0.07371, 0.07392, 0.06920, 0.09037, 0.08141, 0.07955, 0.07602, 0.08533, 0.08301, 0.08712, 0.09590, 0.07202, 0.07367, 0.07795, 0.07841, 0.06092, 0.10510, 0.07587, 0.06851, 0.06828, 0.06611, 0.06732, 0.08425, 0.06587, 0.07345, 0.10706, 0.10716, 0.11542, 0.07150, 0.06869, 0.07894, 0.07549, 0.08319, 0.07830, 0.07086, 0.06226, 0.07665, 0.08001, 0.08025, 0.08466, 0.09086, 0.07473, 0.09789, 0.07265, 0.08219, 0.08798, 0.08156, 0.06926, 0.07955, 0.08714, 0.08606, 0.09270, 0.10937, 0.07727, 0.08615, 0.07501, 0.07127, 0.08125, 0.13069, 0.08291, 0.07161, 0.07507, 0.07096, 0.09923, 0.08771, 0.09869, 0.07791, 0.08731, 0.09547, 0.08645, 0.10189, 0.06933, 0.08714, 0.08141, 0.06289};
double ErrorDownB11[nsb] = {0.07959, 0.08408, 0.07864, 0.08333, 0.08489, 0.08099, 0.08153, 0.07985, 0.08295, 0.08005, 0.08122, 0.07451, 0.10253, 0.08599, 0.08488, 0.08215, 0.09500, 0.09064, 0.09420, 0.10896, 0.07777, 0.08102, 0.08603, 0.08431, 0.06574, 0.11179, 0.08349, 0.07425, 0.07313, 0.07027, 0.07291, 0.09412, 0.07103, 0.08038, 0.11816, 0.11431, 0.12616, 0.07884, 0.07516, 0.08840, 0.08293, 0.09104, 0.08773, 0.07569, 0.06571, 0.08134, 0.08644, 0.08483, 0.09200, 0.09736, 0.08033, 0.10348, 0.07941, 0.08991, 0.09934, 0.08808, 0.07481, 0.08555, 0.09408, 0.09339, 0.09910, 0.11998, 0.08566, 0.09527, 0.08121, 0.07777, 0.09217, 0.14833, 0.09139, 0.07771, 0.08242, 0.07691, 0.10938, 0.09342, 0.10859, 0.08297, 0.09227, 0.10075, 0.09280, 0.11040, 0.07386, 0.09414, 0.08995, 0.06653};

double ErrorUpB12[nsb] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
double ErrorDownB12[nsb] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};

if(print_plot)
{
double ratio[nsb];
double ratio_error[nsb];
double ErrorA[nsb];
double ErrorB[nsb];
double sys_unc_up[nsb];
double sys_unc_down[nsb];

TCanvas *canv = new TCanvas();
//canv -> SetLogy();
Double_t x[nsb];
Double_t exl[nsb];
Double_t exh[nsb];
for (int i = 0; i < nsb; i++)
{
  x[i]= i + 0.5;
  exl[i] = 0.5;
  exh[i] = 0.5;

  if (PredA[i] == 0) ErrorUpA[i] = stat_unc_abs_up[i];
  else ErrorUpA[i] = sqrt(ErrorUpA1[i]*ErrorUpA1[i] + ErrorUpA2[i]*ErrorUpA2[i] + ErrorUpA3[i]*ErrorUpA3[i] + ErrorUpA4[i]*ErrorUpA4[i] + ErrorUpA5[i]*ErrorUpA5[i] + ErrorUpA6[i]*ErrorUpA6[i] + ErrorUpA7[i]*ErrorUpA7[i] + ErrorUpA8[i]*ErrorUpA8[i] + ErrorUpA9[i]*ErrorUpA9[i]) * PredA[i];
  
  ErrorDownA[i] = sqrt(ErrorDownA1[i]*ErrorDownA1[i] + ErrorDownA2[i]*ErrorDownA2[i] + ErrorDownA3[i]*ErrorDownA3[i] + ErrorDownA4[i]*ErrorDownA4[i] + ErrorDownA5[i]*ErrorDownA5[i] + ErrorDownA6[i]*ErrorDownA6[i] + ErrorDownA7[i]*ErrorDownA7[i] + ErrorDownA8[i]*ErrorDownA8[i] + ErrorDownA9[i]*ErrorDownA9[i]) * PredA[i];

  sys_unc_up[i] = sqrt(ErrorUpB2[i]*ErrorUpB2[i] + ErrorUpB3[i]*ErrorUpB3[i] + ErrorUpB4[i]*ErrorUpB4[i] + ErrorUpB5[i]*ErrorUpB5[i] + ErrorUpB6[i]*ErrorUpB6[i] + ErrorUpB7[i]*ErrorUpB7[i] + ErrorUpB8[i]*ErrorUpB8[i] + ErrorUpB9[i]*ErrorUpB9[i] + ErrorUpB10[i]*ErrorUpB10[i] + ErrorUpB11[i]*ErrorUpB11[i] + ErrorUpB12[i]*ErrorUpB12[i]);

  sys_unc_down[i] = sqrt(ErrorDownB2[i]*ErrorDownB2[i] + ErrorDownB4[i]*ErrorDownB4[i] + ErrorDownB4[i]*ErrorDownB4[i] + ErrorDownB5[i]*ErrorDownB5[i] + ErrorDownB6[i]*ErrorDownB6[i] + ErrorDownB7[i]*ErrorDownB7[i] + ErrorDownB8[i]*ErrorDownB8[i] + ErrorDownB9[i]*ErrorDownB9[i] + ErrorDownB10[i]*ErrorDownB10[i] + ErrorDownB11[i]*ErrorDownB11[i] + ErrorDownB12[i]*ErrorDownB12[i]);

  if (PredB[i] == 0) ErrorUpB[i] = 1.8 * avg_weight[i];
  else ErrorUpB[i] = sqrt(ErrorUpB1[i]*ErrorUpB1[i] + sys_unc_up[i]*sys_unc_up[i]) * PredB[i];
  
  ErrorDownB[i] = sqrt(ErrorDownB1[i]*ErrorDownB1[i] + sys_unc_down[i]*sys_unc_down[i]) * PredB[i];
  
  ErrorA[i] = (ErrorUpA[i] + ErrorDownA[i])/2;
  ErrorB[i] = (ErrorUpB[i] + ErrorDownB[i])/2;

  if (PredA[i] == 0 && PredB[i] == 0)
  {ratio[i] = 1; ratio_error[i] = 0;}

  else if (PredA[i] == 0 && PredB[i] != 0)
  {ratio[i] = 1; ratio_error[i] = ErrorB[i]/PredB[i];}

  else if (PredA[i] != 0 && PredB[i] == 0)
  {ratio[i] = 1; ratio_error[i] = ErrorA[i]/PredA[i];}

  else
  {ratio[i] = PredA[i] / PredB[i];
  ratio_error[i] = ratio[i] * sqrt(ErrorA[i]*ErrorA[i]/PredA[i]/PredA[i] + ErrorB[i]*ErrorB[i]/PredB[i]/PredB[i]);}
}

  std::cout << "double psystup[" << nsb << "] = {";
  for (int i = 0; i < nsb; i++)
  {
  std::cout << sys_unc_up[i];
  if (i != nsb-1) std::cout << ", ";
  if (i == nsb-1) std::cout << "}" << std::endl;
  }

  std::cout << "double psystdown[" << nsb << "] = {";
  for (int i = 0; i < nsb; i++)
  {
  std::cout << sys_unc_down[i];
  if (i != nsb-1) std::cout << ", ";
  if (i == nsb-1) std::cout << "}" << std::endl;
  }

TPad *padup = new TPad("padup", "padup", 0, 0.3, 1, 1.0);
padup -> SetBottomMargin(0); // Upper and lower plot are joined                                                                      
padup -> Draw();             // Draw the upper pad: pad1                                                                    
padup -> cd();

TGraphAsymmErrors* A = new TGraphAsymmErrors(nsb,x,PredA,exl,exh,ErrorDownA,ErrorUpA);
A->SetTitle("");
A->SetLineWidth(2);
A->SetLineColor(kBlue);
//A->SetMarkerStyle(21);
A->GetXaxis()->SetTitle("Search region bin number");
A->GetXaxis()->SetTitleSize(0.045);
A->GetYaxis()->SetTitle("Events");
A->GetYaxis()->SetTitleSize(0.045);
A->GetXaxis()->SetRangeUser(0,84);
A->GetYaxis()->SetRangeUser(0,700);
A -> Draw("AP0");

TGraphAsymmErrors* B = new TGraphAsymmErrors(nsb,x,PredB,exl,exh,ErrorDownB,ErrorUpB);
B -> SetLineWidth(2);
B -> SetLineColor(kRed);
B -> Draw("SAME P0");

TLegend* leg = new TLegend(0.75,0.7,0.9,0.85);
leg->SetBorderSize(0);
leg->SetTextFont(42);
leg->SetFillColor(0);
leg->AddEntry(A,"TF method","l");
leg->AddEntry(B,"LL method","l");
leg->Draw("same");

padup->SetLogy();

canv -> cd();

TPad *paddown = new TPad("paddown", "paddown", 0, 0, 1, 0.3);
paddown -> SetTopMargin(0);
paddown -> SetBottomMargin(0.3);
paddown -> Draw();
paddown -> cd();

TGraphErrors* C = new TGraphErrors(nsb, x, ratio, exl, ratio_error);
C->SetTitle("");
C->SetLineWidth(2);
C->SetLineColor(kBlue);
//C->SetMarkerStyle(7);
C->GetXaxis()->SetTitle("Search region bin number");
C->GetXaxis()->SetTitleSize(0.1);
C->GetXaxis()->SetLabelSize(0.1);
C->GetXaxis()->SetRangeUser(0,84);
C->GetYaxis()->SetTitle("TF / LL");
C->GetYaxis()->SetTitleOffset(0.4);
C->GetYaxis()->SetTitleSize(0.1);
C->GetYaxis()->SetRangeUser(0,2);

C -> Draw("AP0");

TLine* line = new TLine(0, 1, 84, 1);
line->SetLineColor(kRed);
line->Draw("SAME");

canv->SaveAs( "DataCardCampare_0_700_el_cs_log.pdf" );
}

if(print_min_and_max)
{
std::cout << "closure ";
min_and_max(nsb, PredB, ErrorUpB2);

std::cout << "dimuon ";
min_and_max(nsb, PredB, ErrorUpB3);

std::cout << "diele ";
min_and_max(nsb, PredB, ErrorUpB4);

std::cout << "purity ";
min_and_max(nsb, PredB, ErrorUpB5);

std::cout << "mt ";
min_and_max(nsb, PredB, ErrorUpB6);

std::cout << "acc ";
min_and_max(nsb, PredB, ErrorUpB7);

std::cout << "muiso ";
min_and_max(nsb, PredB, ErrorUpB8);

std::cout << "eiso ";
min_and_max(nsb, PredB, ErrorUpB9);

std::cout << "mureco ";
min_and_max(nsb, PredB, ErrorUpB10);

std::cout << "ereco ";
min_and_max(nsb, PredB, ErrorUpB11);

std::cout << "isotrk ";
min_and_max(nsb, PredB, ErrorUpB12);
}

return 0;
}

