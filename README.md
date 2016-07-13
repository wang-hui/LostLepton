# AllHadronicSUSY
# LostLepton Background Estimation

1.Set CMS Environment:

setenv SCRAM_ARCH slc6_amd64_gcc530

cmsrel CMSSW_8_0_12

cd CMSSW_8_0_12/src

cmsenv

git cms-merge-topic -u cms-met:CMSSW_8_0_X-METFilterUpdate

git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

git clone -b Ana_V7_new_JEC_new_filters_4ifb git@github.com:susy2015/SusyAnaTools.git

git clone -b master https://github.com/susy2015/LostLepton.git

scram b -j9

cd LostLepton/Tool

mkdir obj

make


4.Run the LostLepton code:

./LostLepton_MuCS_TTbar runList.txt

5.Make Closure Plots (not tested):

./ClosurePlots

