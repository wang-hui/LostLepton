# AllHadronicSUSY
# LostLepton Background Estimation

1.Set CMS Environment:

setenv SCRAM_ARCH slc6_amd64_gcc491

cmsrel CMSSW_7_4_15

cd CMSSW_7_4_15/src/

cmsenv

git cms-merge-topic -u kpedro88:METfix7415

git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

git clone -b Ana_mergeBins_LE3_74X_12Feb2016_v5.0_PreApproval git@github.com:susy2015/SusyAnaTools.git

git clone -b v160224_postpreapp https://github.com/susy2015/LostLepton.git

scram b -j 4

cd LostLepton/Tool

mkdir obj

make


4.Run the LostLepton code:

./LostLepton_MuCS_TTbar runList.txt

5.Make Closure Plots (not tested):

./ClosurePlots

