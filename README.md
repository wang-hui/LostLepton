# AllHadronicSUSY
# LostLepton Background Estimation

1.Set CMS Environment:

setenv SCRAM_ARCH slc6_amd64_gcc491
cmsrel CMSSW_7_4_15
cd CMSSW_7_4_15/src/
cmsenv

git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git
git clone -b Ana_74X_17Nov2015_v3.0 git@github.com:susy2015/SusyAnaTools.git

git clone -b hua_change_structure https://github.com/susy2015/LostLepton.git

scram b -j 6

cd LostLepton/Tool
mkdir obj
make


4.Run the LostLepton code:

./LostLepton_MuCS_TTbar runList.txt 1

5.Make Closure Plots (not tested):

./ClosurePlots

6.Make MtW Plots on TTbar/T2tt samples (not tested):

./MtW runList_ttbar_skimmed_flattree.txt runList_t2tt_425_325.txt runList_t2tt_500_325.txt runList_t2tt_650_325.txt runList_t2tt_850_100.txt 123.root


