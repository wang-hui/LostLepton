# AllHadronicSUSY
# LostLepton Background Estimation

1.Set CMS Environment:

setenv SCRAM_ARCH slc6_amd64_gcc530

cmsrel CMSSW_8_0_10

cd CMSSW_8_0_10/src

cmsenv

git clone -b TestMiniAOD git@github.com:susy2015/recipeAUX.git

git clone Ana_June14_2016_fix_event_filter_bugs git@github.com:susy2015/SusyAnaTools.git

cd SusyAnaTools/Tools

git checkout remotes/origin/FromKash samples.cc

git checkout remotes/origin/FromKash samples.h

cd ../..

git clone -b master https://github.com/susy2015/LostLepton.git

scram b -j9

cd LostLepton/Tool

mkdir obj

make


4.Run the LostLepton code:

./LostLepton_MuCS_TTbar runList.txt

5.Make Closure Plots (not tested):

./ClosurePlots

