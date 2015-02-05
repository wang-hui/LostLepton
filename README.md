# AllHadronicSUSY
# LostLepton Background Estimation

setenv SCRAM_ARCH slc6_amd64_gcc481(export SCRAM_ARCH = slc6_amd64_gcc481)
cmsrel CMSSW_7_2_0
cd CMSSW_7_2_0/src
cmsenv

git clone -b TestMiniAOD https://github.com/lihux25/recipeAUX.git
git clone -b baseline_def_Jan31_2015 https://github.com/susy2015/SusyAnaTools.git
scram b -j9

Then, you need to copy the lost lepton code.
For now it is in my public area, but soon it will be in git.

cd SusyAnaTools
mkdir LostLepton
cd LostLepton
cp /afs/cern.ch/work/l/lacroix/public/stop_v150204 .
mkdir obj
make

Then, run the code:
./LostLepton
