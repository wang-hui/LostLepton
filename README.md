# AllHadronicSUSY
# LostLepton Background Estimation

1.Set CMS Environment:

setenv SCRAM_ARCH slc6_amd64_gcc481(export SCRAM_ARCH = slc6_amd64_gcc481)

cmsrel CMSSW_7_2_0

cd CMSSW_7_2_0/src

cmsenv

2.Download source code from github and compile plugins:

git clone -b TestMiniAOD https://github.com/lihux25/recipeAUX.git

git clone -b baseline_def_Jan31_2015 https://github.com/susy2015/SusyAnaTools.git

git clone https://github.com/susy2015/LostLepton.git

scram b -j9

3.Go to LostLepton directory and then compile the code

cd LostLepton/Tool

mkdir obj

make

PS: We have Error like "cannot find libtbb.so.2", to solve this problem, we want to build a soft link to this lib for temporary solution. To do this:

Go to LostLepton/Tool

and then do:

ln -s /cvmfs/cms.cern.ch/slc6_amd64_gcc481/external/tbb/42_20131118oss/lib/libtbb.so.2

This will generate a soft link : libtbb.so.2

And it will be involved into make command automatically and we do not have error report anymore

4.Run the LostLepton code:

./LostLepton_MuCS_TTbar runList_inputfile.txt outputfile.root

./LostLepton_ElCS_TTbar runList_inputfile.txt outputfile.root

5.Make Closure Plots:

./ClosurePlots
