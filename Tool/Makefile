IDIR       = .
ODIR       = obj
SDIR       = ../../SusyAnaTools/Tools/

#CXX       = ~/.vim/bundle/ClangComplete/bin/cc_args.py g++
CXX        = g++

CXXFLAGS  += -Iinc -std=c++0x -I${CMSSW_BASE}/src/ -I$(LHAPDF_DATA_PATH)/../../include
## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

OPENCV_DIRECTORY=$(CMSSW_BASE)/src/opencv
LIBS       = $(shell root-config --glibs)
LIBS      += -L$(CMSSW_BASE)/src/SusyAnaTools/Tools/obj -lSusyAnaTools
LIBS      += -L$(CMSSW_BASE)/src/TopTagger/TopTagger/test -lTopTagger
LIBS      += -L$(OPENCV_DIRECTORY)/lib/ -lopencv_ml -lopencv_core
LIBS      += -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2
LIBS      += -L$(CMSSW_RELEASE_BASE)/lib/${SCRAM_ARCH}/ -lPhysicsToolsHeppy

LHAPDFLIB = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF

PROGRAMS = tupleTest LostLepton_MuCS_TTbar ClosurePlots TagAndProbeAnalysis MtW

all: mkobj $(PROGRAMS)

mkobj:
	@mkdir -p obj

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

tupleTest: $(ODIR)/NTupleReader.o $(ODIR)/tupleReadTest.o
	$(LD) $^ $(LIBS) $(LHAPDFLIB) -o $@

#makePlots: $(ODIR)/NTupleReader.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o
#	$(LD) $^ $(LIBS) -L .TupleDict_h.so -o $@

LostLepton_MuCS_TTbar: $(ODIR)/LostLepton_MuCS_TTbar.o $(ODIR)/NTupleReader.o $(ODIR)/searchBins.o $(ODIR)/PDFUncertainty.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 $(LHAPDFLIB) -o $@

ClosurePlots: $(ODIR)/ClosurePlots.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS)  -o $@

TagAndProbeAnalysis:$(ODIR)/TagAndProbeAnalysis.o $(ODIR)/NTupleReader.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

MtW: $(ODIR)/MtW.o $(ODIR)/NTupleReader.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
