IDIR       = .
ODIR       = obj
SDIR       = ../../SusyAnaTools/Tools/

#CXX       = ~/.vim/bundle/ClangComplete/bin/cc_args.py g++
CXX        = g++

CXXFLAGS  += -Iinc -std=c++0x -I${CMSSW_BASE}/src/
## Optimization flag
CXXFLAGS += -g #-O3
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

LIBS       = $(shell root-config --glibs)

PROGRAMS = tupleTest LostLepton_MuCS_TTbar ClosurePlots SignalScan

all: $(PROGRAMS)

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

tupleTest: $(ODIR)/NTupleReader.o $(ODIR)/tupleReadTest.o
	$(LD) $^ $(LIBS)  -o $@

#makePlots: $(ODIR)/NTupleReader.o $(ODIR)/Plotter.o $(ODIR)/MakePlots.o
#	$(LD) $^ $(LIBS) -L .TupleDict_h.so -o $@

LostLepton_MuCS_TTbar: $(ODIR)/LostLepton_MuCS_TTbar.o $(ODIR)/NTupleReader.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

#LostLepton_ElCS_TTbar: $(ODIR)/LostLepton_ElCS_TTbar.o $(ODIR)/NTupleReader.o
#	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

#LostLepton_MuCS_SingleTop: $(ODIR)/LostLepton_MuCS_SingleTop.o $(ODIR)/NTupleReader.o
#	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

ClosurePlots: $(ODIR)/ClosurePlots.o
	$(LD) $^ $(LIBS)  -o $@

SignalScan:$(ODIR)/SignalScan.o $(ODIR)/NTupleReader.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
