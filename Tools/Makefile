ROOTCONFIG   := root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)

CXX           = g++
CXXFLAGS      = -O -Wall -fPIC $(ROOTCFLAGS)
LD            = g++
LDFLAGS       = -O $(ROOTLDFLAGS)
SOFLAGS       = -shared
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

OBJS          = skimtree.o
EXES          = $(OBJS:.o=)

all:    $(EXES)

skimtree: skimtree.o
	$(LD) $(LDFLAGS) $^ $(GLIBS) -lTreePlayer -lGenVector -o $@

.cc.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	@rm -f core $(OBJS) $(EXES)
