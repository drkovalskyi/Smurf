# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else

all: 	help
help:
	@echo "Available Targets:";\
	cat Makefile | perl -ne 'printf("\t%-15s %s\n",$$1,$$2) if(/^(\S+):[^#]+(#.*)$$/)';\
	echo "Example:";\
	echo "\tmake cards plots"
ifndef VERBOSE
  QUIET := @
endif

LANDSPATH = ../../LandS
LANDSLIB = lands.so

LIB = $(LANDSPATH)/$(LANDSLIB)
lands: $(LANDSPATH)/makefile # compile LandS
	$(QUIET) echo "compiling LandS"; \
	ln -s $(LANDSPATH)/$(LANDSLIB); \
	cd $(LANDSPATH); \
	$(MAKE) -f makefile;
data:
	@echo 'Make a directory or a symbolic link called "data" pointing to sample location'; echo; exit 1

cards: 	CutBasedAnalysis.C data # make LandS cards
	$(QUIET) root -b -q CutBasedAnalysis.C+

plots:  $(LANDSLIB) # make plots based on available cards
	$(QUIET) root -l PlotExpectedLimits.C

clean: # clean up 
	$(QUIET) rm -v -f \
	$(LANDSLIB) \
	*.card *.o *.d *.so LinkDef_out*
cleanall: clean # complete clean up
	$(QUIET) rm -v -f \
	limits* *~
	

$(LANDSLIB): lands

CC = g++
CFLAGS = -Wall -Wno-unused-function -O2 -fPIC $(shell root-config --cflags) 

LINKER = g++
LINKERFLAGS = $(shell root-config --ldflags --libs)

ifeq ($(shell root-config --platform),macosx)
	MACOSXFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace
endif

# General rule for making object files
%.d:	%.cc
	$(QUIET) echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.o: 	%.cc 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx 
	$(QUIET) echo "Compiling $<"; \
	$(CC) $(CFLAGS) $< -c -o $@

endif
