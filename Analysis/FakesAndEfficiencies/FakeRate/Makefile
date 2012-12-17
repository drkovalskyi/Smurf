CC = g++
ROOFITINCLUDE = $(shell scramv1 tool info roofitcore | grep INCLUDE | sed 's/^INCLUDE=/-I/')
INCLUDE = -I../ -I./ $(ROOFITINCLUDE)
CFLAGS = -Wall -g -fPIC $(shell root-config --cflags) $(INCLUDE) $(EXTRACFLAGS) -DTOOLSLIB
LINKER = g++

LINKERFLAGS = $(shell root-config --ldflags)
ifeq ($(shell root-config --platform),macosx)
	LINKERFLAGS = -dynamiclib -undefined dynamic_lookup -Wl,-x -O -Xlinker -bind_at_load -flat_namespace $(shell root-config --libs) -lEG -lGenVector
endif

SOURCES = FakeLooper.cc SmurfPlotUtilities.cc ProcessFakeRate.cc
OBJECTS = $(SOURCES:.cc=.o) LinkDef_out.o
LIB = libSmurfFakeLooper.so

$(LIB):	$(OBJECTS) 
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) -o $@ 

LIBS = $(LIB)

LinkDef_out.cxx: LinkDef.h FakeLooper.h SmurfPlotUtilities.h ProcessFakeRate.h
	rootcint -f $@ -c $(INCLUDE)  FakeLooper.h SmurfPlotUtilities.h ProcessFakeRate.h $<

# General rule for making object files
%.d:	%.cc
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@
%.d:	%.cxx
	$(CC) -MM -MT $@ -MT ${@:.d=.o} $(CFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

%.o: 	%.cc 
	$(CC) $(CFLAGS) $< -c -o $@

%.o: 	%.cxx
	$(CC) $(CFLAGS) $< -c -o $@

.PHONY: clean all dy
all:  $(LIBS)

fakes:  $(LIBS) 
	root -b -q doFakeRate.C

clean:  
	rm -f *.d \
	rm -f *.o \
	rm -f *.so \
	rm -f *.cxx \

-include $(SOURCES:.cc=.d)
-include $(LIBDIR)/LinkDef_out.d

