CXX=g++
CC=gcc
CFLAGS=-O2 -Wall
LDFLAGS=-Llib
PRFFLAGS=-lProof
THRDFLAGS=-lThread
IN=-I$(ROOTSYS)/include/root
INS=-I../../Higgs_CS_and_Width_SM4/include
INSS=-I../../Higgs_CS_and_Width_SM4/src

CFLAGS += `root-config --cflags`
LIBS += `root-config --glibs`


OBJSS=validationTest.o
OBJSSS=validationSM4_Err.o
OBJSSSS=validationTest_scale.o
OBJSSSSS=validationTest_BR.o
OBJ4=validationTest_BRerr.o

.PHONY: clean all main test

all: validate validatePdf validateScale validateBR validationTest_BRerr

validate: validationTest.o
	$(CXX) -o validationTest.exe $(OBJSS) $(LIBS)

validatePdf: validationSM4_Err.o
	$(CXX) -o validationSM4_Err.exe $(OBJSSS) $(LIBS)

validateScale: validationTest_scale.o
	$(CXX) -o validationTest_scale.exe $(OBJSSSS) $(LIBS)

validateBR: validationTest_BR.o
	$(CXX) -o validationTest_BR.exe $(OBJSSSSS) $(LIBS)

validationTest_BRerr: validationTest_BRerr.o
	$(CXX) -o validationTest_BRerr.exe $(OBJ4) $(LIBS)

clean:
	@rm *.o *.exe *~ 

new:
	@rm *.root *.eps *.gif *~

##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(IN) $(INS) $(INSS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(IN) $(INS) $(INSS) -c $<


