CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
    CC := $(shell brew list --versions gcc | awk '{print $$2}' | cut -d'.' -f1 | awk '{print "g++-"$$1; exit}')
endif
INCLUDEFLG=-I ./include/
OMPFLG=-fopenmp
HASHFLG=-Wno-deprecated
OPTFLAG ?= -O3
BUILDFLG=$(OPTFLAG) -w -ffunction-sections -fdata-sections -fmodulo-sched
EXE_MIP=bin/Mt-tracker

all:
	$(CC) -o $(EXE_MIP) src/MT_tracker.cpp $(INCLUDEFLG) $(HASHFLG) $(BUILDFLG) $(OMPFLG)

clean:
	rm -rf bin/* src/*.o