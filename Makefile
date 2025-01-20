CC:=g++
ifneq (,$(findstring Darwin,$(shell uname)))
      CC:=g++-9
endif
INCLUDEFLG=-I ./include/
OMPFLG=-fopenmp
HASHFLG=-Wno-deprecated
BUILDFLG=-w -ffunction-sections -fdata-sections -fmodulo-sched -msse
EXE_MIP=bin/PM-Mt-tracker

all:
	$(CC) -o $(EXE_MIP) src/MT_tracker.cpp $(INCLUDEFLG) $(HASHFLG) $(BUILDFLG) $(OMPFLG)

clean:
	rm -rf bin/* src/*.o