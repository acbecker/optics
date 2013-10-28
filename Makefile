#################################################################
CC          = gcc #-g3
CFLAGS      = -Wall -Wno-deprecated -funroll-loops -O3 
CLIBS       = -lm #-lefence 

CPPC        = g++ #-g3
CPPCOPTS    = -c -Wall -Wno-deprecated -funroll-loops -O3 
CPPLINK     = g++ #-g3 
CPPLIBS     = -lm # -lwcs

#################################################################

all: clean run_optics

htm_aux.o: htm_aux.cc
	$(CPPC) $(CPPCOPTS) htm_aux.cc

optics.o: optics.cc 
	$(CPPC) $(CPPCOPTS) optics.cc 

run_optics.o: run_optics.cc 
	$(CPPC) $(CPPCOPTS) run_optics.cc 

run_optics: run_optics.o htm_aux.o optics.o 
	$(CPPLINK) run_optics.o htm_aux.o optics.o $(CPPLIBS) -o run_optics

test_data:
	python maketestdata.py 10 100 1

clean :
	rm -f *.o
	rm -f *~ .*~
	rm -f run_optics
