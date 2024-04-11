CC              = g++
CC_FLAGS        = -g3 -O3 -Wall -D_GLIBCXX_DEBUG -I  /opt/homebrew/Cellar/gsl/2.7.1/include/
LD_FLAGS        = -L/opt/homebrew/Cellar/gsl/2.7.1/lib  -lgsl -lgslcblas -lm -lstdc++ 
REC_OBJECTS	= mapping.o io.o utilities.o 

align: $(REC_OBJECTS)
	$(CC) $(CC_FLAGS) $(REC_OBJECTS) -o blanche $(LD_FLAGS)
mapping.o: mapping.cpp
	$(CC) $(CC_FLAGS) -c mapping.cpp
io.o: io.cpp
	$(CC) $(CC_FLAGS) -c io.cpp
utilities.o: utilities.cpp
	$(CC) $(CC_FLAGS) -c utilities.cpp

