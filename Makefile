CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: direct_v9-ral

direct_v7-opt: direct_v9-ral.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf direct_v9-ral *.o *.out *.err *.prv *.pcf *.row *.sym

