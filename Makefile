CC 		= g++

OPT3 	= -O3
CFLAGS 	= -Wall

LDFLAGS = -lm -lquadmath

all: direct_v9-rkp

direct_v7-opt: direct_v9-rkp.cpp
	$(CC) $(CFLAGS) $(EXTRA_CFLAGS) $(OPT3) -fopenmp $+ $(LDFLAGS) -o $@

clean:
	rm -rf direct_v9-rkp *.o *.out *.err *.prv *.pcf *.row *.sym

