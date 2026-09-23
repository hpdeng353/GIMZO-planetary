objects = maneos.o #aneos.o
fortran_objects = libaneos.o

FC := gfortran

exe = testlibmaneos testcallaneos maketable

CFLAGS ?= -O3 -march=native
FFLAGS ?= $(CFLAGS)
LIBS ?= -lm -lgfortran

default: testlibmaneos

all: testlibmaneos maketable

# A little program to test the MANEOS wrapper
testlibmaneos: testlibmaneos.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

# A little program to test the function callaneos() of the MANEOS wrapper
testcallaneos: testcallaneos.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

maketable: maketable.o $(objects) $(fortran_objects)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(exe) $(objects) 

cleanall:
	rm -f $(exe) $(objects) $(fortran_objects)
