_OBJS =	boxchain.o dataset.o funclmdif.o funsupport.o globalVariables.o \
        modelcalculator.o nrutil.o Para.o toad.o toadcmd.o toadmisc.o \
        tvds.o tvDSfit.o tvImg.o tvLinfitDriver.o utable.o \
        fileTools.o

OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

CC = g++
ODIR = obj
SDIR = src
IDIR = inc

CFLAGS = -Wall -O3 -ffast-math -fPIC -lpthread -g
CLINK = -lm -ltcl8.4 -ltk8.4 -ltiff -lBLT24 -lpthread -lmydll \
        -lgsl -lgslcblas -Wl,-rpath,.
CPATH = -I/usr/include/tcl8.4 -I./$(IDIR) -L. -I.

all: libtoad toad

toad: libtoad
	gcc $(ODIR)/toad.o $(CPATH) -ltoad_threaded $(CLINK) -o toad 

libtoad: $(OBJS)
	$(CC) $(OBJS) $(CFLAGS) $(CPATH) $(CLINK) -shared -o libtoad_threaded.so 

$(ODIR)/%.o: $(SDIR)/%.cpp $(IDIR)/%.h
	$(CC) -c $(CPATH) $(CFLAGS) -o $@ $< 

clean:
	rm -f $(ODIR)/*.o libtoad_threaded.so


	

