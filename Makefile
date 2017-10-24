TARGET = TMDB

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
CXXFLAGS   = $(ROOTCFLAGS) -I$(ROOTSYS)/include -O -Wall -fPIC
CXXLIBS    = $(ROOTLIBS)
CC = g++

SRC  = main.C TMDB.C
OBJS = main.o TMDB.o



all: $(TARGET)

$(TARGET) : $(OBJS)
	$(CC) $(OBJS) $(CXXFLAGS) $(CXXLIBS) -o $@ \

.C.o :
	$(CC) $(CXXFLAGS) -c $<

main.o : TMDB.h
TMDB.o : TMDB.h functions_mf.h

clean:
	rm -f $(TARGET) $(OBJS) core.* *~
