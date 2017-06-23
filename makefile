CC=g++
SOURCES=FitTracks.cpp linearRegressionFit.cpp ResidualHistogrammer.cpp
FASTJETPATH=/localstore/student/wpolarised/branches/branchkees/fastjet313
EXECUTABLE=trackFitter
LDFLAGS=`root-config --glibs`
DEBUGFLAGS=-Wall -g
CFLAGS=`root-config --cflags`
#INCLUDES=


all: $(EXECUTABLE) $(EXECUTABLE2)

#LHEFrunner
$(EXECUTABLE): $(SOURCES) 
	$(CC) -O3 -o $(EXECUTABLE) $(SOURCES) $(CFLAGS) $(DEBUGFLAGS) $(INCLUDES) $(LDFLAGS) -std=c++11


