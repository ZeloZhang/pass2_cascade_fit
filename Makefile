CC=g++
CFLAGS=-c -g -fPIC -O3 -Wall -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/boost -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/python2.7 `root-config --cflags`

LDFLAGS=-lboost_python -lpython2.7 `root-config --glibs`
LDIR=/data/user/zzhang1/ROOT/build/lib/
LIBS=-lMathMore
SOURCES=main.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard src/bootstrap/*.cpp)
SOURCES_INJECT=inject.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard src/bootstrap/*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=main

OBJECTS_INJECT=$(SOURCES_INJECT:.cpp=.o)
    EXECUTABLE_INJECT=inject

all: $(SOURCES) $(EXECUTABLE)

inject: $(OBJECTS_INJECT)
	   $(CC) $(OBJECTS_INJECT) -o $@ $(LDFLAGS) $(LIBS)

$(EXECUTABLE): $(OBJECTS)
	   $(CC) $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS)

.cpp.o:
	   $(CC) $(CFLAGS) $< -o $@

clean:
	   rm ./*.o ./main src/*.o src/*/*.o
