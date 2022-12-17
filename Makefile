#CC=g++
#CC=`root-config --cxx`
CC=/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/bin/gcc
#CFLAGS=-c -g -fPIC -O3 -Wall -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/boost -I/cvmfs/icecube.opensciencegrid.org/py2-v3/RHEL_7_x86_64/include/python2.7 `root-config --cflags`
#CFLAGS=-c -g -fPIC -O3 -Wall -I/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/include/boost -I/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/include/python3.6m/ `root-config --cflags`
CFLAGS=-c -g -fPIC -O3 -Wall -I/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/include/boost `root-config --cflags`

LDIR1=/data/user/zzhang1/ROOT/build/lib/
LDIR2=/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/RHEL_7_x86_64/lib/
#LDFLAGS=-lboost_python36 -lpython3.6m -lstdc++ `root-config --glibs` -L$(LDIR2)
LDFLAGS=-lstdc++ `root-config --glibs` -L$(LDIR2)
#LIBS=-lMathMore
LIBS=-lMathMore -lASImage	
SOURCES=main.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard src/bootstrap/*.cpp)
SOURCES_INJECT=inject.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard src/bootstrap/*.cpp)

SOURCES_PROFILE=profile_scan.cpp $(wildcard src/*.cpp) $(wildcard src/systematics/*.cpp) $(wildcard src/models/*.cpp) $(wildcard sr    c/bootstrap/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o)
	EXECUTABLE=main

OBJECTS_INJECT=$(SOURCES_INJECT:.cpp=.o)
    EXECUTABLE_INJECT=inject

OBJECTS_PROFILE=$(SOURCES_PROFILE:.cpp=.o)
    EXECUTABLE_PROFILE=profile_scan

all: $(SOURCES) $(EXECUTABLE)

inject: $(OBJECTS_INJECT)
	$(CC) -pg $(OBJECTS_INJECT) -o $@ $(LDFLAGS) $(LIBS)

profile_scan: $(OBJECTS_PROFILE)
	$(CC) -pg $(OBJECTS_PROFILE) -o $@ $(LDFLAGS) $(LIBS)

$(EXECUTABLE): $(OBJECTS)
	   $(CC) -pg $(OBJECTS) -o $@ $(LDFLAGS) $(LIBS) 

.cpp.o:
	   $(CC) -pg $(CFLAGS) $< -o $@

clean:
	   rm ./*.o ./main src/*.o src/*/*.o
