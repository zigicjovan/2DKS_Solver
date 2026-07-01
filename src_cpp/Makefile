CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall
LDFLAGS =
LDLIBS = -lfftw3 -lm
INCLUDES = -I./mainFunctions -I./mainObjects -I./solverFunctions -I./solverObjects -I./dataFunctions -I./dataObjects

TARGET = solver

SOURCES = run_2DKS.cpp
SOURCES += $(wildcard mainFunctions/*.cpp)
SOURCES += $(wildcard mainObjects/*.cpp)
SOURCES += $(wildcard solverFunctions/*.cpp)
SOURCES += $(wildcard solverObjects/*.cpp)
SOURCES += $(wildcard dataFunctions/*.cpp)
SOURCES += $(wildcard dataObjects/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run:
	./runscript.sh

clean:
	rm -f $(OBJECTS) $(TARGET)