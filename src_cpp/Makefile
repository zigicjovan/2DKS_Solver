CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall
LDFLAGS =
LDLIBS = -lfftw3 -lm
INCLUDES = -I./dataProcessing -I./problemSolver -I./utilities

TARGET = solver

SOURCES = run_2DKS.cpp
SOURCES += $(wildcard dataProcessing/*.cpp)
SOURCES += $(wildcard problemSolver/*.cpp)
SOURCES += $(wildcard utilities/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run:
	./runscript.sh

clean:
	rm -f $(OBJECTS) $(TARGET)