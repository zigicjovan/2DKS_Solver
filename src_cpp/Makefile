CXX = g++
CXXFLAGS = -std=c++17 -Wall
INCLUDES = -I./main_functions -I./main_objects -I./solver_functions -I./solver_objects -I./processing_functions

TARGET = solver

SOURCES = run_2DKS.cpp
SOURCES += $(wildcard main_functions/*.cpp)
SOURCES += $(wildcard main_objects/*.cpp)
SOURCES += $(wildcard solver_objects/*.cpp)
SOURCES += $(wildcard solver_functions/*.cpp)
SOURCES += $(wildcard processing_functions/*.cpp)

OBJECTS = $(SOURCES:.cpp=.o)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

run:
	./runscript.sh

clean:
	rm -f $(OBJECTS) $(TARGET)