CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall
LDFLAGS =
LDLIBS = -lfftw3 -lm
INCLUDES = -I./dataProcessing -I./problemSolver -I./utilities

TARGET = solver
BUILD_DIR := build

SOURCES = run_2DKS.cpp
SOURCES += $(wildcard dataProcessing/*.cpp)
SOURCES += $(wildcard problemSolver/*.cpp)
SOURCES += $(wildcard utilities/*.cpp)

OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPENDS := $(OBJECTS:.o=.d)

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

-include $(DEPENDS)

RUNSCRIPT ?= runscripts/runscript.sh

run: $(TARGET)
	./$(RUNSCRIPT)

clean:
	rm -f $(OBJECTS) $(TARGET)