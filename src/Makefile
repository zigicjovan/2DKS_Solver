CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall
LDFLAGS =
LDLIBS = -lfftw3_mpi -lfftw3 -lm
INCLUDES = -I./dataProcessing -I./problemSolver -I./utilities

TARGET := solver
BUILD_DIR := build

AGGREGATOR := aggregateResults
AGGREGATE_SOURCE := utilities/aggregateResults.cpp
AGGREGATE_INPUT_DIR ?= Data
AGGREGATE_OUTPUT_DIR ?= Data/SolutionBranches

SOURCES := run_2DKS.cpp
SOURCES += $(wildcard dataProcessing/*.cpp)
SOURCES += $(wildcard problemSolver/*.cpp)
SOURCES += $(filter-out $(AGGREGATE_SOURCE),$(wildcard utilities/*.cpp))

OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPENDS := $(OBJECTS:.o=.d)

RUNSCRIPT ?= runscripts/runscript.sh

.PHONY: all run aggregate clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

-include $(DEPENDS)

run: $(TARGET)
	./$(RUNSCRIPT)

$(AGGREGATOR): $(AGGREGATE_SOURCE)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@

aggregate: $(AGGREGATOR)
	./$(AGGREGATOR) $(AGGREGATE_INPUT_DIR) $(AGGREGATE_OUTPUT_DIR)

clean:
	rm -f $(OBJECTS) $(DEPENDS) $(TARGET) $(AGGREGATOR)