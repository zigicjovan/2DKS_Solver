CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall
LDFLAGS =
LDLIBS = -lfftw3_mpi -lfftw3 -lm

SRC_DIR := src
BUILD_DIR := $(SRC_DIR)/build

INCLUDES = -I$(SRC_DIR)/dataProcessing \
           -I$(SRC_DIR)/problemSolver \
           -I$(SRC_DIR)/utilities

TARGET := solver

AGGREGATOR := aggregateResults
AGGREGATE_SOURCE := $(SRC_DIR)/utilities/aggregateResults.cpp
AGGREGATE_INPUT_DIR ?= Data
AGGREGATE_OUTPUT_DIR ?= Data/SolutionBranches

SOURCES := $(SRC_DIR)/run_2DKS.cpp
SOURCES += $(wildcard $(SRC_DIR)/dataProcessing/*.cpp)
SOURCES += $(wildcard $(SRC_DIR)/problemSolver/*.cpp)
SOURCES += $(filter-out $(AGGREGATE_SOURCE),$(wildcard $(SRC_DIR)/utilities/*.cpp))

OBJECTS := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPENDS := $(OBJECTS:.o=.d)

RUNSCRIPT ?= runscripts/runscript.sh

.PHONY: all run aggregate clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(TARGET) $(LDFLAGS) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
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