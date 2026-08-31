# -----------------------------------------------------------------------------
# Select exactly one cluster. Keep this consistent with CLUSTER in
# param_driver.py. It may also be overridden on the command line, for example:
#     make CLUSTER=anvil
# -----------------------------------------------------------------------------
CLUSTER = nibi
# CLUSTER = anvil
# CLUSTER = bridges2
# CLUSTER = stampede3

# mpicxx is available on all four systems; mpic++ is absent on Stampede3.
CXX = mpicxx

ARCHFLAGS =
FFTW_INCLUDES =
FFTW_LDFLAGS =

ifeq ($(CLUSTER),nibi)
    ARCHFLAGS = -march=native
else ifeq ($(CLUSTER),anvil)
    ARCHFLAGS = -march=znver3
else ifeq ($(CLUSTER),bridges2)
    ARCHFLAGS = -march=znver2
else ifeq ($(CLUSTER),stampede3)
    ARCHFLAGS = $(TACC_VEC_FLAGS)
    FFTW_INCLUDES = -I$(TACC_FFTW3_INC)
    FFTW_LDFLAGS = -L$(TACC_FFTW3_LIB) -Wl,-rpath,$(TACC_FFTW3_LIB)
else
    $(error Unsupported CLUSTER '$(CLUSTER)'; use nibi, anvil, bridges2, or stampede3)
endif

CXXFLAGS = -std=c++17 -O3 $(ARCHFLAGS) -Wall
LDFLAGS = $(FFTW_LDFLAGS)
LDLIBS = -lfftw3_mpi -lfftw3 -lm

SRC_DIR := src
BUILD_DIR := $(SRC_DIR)/build

INCLUDES = $(FFTW_INCLUDES) \
           -I$(SRC_DIR)/dataProcessing \
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
