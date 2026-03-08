CXX      ?= g++
CXXFLAGS ?= -O2 -Wall -Wextra -std=c++17

ROOTCLING ?= rootcling
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)
ROOTINCDIR := $(shell root-config --incdir)

INCLUDES := -Iinclude

CXXFLAGS += $(ROOTCFLAGS) $(INCLUDES)
LDFLAGS  += $(ROOTLIBS)

#
# ROOT dictionary (needed for std::vector<std::vector<int>> branches)
#
DICT_NAME := dict
DICT_SRC  := src/$(DICT_NAME).cxx
DICT_OBJ  := src/$(DICT_NAME).o
DICT_PCM  := src/$(DICT_NAME)_rdict.pcm
DICT_PCM_LOCAL := $(DICT_NAME)_rdict.pcm

# Dictionary inputs (human-written)
LINKDEF           := dict/LinkDef.h
DICT_HEADER_NAME  := DictTypes.h           # passed to rootcling (no include/ prefix)
DICT_HEADER_PATH  := include/DictTypes.h   # used only for make dependency tracking

# rootcling must NOT receive -m64/-pthread/etc.  Pass only include paths.
ROOTCLING_INC := -I$(ROOTINCDIR) $(INCLUDES)

SRC_COMMON := \
	src/FROSTConstants.cc \
	src/MultiHitAlgo.cc \
	src/TreeReader.cc

BIN := mapfunction_tool FROST_reconstruction

all: $(BIN)

mapfunction_tool: src/mapfunction_tool.cc $(SRC_COMMON) $(DICT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	@test -f $(DICT_PCM) && cp -f $(DICT_PCM) $(DICT_PCM_LOCAL) || true

FROST_reconstruction: src/FROST_reconstruction.cc $(SRC_COMMON) $(DICT_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	@test -f $(DICT_PCM) && cp -f $(DICT_PCM) $(DICT_PCM_LOCAL) || true

# Generate dictionary source (rootcling will also create dict_rdict.pcm next to dict.cxx)
$(DICT_SRC): $(LINKDEF) $(DICT_HEADER_PATH)
	$(ROOTCLING) -f $(DICT_SRC) $(ROOTCLING_INC) $(DICT_HEADER_NAME) $(LINKDEF)

# Compile dictionary object
$(DICT_OBJ): $(DICT_SRC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(BIN) *.o src/*.o $(DICT_SRC) $(DICT_OBJ) $(DICT_PCM) $(DICT_PCM_LOCAL)

.PHONY: all clean
