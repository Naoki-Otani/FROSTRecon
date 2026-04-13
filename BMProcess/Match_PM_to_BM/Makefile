CXX := g++
CXXFLAGS := -O2 -std=c++17 $(shell root-config --cflags)
LDLIBS := $(shell root-config --libs)

TARGET := match_pm_to_bm
SRC := match_pm_to_bm.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(TARGET)
