# Makefile for mppc_correction
# Usage:
#   make
#   make clean

CXX      := g++
CXXFLAGS := -O2 -std=c++17 -Wall -Wextra -Wshadow -Wformat=2
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

TARGET := mppc_correction
SRCS   := mppc_correction.cpp
OBJS   := $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $^ $(ROOTLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
