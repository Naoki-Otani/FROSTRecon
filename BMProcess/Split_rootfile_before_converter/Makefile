CXX := g++

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)

CXXFLAGS := -O2 -Wall -Wextra -std=c++17 $(ROOT_CFLAGS)
LDFLAGS :=
LDLIBS := $(ROOT_LIBS)

TARGETS := split_daily split_at_time

.PHONY: all clean

all: $(TARGETS)

split_daily: split_daily.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

split_at_time: split_at_time.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGETS) *.o
