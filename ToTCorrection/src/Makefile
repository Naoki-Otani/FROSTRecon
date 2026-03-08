CXX      := g++
CXXFLAGS := -O2 -Wall -Wextra -std=c++17 $(shell root-config --cflags)
LDFLAGS  := $(shell root-config --libs)

TARGET := ToTCorrection
SRC    := ToTCorrection.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f $(TARGET)
