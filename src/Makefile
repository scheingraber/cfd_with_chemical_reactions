CXX:=g++
DEPEND:=$(CXX) -MM

CXXFLAGS:=-std=gnu++0x -c -Wall -pedantic # -O2 -g -Werror
LDFLAGS:=

CXX_TOO_OLD:=$(shell expr `$(CXX) -dumpversion` \< 4.6)

INCLUDES:=-I.

LDIR:=-L/usr/lib
LIBS:= # e.g. -lwt

CXX_SOURCES=$(wildcard *.cpp)
CXX_OBJECTS=$(CXX_SOURCES:.cpp=.o)
CXX_DEPS=$(CXX_OBJECTS:.o=.d)

NODEPS:=clean version_check

.DEFAULT_GOAL:=all

TARGET=sim

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
	-include $(CXX_DEPS)
endif

.PHONY: all clean version_check

version_check:
ifeq ("$(CXX_TOO_OLD)", "1")
	@echo "\nHowdy!\nYou obviously didn't read the README?! Let me help you with that:"
	@/bin/cat README
	@false
endif

all: $(TARGET)

$(TARGET): version_check $(CXX_DEPS) $(CXX_OBJECTS)
	$(CXX) $(LDFLAGS) $(LDIR) $(LIBS) $(CXX_OBJECTS) -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@

clean:
	rm -f $(CXX_OBJECTS) $(CXX_DEPS) $(TARGET) *~

%.d: %.cpp
	$(DEPEND) $< >> $@

