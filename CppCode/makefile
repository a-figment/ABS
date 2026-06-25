# Variables to control Makefile operation
#' @source: 1. https://stackoverflow.com/questions/52034997/how-to-make-makefile-recompile-when-a-header-file-is-changed
CXX = g++
CXXFLAGS = -Wall -Wextra -g -Wno-parentheses -Ofast

SOURCES := $(wildcard *.cpp)
OBJECTS := $(patsubst %.cpp,%.o,$(SOURCES))
DEPENDS := $(patsubst %.cpp,%.d,$(SOURCES))

# ****************************************************
# .PHONY means these rules get executed even if
# files of those names exist.
.PHONY: all clean

# The first rule is the default, ie. "make",
# "make all" and "make parking" mean the same
all: main
	
clean:
	$(RM) $(OBJECTS) $(DEPENDS) main

# Linking the executable from the object files
main: $(OBJECTS)
	$(CXX) $(WARNING) $(CXXFLAGS) $^ -o $@

-include $(DEPENDS)

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@


