# Makefile for Screensaver
# Copyright (c) 2012 the Massachusetts Institute of Technology
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# This Makefile is set up to build Screensaver in one of two modes.
# In debug mode, the C compiler does virtually no optimizations
# (-O0) and embeds a bunch of debugging symbols to make work with GDB easier.
# In release mode, the compiler will do every optimization it can (-O3) and keep
# debugging symbols out; this should decrease the size of the final binary and
# make it perform better.
#
# To compile in debug mode, type "make DEBUG=1".  To to compile in release
# mode, type "make DEBUG=0" or simply "make".
#
# If you type "make prof", Make will instrument the output for profiling with
# gprof.  Be sure you run "make clean" first!
#
# If everything gets wacky and you need a sane place to start from, you can
# type "make clean", which will remove all compiled code.
#
# If you want to do something wacky with your compiler flags--like enabling
# debug symbols but keeping optimizations on--you can specify CXXFLAGS or
# LDFLAGS on the command line.  If you want to use a predefined mode but augment
# the predefined CXXFLAGS or LDFLAGS, you can specify EXTRA_CXXFLAGS or
# EXTRA_LDFLAGS on the command line.


# The sources we're building
HEADERS = $(wildcard *.h)
PRODUCT_SOURCES = $(filter-out GraphicStuff.c, $(wildcard *.c))

# What we're building
PRODUCT_OBJECTS = $(PRODUCT_SOURCES:.c=.o)
PRODUCT = Screensaver
PROFILE_PRODUCT = $(PRODUCT:%=%.prof) #the product, instrumented for gprof

# What we're building with
CXX = gcc
CXXFLAGS = -std=gnu99 -Wall -fcilkplus
LDFLAGS = -lrt -lm -lcilkrts


# Determine which profile--debug or release--we should build against, and set
# CFLAGS appropriately.
#
# At this time, we also update the .buildmode stamp, which records if the last
# build was in debug or release mode.  This is a bit hackish--we set all C
# compiler outputs to depend on .buildmode, and then we touch .buildmode if it
# should change.  Touching .buildmode invalidates all the compiler outputs, so
# they all get built again in the correct mode.  Credit to Ceryen Tan and Marek
# Olszewski, who developed this code for 6.197 back in the day.
OLD_MODE = $(shell cat .buildmode 2> /dev/null)
ifeq ($(DEBUG),1)
# We want debug mode.
CXXFLAGS += -g -O0 -gdwarf-3
ifneq ($(OLD_MODE),debug)
$(shell echo debug >.buildmode)
endif
else
# We want release mode.
CXXFLAGS += -O3 -DNDEBUG
ifneq ($(OLD_MODE),release)
$(shell echo release >.buildmode)
endif
endif


# By default, make the product.
all:		$(PRODUCT)

# How to build for profiling
prof:		$(PROFILE_PRODUCT)

# How to clean up
clean:
	$(RM) $(PRODUCT) $(PROFILE_PRODUCT) *.o *.out


# How to compile a C file
%.o:		%.c $(HEADERS) .buildmode
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -o $@ -c $<

# How to link the product
$(PRODUCT): LDFLAGS += -lXext -lX11
$(PRODUCT):	$(PRODUCT_OBJECTS) GraphicStuff.o .buildmode
	$(CXX) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ $(PRODUCT_OBJECTS) GraphicStuff.o

# How to build the product, instrumented for profiling
$(PROFILE_PRODUCT): CXXFLAGS += -DPROFILE_BUILD -pg
$(PROFILE_PRODUCT): LDFLAGS += -pg
$(PROFILE_PRODUCT): $(PRODUCT_OBJECTS) .buildmode
	$(CXX)  $(PRODUCT_OBJECTS) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $(PROFILE_PRODUCT)