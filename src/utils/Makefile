# 
#  Copyright (C) 2011 University of Southern California
#                     and Andrew D. Smith
#
#  Authors: Fang Fang and Andrew D. Smith
# 
#  This software is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of the
#  License, or (at your option) any later version.
# 
#  This software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this software; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
#  02110-1301 USA
ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS = methstates cpg2base sam2mr

SOURCES = $(wildcard *.cpp)

INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas 

CXX = g++
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o MappedRead.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

install: all
	@mkdir -p $(ROOT)/bin
	@install -m 755 $(PROGS) $(ROOT)/bin

%: %.cpp
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

test_%: %
	@$(TEST_DIR)/$@ $(TEST_DIR)

test:   $(addprefix test_, $(PROGS))

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
