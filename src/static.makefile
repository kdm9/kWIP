###############################################################################
# Note: This makefile is used by the kWIP developers to compile the static    #
# binaries provided with kWIP. Please use the CMake-based build system unless #
# you know what you're doing. See README.md for install instructions.         #
###############################################################################

# Copyright 2015 Kevin Murray <spam@kdmurray.id.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


CXX 	 := g++
CXXFLAGS += -static -static-libstdc++ -std=c++11 -fopenmp  -DSEQAN_HAS_GZIP=1 -DSEQAN_HAS_BZIP2=1 -I src -I src/ext/
LDFLAGS  += -lz -lbz2

VER=$(shell git describe --always)

OXLI_OBJ=src/ext/oxli/read_aligner.o \
         src/ext/oxli/labelhash.o \
         src/ext/oxli/MurmurHash3.o \
         src/ext/oxli/hashbits.o \
         src/ext/oxli/counting.o \
         src/ext/oxli/kmer_hash.o \
         src/ext/oxli/hashtable.o \
         src/ext/oxli/traversal.o \
         src/ext/oxli/hllcounter.o \
         src/ext/oxli/subset.o \
         src/ext/oxli/read_parsers.o

KWIP_OBJ=src/countmin.o \
		 src/kernel.o \
		 src/kernels/ip.o \
		 src/kernels/wip.o \
		 src/kwip-cli.o \
		 src/kwip-utils.o \
		 src/population.o \
		 $(OXLI_OBJ)


kwip: $(KWIP_OBJ) src/kwip-config.hh
	$(CXX) $(CXXFLAGS) -o kwip $(KWIP_OBJ) $(LDFLAGS)

$(KWIP_OBJ): src/kwip-config.hh

src/kwip-config.hh: src/kwip-config.hh.in
	sed -e 's/^#define KWIP_VERSION.*/#define KWIP_VERSION "$(VER)"/' \
		-e 's/^#cmakedefine ENABLE_MULTITABLE.*/#undef ENABLE_MULTITABLE/' \
		$< >$@

clean:
	rm -f $(KWIP_OBJ) kwip src/kwip-config.hh
