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


GPGKEY 		= "E56C E72C 8060 C9F9 5F6B  5E72 0FEB 966A 205B 4780"
HBB_DOCKER 	= docker run -t -i --rm -v $(shell pwd):/io \
			  		phusion/holy-build-box-64:latest bash

DEP_TARS_	= libbz2.tar.gz
DEP_TARS 	= $(foreach dep,$(DEP_TARS_),deps/$(dep))

VERSION		= $(shell git describe)
PREFIXDIR	= kwip_$(VERSION)_amd64
TARS		= $(PREFIXDIR).tar.gz
TARSUMS 	= $(foreach tar,$(TARS),$(tar).sha512sums)
SIGS   		= $(foreach tar,$(TARS),$(tar).asc)

.PHONY: all clean cleandep sign
all: $(TARS) $(TARSUMS)

sign: $(SIGS)

clean:
	rm -rf kwip_*_amd64*

cleandep:
	rm -rf deps

$(PREFIXDIR): src/hbb-build.sh $(DEP_TARS) .git/index
	$(HBB_DOCKER) /io/$< $(VERSION)


%.sha512sums: %
	sha512sum $< >$@

%.asc: %
	gpg --armour --local-user $(GPGKEY) --detach-sign $<


%.tar.gz: %
	tar czf $@ $<

deps/libbz2.tar.gz:
	mkdir -p deps
	wget -O $@ http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
