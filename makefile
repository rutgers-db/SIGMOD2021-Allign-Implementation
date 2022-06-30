#
# Author: Dong Deng
#

COMPILER = g++
CCFLAGS = -O3 -std=c++11
METHOD = allign allign-window allign-1vn-sentence allign-1vn

all: $(METHOD)

allign: allign.cc segtree.h utils.hpp
	${COMPILER} ${CCFLAGS} -DNORMAL -pg -o $@ $<

clean:
	rm -f $(METHOD) *.o
