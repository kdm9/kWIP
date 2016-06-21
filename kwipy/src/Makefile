CFLAGS ?= -O3
CFLAGS += -fPIC -pthread -I. $(shell pkg-config --cflags hdf5)
LDFLAGS ?=
LDFLAGS += $(shell pkg-config --libs-only-L hdf5)
LIBS += -llz4 $(shell pkg-config --libs-only-l hdf5) -lm -lz -lsz -ldl
BLOSC_CFLAGS = -DHAVE_LZ4 -DHAVE_ZLIB -Iblosc -DSHUFFLE_SSE2_ENABLED
BLOSC_OBJS = blosc/bitshuffle-sse2.o \
			 blosc/bitshuffle-generic.o \
			 blosc/blosc.o \
			 blosc/blosc_filter.o \
			 blosc/blosc_plugin.o \
			 blosc/blosclz.o \
			 blosc/shuffle-sse2.o \
			 blosc/shuffle-generic.o \
			 blosc/shuffle.o
AVX ?= true
ifneq ($(strip $(AVX)),false)
	BLOSC_CFLAGS += -DSHUFFLE_AVX2_ENABLED
	BLOSC_OBJS += blosc/bitshuffle-avx2.o blosc/shuffle-avx2.o
endif

LIB_OBJS = xxhash.o kc_array.o kmerhash.o kmercount.o $(BLOSC_OBJS)

.PHONY: all
all: libkwipy.so kmercount

libkwipy.so: $(LIB_OBJS)
	$(CC) -shared $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

kmercount: main.o $(LIB_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

blosc/bitshuffle-avx2.o: CFLAGS += -mavx2
blosc/shuffle-avx2.o: CFLAGS += -mavx2
blosc/bitshuffle-sse2.o: CFLAGS += -msse2
blosc/shuffle-sse2.o: CFLAGS += -msse2

blosc/%.o: blosc/%.c | $(wildcard blosc/*.h)
	$(CC) $(BLOSC_CFLAGS) $(CFLAGS) -c -o $@ $^


%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f blosc/*.o *.o kmercount libkwipy.so
