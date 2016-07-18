CeFLAGS ?= -O3
CFLAGS += -fPIC -std=c99 -pthread -I. -Iext $(shell pkg-config --cflags hdf5)
LDFLAGS ?=
LDFLAGS += $(shell pkg-config --libs-only-L hdf5)
LIBS += -llz4 $(shell pkg-config --libs-only-l hdf5) -lm -lz -lsz -ldl -lpthread
BLOSC_CFLAGS = -DHAVE_LZ4 -DHAVE_ZLIB -Iext/blosc -DSHUFFLE_SSE2_ENABLED
BLOSC_OBJS = ext/blosc/bitshuffle-sse2.o \
			 ext/blosc/shuffle-sse2.o \
			 ext/blosc/bitshuffle-generic.o \
			 ext/blosc/blosc.o \
			 ext/blosc/blosc_filter.o \
			 ext/blosc/blosc_plugin.o \
			 ext/blosc/blosclz.o \
			 ext/blosc/shuffle-generic.o \
			 ext/blosc/shuffle.o
AVX ?= true
ifneq ($(strip $(AVX)),false)
	BLOSC_CFLAGS += -DSHUFFLE_AVX2_ENABLED
	BLOSC_OBJS += ext/blosc/bitshuffle-avx2.o ext/blosc/shuffle-avx2.o
endif

LIB_OBJS = ext/xxhash.o kc_array.o kmerhash.o kmercount.o $(BLOSC_OBJS)

.PHONY: all
all: libkwipy.so kmercount test.elf

.PHONY: test
test: test.elf
	./test.elf

libkwipy.so: $(LIB_OBJS)
	$(CC) -shared $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

kmercount: main.o $(LIB_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

test.o: test.c $(wildcard tests/*.c)
	$(CC) -Wall $(CFLAGS) -c -o $@ $<

test.elf: test.o $(LIB_OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

ext/blosc/bitshuffle-avx2.o: CFLAGS += -mavx2
ext/blosc/shuffle-avx2.o: CFLAGS += -mavx2
ext/blosc/bitshuffle-sse2.o: CFLAGS += -msse2
ext/blosc/shuffle-sse2.o: CFLAGS += -msse2

ext/blosc/%.o: ext/blosc/%.c | $(wildcard ext/blosc/*.h)
	$(CC) $(BLOSC_CFLAGS) $(CFLAGS) -c -o $@ $^


%.o: %.c %.h
	$(CC) -Wall $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f ext/blosc/*.o ext/*.o *.o kmercount libkwipy.so test.elf
