CFLAGS ?= -O3 -fPIC
CFLAGS += -std=c11 -I.
LDFLAGS ?=
LIBS += -lz

kmercount: main.o xxhash.o kmerhash.o kmercount.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o kmercount
