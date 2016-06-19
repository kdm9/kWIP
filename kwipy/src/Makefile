CFLAGS ?= -O3
CFLAGS += -std=c11 -I.

kmerhash: main.o xxhash.o kmerhash.o
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean
clean:
	rm -f *.o kmerhash
