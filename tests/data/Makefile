
all: defined-1.ct defined-2.ct defined-3.ct defined-4.ct

defined-%.ct: defined-%.fa
	rm -f $@*
	load-into-countgraph.py -f -x 100 -N 1 -k 5 $@ $<
	rm -f $@.info

