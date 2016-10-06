bin_PROGRAMS = kwip2

kwip2_CFLAGS = $(PTHREAD_CFLAGS) $(zlib_CFLAGS) $(hdf5_CFLAGS) $(gsl_CFLAGS) $(SIMD_FLAGS) $(CPUEXT_FLAGS)
kwip2_CPPFLAGS = -Iext -Iext/blosc/blosc -Iext/blosc-h5/src
kwip2_LDFLAGS = $(PTHREAD_CFLAGS)
kwip2_LDADD = $(zlib_LIBS) $(hdf5_LIBS) $(gsl_LIBS) $(PTHREAD_LIBS)
kwip2_SOURCES = kc_array.c 								\
				kernelcalc.c 							\
				kmercount.c								\
				kmerhash.c								\
				main.c									\
				\
				ext/blosc.c								\
				ext/blosc-h5/src/blosc_filter.c			\
				ext/blosc-h5/src/blosc_plugin.c			\
				\
				ext/xxhash.c

#
#				\
#				ext/blosc/blosc/bitshuffle-generic.c 	\
#				ext/blosc/blosc/blosc.c					\
#				ext/blosc/blosc/blosclz.c				\
#				ext/blosc/blosc/shuffle.c				\
#				ext/blosc/blosc/shuffle-generic.c		\
#if HAVE_SSE2
#kwip2_SOURCES += ext/blosc/blosc/bitshuffle-sse2.c		\
#				 ext/blosc/blosc/shuffle-sse2.c
#kwip2_CFLAGS += -msse2
#endif
#if HAVE_AVX2
#kwip2_SOURCES += ext/blosc/blosc/bitshuffle-avx2.c		\
#			 	 ext/blosc/blosc/shuffle-avx2.c
#kwip2_CFLAGS += -mavx2
#endif
#
