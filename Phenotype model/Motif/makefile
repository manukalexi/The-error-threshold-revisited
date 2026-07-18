CC = gcc

CFLAGS = -I/usr/local/include \
          -I/usr/local/include/ViennaRNA \
          -pthread \
          -fno-lto \
          -fopenmp \
          -Ofast

LDFLAGS = -L/usr/local/lib \
           -Wl,-fno-lto

LIBS = -lRNA \
       -lgsl \
       -lgslcblas \
       -lpthread \
       -lstdc++ \
       -lm

motif: motif.c
	$(CC) motif.c $(CFLAGS) $(LDFLAGS) $(LIBS) -o motif

clean:
	rm -f motif