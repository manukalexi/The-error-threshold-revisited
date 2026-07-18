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

large_loop: large_loop.c
	$(CC) large_loop.c $(CFLAGS) $(LDFLAGS) $(LIBS) -o large_loop

clean:
	rm -f large_loop