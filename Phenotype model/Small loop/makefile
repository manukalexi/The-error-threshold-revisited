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

small_loop: small_loop.c
	$(CC) small_loop.c $(CFLAGS) $(LDFLAGS) $(LIBS) -o small_loop

clean:
	rm -f small_loop