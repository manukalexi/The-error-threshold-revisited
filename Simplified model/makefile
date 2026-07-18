CC = gcc

CFLAGS = -I/usr/include

LDFLAGS = -L/usr/lib

LIBS = -lsundials_cvode \
       -lsundials_nvecserial \
       -lsundials_sunlinsoldense \
       -lsundials_sunmatrixdense \
       -lsundials_core \
       -lm

simplified_model: simplified_model.c
	$(CC) simplified_model.c $(CFLAGS) $(LDFLAGS) $(LIBS) -o simplified_model

clean:
	rm -f simplified_model