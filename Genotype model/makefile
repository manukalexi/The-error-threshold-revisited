CC = gcc

CFLAGS = -I/usr/include
LDFLAGS = -L/usr/lib
LIBS = -lsundials_cvode \
       -lsundials_nvecserial \
       -lsundials_sunlinsoldense \
       -lsundials_sunmatrixdense \
       -lsundials_core \
       -lm

genotype_model: genotype_model.c
	$(CC) genotype_model.c $(CFLAGS) $(LDFLAGS) $(LIBS) -o genotype_model

clean:
	rm -f genotype_model