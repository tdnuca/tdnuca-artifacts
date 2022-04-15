
CC=mcc
CFLAGS += --ompss --no-copy-deps
LDFLAGS = -L/usr/lib64/atlas 
LIBS = -lm -ltatlas #-lblas -lcblas -llapack

all: cholesky

cholesky: nonnested.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o cholesky nonnested.c $(LIBS)
clean:
	rm -rf cholesky

