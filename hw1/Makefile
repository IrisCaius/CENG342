CC=mpicc
CFLAGS=-I.
DEPS = 19050111032.h
OBJ = 19050111032.o hellofunc.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

19050111032: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm *.o