CC = g++
CFLAGS = -g -Wall

# targets to bring executable up to date
main: main.o sparse_linearsys.o sparse_matrix.o
	$(CC) $(CFLAGS) -o main main.o sparse_linearsys.o sparse_matrix.o

main.o: main.cpp sparse_matrix.h sparse_linearsys.h
	$(CC) $(CFLAGS) -c main.cpp

sparse_linearsys.o: sparse_linearsys.h

sparse_matrix.o: sparse_matrix.h
