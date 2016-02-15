CC = g++
CFLAGS = -c -std=c++11 -O2 -lm -Wall

all: macroscopic

macroscopic: main.o monodomain.o purkinje.o queue.o
	$(CC) -o macroscopic main.o monodomain.o purkinje.o queue.o

monodomain.o: monodomain.cpp
	$(CC) $(CFLAGS) monodomain.cpp

purkinje.o: purkinje.cpp
	$(CC) $(CFLAGS) purkinje.cpp

queue.o: queue.cpp
	$(CC) $(CFLAGS) queue.cpp

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

clean:
	rm -r VTK
	rm -rf *o macroscopic
	rm graph.plt
	rm data.dat
	rm monodomain.png
