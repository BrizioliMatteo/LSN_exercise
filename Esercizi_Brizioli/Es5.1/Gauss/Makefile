CC = g++
CFLAGS = -Wall -O3 --std=c++11

Hydrogen_atom_gauss.exe : Hydrogen_atom_gauss.o random.o
	$(CC) random.o Hydrogen_atom_gauss.o -o Hydrogen_atom_gauss.exe
Hydrogen_atom_gauss.o : Hydrogen_atom_gauss.cpp Hydrogen_atom_gauss.h
	$(CC) -c Hydrogen_atom_gauss.cpp -o Hydrogen_atom_gauss.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Hydrogen_atom.exe
