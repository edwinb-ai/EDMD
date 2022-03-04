# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdSingle.c mdSingle.h xoshiro256plus.c splitmix64.c
	$(CC) $(CFLAGS) mdSingle.c -o md  $(LIBS)
