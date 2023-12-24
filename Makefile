CC=gcc
CFLAGS=-Wall -Wextra -lglfw -lGL -lm -lOpenCL

all: maxwell

maxwell: maxwell.o
	$(CC) $(CFLAGS) -o maxwell maxwell.o

maxwell.o: maxwell.c 
	$(CC) $(CFLAGS) -c maxwell.c

clean:
	rm -f maxwell maxwell.o

.PHONY: all clean
