OBJS =$(patsubst %.cpp, %.o, $(wildcard *.cpp))
CC = g++
CFLAGS = -std=c++11 -Wall -pedantic

all: main run

main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o main

run:
	./main

.PHONY: clean

clean:
	rm *.o main

check:
	valgrind ./main
