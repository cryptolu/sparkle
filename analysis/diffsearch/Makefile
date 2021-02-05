CC = g++
CFLAGS = -std=c++11 -O3
DEPS = defs.h xdp_add.h trail.h sparkle_best_trail_search.h
OBJ = trail.o sparkle_best_trail_search.o main_serial.o
TARGET = diffsearch

all: compile

compile: $(OBJ)
	$(CC) -o $(TARGET) $^ $(CFLAGS)

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	/bin/rm -f *.o $(TARGET) 
