CFLAGS += -I . -march=native -mtune=native -O3

all: sparkle-test.exe

sparkle-test.exe: sparkle_vectors.h sparkle-test.c 
	$(CC) $(CFLAGS) $^ -o $@
	./$@

.PHONY : clean
clean: 
	rm -rf *.o *.exe

include $(patsubst %.c,%.d,$(SOURCES))
	