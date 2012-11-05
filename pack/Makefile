CC = gcc -Wall -Wextra -O0 -g3 -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
# CC = icc -Wall -O3 -g -vec_report0
dep = zcom.h avb.h avp.h Makefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) -o $@ $< -lm

clean:
	rm -f *~ *.out MTSEED $(bin) *.pos zcom.h
	rm -rf bin obj *.depend *.layout *.dat

.PHONY: clean all

