CC = gcc -Wall -Wextra -O0 -g3 -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
#CC = icc -Wall -O3 -g -vec_report0
dep = zcom.h avb.h avp.h GNUmakefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) -o $@ $< -lm

zip: zcom.h
	zip tstatprog *.c zcom.h avb.h avp.h Makefile README license.txt

zcom.h: zcom.full.h GNUmakefile
	./zcompick.py -i zcom.full.h -o zcom.h -k "def util ss rng opt argopt av rc rv2 rv3 md specfunc hist lj rotfit pdb cago" 

clean:
	rm -f *~ *.out MTSEED $(bin) *.pos tstatprog.zip zcom.h
	rm -rf bin obj *.depend *.layout *.dat

.PHONY: clean all

