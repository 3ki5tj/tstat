CC=icc -O3 -Wall -Wunused

APP=ljmd lju ljmc ljvol gomd

all: $(APP)

$(APP): Makefile zcom.h

ljmd: avb.h

lju: avb.h

ljvol: avp.h

%:: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(APP) *~ a.out

