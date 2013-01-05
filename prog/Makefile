CC=icc -O3 -Wall -Wunused

APP=ljmd lju ljmc ljvol gomd enerw volrw

all: $(APP)

$(APP): Makefile zcom.h

ljmd: avb.h

lju: avb.h

gomd: avb.h

ljmc: avb.h

ljvol: avp.h

enerw: avb.h

volrw: avp.h

%:: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(APP) *~ a.out  *avb.dat avp.dat epot*.his vol*.his *.pos

