CC= bspcc
CFLAGS= -std=c99 -Wall -O3 -I/home/rick/Documents/BSPedupack2.0
LFLAGS= -lm
EDUPACKPATH = /home/rick/Documents/BSPedupack2.0

OBJS = bspmattoicrs.o $(EDUPACKPATH)/bspsparse_input.o $(EDUPACKPATH)/bspedupack.o 

mattoicrs: $(OBJS)
	 $(CC) $(CFLAGS) -o mattoicrs $(OBJS) $(LFLAGS)


bspmattoicrs.o: $(EDUPACKPATH)/bspedupack.h $(EDUPACKPATH)/bspsparse_input.h
$(EDUPACKPATH)/bspsparse_input.o: $(EDUPACKPATH)/bspedupack.h $(EDUPACKPATH)/bspsparse_input.h
$(EDUPACKPATH)/bspedupack.o: $(EDUPACKPATH)/bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJS) mattoicrs
