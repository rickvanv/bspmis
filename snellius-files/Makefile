CC= bspcc
EDUPACKPATH = /home/bissstud/Students22/RickLisa/BSPedupack2.0
CFLAGS= -std=c99 -Wall -O3 -I/home/bissstud/BSPonMPI/include -I$(EDUPACKPATH)
LFLAGS= -lm


OBJS = bspmis.o $(EDUPACKPATH)/bspsparse_input.o $(EDUPACKPATH)/bspedupack.o 

bspmis: $(OBJS)
	 $(CC) $(CFLAGS) -o bspmis $(OBJS) $(LFLAGS)


bspmis.o: $(EDUPACKPATH)/bspedupack.h $(EDUPACKPATH)/bspsparse_input.h
$(EDUPACKPATH)/bspsparse_input.o: $(EDUPACKPATH)/bspedupack.h $(EDUPACKPATH)/bspsparse_input.h
$(EDUPACKPATH)/bspedupack.o: $(EDUPACKPATH)/bspedupack.h

.PHONY: clean
clean:
	rm -f $(OBJS) bspmis
