#	Generic Makefile - Probably requires GnuMake
#	Erik Luijten, original May 1998
#	$Id$

CC	= gcc

CFLAGS	= -Wall -ansi -pedantic

#OPT	= -g
OPT	= -O2 -s

LDFLAGS = -lm

OBJECTS  = alloc2d.o alloc3d.o matrix59_3.o

PROGRAMS = matrix

default: matrix

all: $(PROGRAMS)

%.o : %.c
	$(CC) $(CFLAGS) $(OPT) -c $<

matrix: $(OBJECTS)
	$(CC) -o matrix $(CFLAGS) $(OPT) $(OBJECTS) $(LDFLAGS)

### NO CHANGES BELOW THIS LINE ###

.PHONY: cleanall clean remake

cleanall: clean
	@rm -f $(PROGRAMS)

remake: cleanall default

