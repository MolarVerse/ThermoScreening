CC = gcc
CFLAGS = -Wall -Wextra -std=gnu11
LIBS = -lm

.PHONY: all
all: gibbs_dftb

.PHONY: clean
clean:
	$(RM) gibbs_dftb

gibbs_dftb: gibbs_dftb.c external/symmetry_i.c atom_properties.c
	$(CC) $(CFLAGS) $^ $(LIBS) -o $@
