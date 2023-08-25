########################################################################
#                                                                      #
#             Makefile for Scripts                                     #
#                                                                      #
#             03-04-2007                                               #
#                                                                      #
########################################################################
#
# Default target:

all: Cgen_null Cinc_field Cscale_field

# Compiler

CC = cc

# Scripts written in C:

Cgen_null: Cgen_null.c
	$(CC) -o Cgen_null Cgen_null.c -lm

Cinc_field: Cinc_field.c
	$(CC) -o Cinc_field Cinc_field.c -lm

Cscale_field: Cscale_field.c
	$(CC) -o Cscale_field Cscale_field.c -lm

# Special targets:

new:
	touch Cgen_null.c Cinc_field.c Cscale_field.c
	make
