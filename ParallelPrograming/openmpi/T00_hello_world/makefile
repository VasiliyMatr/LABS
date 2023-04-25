CC = mpicc
MPIRUN = mpirun

DEFAULT_FLAGS = -Wall -Wextra -Werror -g
OPT_FLAGS = -O2

SOURCE = main.c
OUT = hello_world

build:
	$(CC) $(DEFAULT_FLAGS) $(OPT_FLAGS) $(SOURCE) -o $(OUT)

run:
	$(MPIRUN) -np 4 ./$(OUT)

