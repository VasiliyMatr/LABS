CC = mpicc
MPIRUN = mpirun

DEFAULT_FLAGS = -Wall -Wextra -Werror -g
OPT_FLAGS = -O2

SOURCE = main.c
OUT = cycle_send

build:
	$(CC) $(DEFAULT_FLAGS) $(OPT_FLAGS) $(SOURCE) -o $(OUT)

run:
	$(MPIRUN) -np $(JOBS) ./$(OUT)

