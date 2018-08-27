SOURCE = a_a.c  \
	cascade.c  \
	do_simulation.c  \
	manage_transitions.c  \
	read_files.c  \
	b_b.c  \
	honl-london.c    \
	manage_rotations.c  \
	rvesim.c 
CFLAGS = -Wall -O0 -lm
NAME = -o rvesim_new

default: $(SOURCE) 
	gcc $(SOURCE)  $(NAME) $(CFLAGS)

gdb: $(SOURCE) 
	gcc  $(SOURCE) $(NAME) $(CFLAGS) -g

