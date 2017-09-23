CC = g++
#CCFLAG = -g3 -std=gnu++11
CCFLAG = -Wall -O3 -std=gnu++11 -fopenmp

all : run_linmod_im

run_linmod_im : run_linmod_im.o ../../linmod.o ../../linmod.h ../../stf_parallel.tpp ../../stf_parallel.h
	$(CC) $(CCFLAG) run_linmod_im.o ../../linmod.o -o run_linmod_im

run_linmod_im.o : run_linmod_im.c ../../linmod.h ../../stf_parallel.tpp ../../stf_parallel.h
	$(CC) $(CCFLAG) -c run_linmod_im.c

../../linmod.o : ../../linmod.c ../../linmod.h ../../stf_parallel.tpp ../../stf_parallel.h
	$(CC) $(CCFLAG) -c ../../linmod.c -o ../../linmod.o

clean :
	rm ../../linmod.o run_linmod_im.o run_linmod_im
