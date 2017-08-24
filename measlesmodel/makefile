CC = g++
#CCFLAG = -g3 -std=gnu++11
CCFLAG = -Wall -O3 -std=gnu++11 -fopenmp

all : run_measlesmod

run_measlesmod : run_measlesmod.o measlesmod.o ../stf_parallel.h ../stf_parallel.tpp measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod.o measlesmod.o -o run_measlesmod

run_measlesmod.o : run_measlesmod.c ../stf_parallel.h ../stf_parallel.tpp measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod.c

measlesmod.o : measlesmod.c ../stf_parallel.h ../stf_parallel.tpp measlesmod.h
	$(CC) $(CCFLAG) -c measlesmod.c

clean :
	rm measlesmod.o run_measlesmod.o run_measlesmod

