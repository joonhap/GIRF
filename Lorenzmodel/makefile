CC = g++
#CCFLAG = -g3 -std=gnu++11
CCFLAG = -Wall -O3 -std=gnu++11 -fopenmp

all : run_Lorenzmod

run_Lorenzmod : run_Lorenzmod.o ../Lorenzmod.o ../../stf_parallel.h ../../stf_parallel.tpp ../Lorenzmod.h
	$(CC) $(CCFLAG) run_Lorenzmod.o ../Lorenzmod.o -o run_Lorenzmod

run_Lorenzmod.o : run_Lorenzmod.c ../../stf_parallel.h ../../stf_parallel.tpp ../Lorenzmod.h
	$(CC) $(CCFLAG) -c run_Lorenzmod.c -o run_Lorenzmod.o

../Lorenzmod.o : ../Lorenzmod.c ../../stf_parallel.h ../../stf_parallel.tpp ../Lorenzmod.h
	$(CC) $(CCFLAG) -c ../Lorenzmod.c -o ../Lorenzmod.o

clean :
	rm run_Lorenzmod.o run_Lorenzmod ../Lorenzmod.o
