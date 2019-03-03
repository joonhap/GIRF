CC = g++
#CCFLAG = -g3 -std=gnu++11
CCFLAG = -Wall -O3 -std=gnu++11 -fopenmp

all : run_measlesmod run_measlesmod_evallik

run_measlesmod : run_measlesmod.o ../../measlesmod_momentmatching.o ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod.o ../../measlesmod_momentmatching.o -o run_measlesmod

run_measlesmod.o : run_measlesmod.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod.c -o run_measlesmod.o

run_measlesmod_evallik : run_measlesmod_evallik.o ../../measlesmod_momentmatching.o ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod_evallik.o ../../measlesmod_momentmatching.o -o run_measlesmod_evallik

run_measlesmod_evallik.o : run_measlesmod_evallik.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod_evallik.c -o run_measlesmod_evallik.o

../../measlesmod_momentmatching.o : ../../measlesmod_momentmatching.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c ../../measlesmod_momentmatching.c -o ../../measlesmod_momentmatching.o

clean :
	rm run_measlesmod.o run_measlesmod run_measlesmod_evallik.o run_measlesmod_evallik
