CC = g++
#CCFLAG = -g3 -std=gnu++11
CCFLAG = -Wall -O3 -std=gnu++11 -fopenmp

all : run_measlesmod_realdata run_measlesmod_evallik run_measlesmod_simulate

run_measlesmod_simulate : run_measlesmod_simulate.o ../../measlesmod_varyingrep.o ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod_simulate.o ../../measlesmod_varyingrep.o -o run_measlesmod_simulate

run_measlesmod_simulate.o : run_measlesmod_simulate.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod_simulate.c -o run_measlesmod_simulate.o

run_measlesmod_evallik : run_measlesmod_evallik.o ../../measlesmod_varyingrep.o ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod_evallik.o ../../measlesmod_varyingrep.o -o run_measlesmod_evallik

run_measlesmod_evallik.o : run_measlesmod_evallik.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod_evallik.c -o run_measlesmod_evallik.o

run_measlesmod_realdata : run_measlesmod_realdata.o ../../measlesmod_varyingrep.o ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) run_measlesmod_realdata.o ../../measlesmod_varyingrep.o -o run_measlesmod_realdata

run_measlesmod_realdata.o : run_measlesmod_realdata.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c run_measlesmod_realdata.c -o run_measlesmod_realdata.o

../../measlesmod_varyingrep.o : ../../measlesmod_varyingrep.c ../../../stf_parallel.h ../../../stf_parallel.tpp ../../measlesmod.h
	$(CC) $(CCFLAG) -c ../../measlesmod_varyingrep.c -o ../../measlesmod_varyingrep.o

clean :
	rm run_measlesmod_realdata.o run_measlesmod_realdata run_measlesmod_evallik.o run_measlesmod_evallik
