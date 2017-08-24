# GIRF
An implementation of the GIRF algorithm proposed in the paper written by J. Park and E. L. Ionides (2017).

The program code written in C++ consists of three parts: 
1) the core engine that encodes the filtering algorithm ["stf_parallel.tpp" and "stf_parallel.h" in the main directory], 
2) definition of the POMP model [such as "linmod.c" under directory "BMmodel" or "measlesmod.c" under directory "measlesmodel"], and
3) specifications for which tasks (parameter estimation or likelihood evaluation, etc.) will be conducted with which model and parameter setting [such as "run_linmod_im.c" or "run_measlesmod.c" in each of the directories named by the date the experiment was performed, e.g. 170301, under "BMmodel" or "measlesmodel" directories].

The makefile in each date-named directory automatically builds a program that conducts a specific task. The programs thus built generate result data files, stored in the same directory. These results are analyzed and plotted using the R codes, such as "kalman.R" for the BM model and "analysis.R" for the measles model.

More specific explanations for what the program in each directory does and how plots were generated can be found in the readme file in the directory.
