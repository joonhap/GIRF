# GIRF
An implementation of the GIRF algorithm proposed in the paper written by J. Park and E. L. Ionides (2019) (https://arxiv.org/abs/1708.08543).

The program codes written in C++ consist of three types of files: 
1) the core engine that encodes the filtering algorithm ["stf_parallel.tpp" and "stf_parallel.h" in the main directory], 
2) definitions of POMP models [e.g., "linmod.c" under directory "BMmodel" or "measlesmod.c" under directory "measlesmodel"], and
3) executable codes for experiments [e.g., "run_linmod_im.c" in "BMmodel/data/181228" or "run_measlesmod.c" in "measlesmodel/data/180917"].

The makefile in each date-named directory (e.g., "measlesmodel/data/180917") automatically builds a program that can be runned for a specific task. The programs thus built generate data files, stored in the same directory. These results are analyzed and plotted using the R codes, such as "kalman.R" for the BM model and "analysis.R" for the measles model.

More detailed explanations may be found in the readme file in each date-named directory.
