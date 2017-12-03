#makefile for compiling and running gaussElimination matrix inversion program
# hprakash@iu.edu

#compiling for GNU
#change the compiler to cc for Cray
#Make sure that before running the program, appropriate module is loaded
#to run on BigRed2 make sure to change the default cray module to gnu 
# and run like ....
# step1: module swap PrgEnv-cary PrgEnv-gnu
# step2: make (i.e compile for gcc)
# step3: export OMP_NUM_THREADS=m ; m is number of threads
# step4: time aprun -d m ./gaussElimination.out <rows> <cols>
#  eg  : time aprun -d 8 ./gaussElimination.out 1000 1000 ; threads=8, row=1000, cols=1000 do not forget to set OMP_NUM_THREADS=8

COMPILER=gcc
CFLAG= -fopenmp -lm

# "make" command will only compile the code
all:
	@$(COMPILER) gaussElimination.c -o gaussElimination.out $(CFLAG)

# "make run" command will compile ad run the code
#  recommended to run manually as suggested in comments given above, make run will run it for 8 threads (provided  OMP_NUM_THREADS=8 is set)
run:
	@$(COMPILER) gaussElimination.c -o gaussElimination.out $(CFLAG)
	@ aprun -d 8 ./gaussElimination.out 1000 1000

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

#also delete result file (not required,when you run the program again it will anyway overwrite the file)	
kleen:
	@rm -rf *.out
	@rm results.txt

