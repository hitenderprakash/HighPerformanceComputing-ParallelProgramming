#makefile for compiling and running MonteCarlo simulation code
# hprakash@iu.edu

COMPILER=cc
CFLAG= 

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) monteCarlo_MPI.c -o monteCarlo_MPI.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) monteCarlo_MPI.c -o monteCarlo_MPI.out $(CFLAG)	
	@ aprun -n 8 ./monteCarlo_MPI.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

