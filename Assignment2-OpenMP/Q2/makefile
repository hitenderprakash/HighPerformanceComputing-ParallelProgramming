#makefile for compiling and running assignment 2: Q2 codes
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp

# "make" command will only compile the code
all:
	@$(COMPILER) Q2_scheduling.c -o Q2_scheduling.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) Q2_scheduling.c -o Q2_scheduling.out $(CFLAG)
	@ ./Q2_scheduling.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

