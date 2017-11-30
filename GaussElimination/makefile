#makefile for compiling and running gaussElimination matrix inversion program
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp -lm

# "make" command will only compile the code
all:
	@$(COMPILER) gaussElimination.c -o gaussElimination.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) gaussElimination.c -o gaussElimination.out $(CFLAG)
	@ ./gaussElimination.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

