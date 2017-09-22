#makefile for compiling and running assignment 2: Q4 codes
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp -lm

# "make" command will only compile the code
all:
	@$(COMPILER) Q4_MatrixVecMultiply_strongScaling.c -o Q4_MatrixVecMultiply_strongScaling.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) Q4_MatrixVecMultiply_strongScaling.c -o Q4_MatrixVecMultiply_strongScaling.out $(CFLAG)
	@ ./Q4_MatrixVecMultiply_strongScaling.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

