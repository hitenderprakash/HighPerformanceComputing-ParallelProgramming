#makefile for compiling and running assignment 2 codes
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp

# "make" command will only compile the code
all:
	@$(COMPILER) dotProduct_q1.c -o dotProduct_q1.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) dotProduct_q1.c -o dotProduct_q1.out $(CFLAG)
	@ ./dotProduct_q1.out

# "make aprun"for compiling and running the code on BigRed2 (not tesed, therefor commenting out)
#aprun:
#	@$(COMPILER) dotProduct_q1.c -o dotProduct_q1.out $(CFLAG)
#	@qsub -I -l nodes=2:ppn=16 -q debug_gpu
#	@aprun -n 32 ./dotProduct_q1.out


# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm dotProduct_q1.out >/dev/null || true 
