#makefile for compiling and running MatrixVector_mul_l2norm.c
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -lm

# "make" command will only compile the code
all:
	@$(COMPILER) MatrixVector_mul_l2norm.c -o MatrixVector_mul_l2norm.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) MatrixVector_mul_l2norm.c -o MatrixVector_mul_l2norm.out $(CFLAG)
	@ ./MatrixVector_mul_l2norm.out

# "make aprun"for compiling and running the code on BigRed2 (not tesed, therefor commenting out)
#aprun:
#	@$(COMPILER) MatrixVector_mul_l2norm.c -o MatrixVector_mul_l2norm.out $(CFLAG)
#	@qsub -I -l nodes=2:ppn=16 -q debug_gpu
#	@aprun -n 32 ./MatrixVector_mul_l2norm.out


# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm MatrixVector_mul_l2norm.out >/dev/null || true 
