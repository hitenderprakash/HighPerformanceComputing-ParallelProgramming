#makefile for compiling and running assignment 3: Q3 code
# hprakash@iu.edu

COMPILER=mpicc
CFLAG= -lm

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) q3_dotprod.c -o q3_dotprod.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) q3_dotprod.c -o q3_dotprod.out $(CFLAG)	
	@ mpiexec -n 8 ./q3_dotprod.out 

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

