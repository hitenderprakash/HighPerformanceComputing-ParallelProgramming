#makefile for compiling and running assignment 3: Q2 code
# hprakash@iu.edu

COMPILER=mpicc
CFLAG= 

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) q2_sendrecv.c -o q2_sendrecv.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) q2_sendrecv.c -o q2_sendrecv.out $(CFLAG)	
	@ mpiexec -n 8 ./q2_sendrecv.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

