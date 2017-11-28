#makefile for compiling and running code for parallel I/O
# hprakash@iu.edu

COMPILER=mpicc
CFLAG= 

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) mpi_fs.c -o mpi_fs.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) mpi_fs.c -o mpi_fs.out $(CFLAG)	
	@ mpirun -n 8 ./mpi_fs.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

