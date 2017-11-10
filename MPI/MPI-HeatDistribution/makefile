#makefile for compiling and running code for parallel quicksort
# hprakash@iu.edu

COMPILER=cc
CFLAG= -lm

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) heat_distribution.c -o heat_distribution.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) heat_distribution.c -o heat_distribution.out $(CFLAG)	
	@ aprun -n 8 ./heat_distribution.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

