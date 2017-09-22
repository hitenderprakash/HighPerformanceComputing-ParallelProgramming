#makefile for compiling and running assignment 2: Q2 codes
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp

# "make" command will only compile the code
all:
	@$(COMPILER) Q3_section.c -o Q3_section.out $(CFLAG)

# "make run" command will compile ad run the code
run:
	@$(COMPILER) Q3_section.c -o Q3_section.out $(CFLAG)
	@ ./Q3_section.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

