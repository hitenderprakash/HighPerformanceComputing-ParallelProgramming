#makefile for compiling and running integration of trignometric function code
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -lm

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) integrationCalc.c -o integrationCalc.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) integrationCalc.c -o integrationCalc.out $(CFLAG)	
	@ ./integrationCalc.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

