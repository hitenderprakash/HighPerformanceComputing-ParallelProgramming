#makefile for compiling and running assignment 3: Q4 code
# hprakash@iu.edu

COMPILER=mpicc
CFLAG= 

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@$(COMPILER) q4_mandelbrot.c -o q4_mandelbrot.out $(CFLAG)


# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@$(COMPILER) q4_mandelbrot.c -o q4_mandelbrot.out $(CFLAG)	
	@ mpiexec -n 8 ./q4_mandelbrot.out 256 256

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

