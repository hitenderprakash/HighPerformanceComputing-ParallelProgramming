#makefile for compiling and running assignment 2: Q1 codes
# hprakash@iu.edu

COMPILER=gcc
CFLAG= -fopenmp

# "make" command will only compile the code
# only one file will compile, @#@ represented commented out lines and will not print to screen
all:
	@# @$(COMPILER) Q1_RaceCondition.c -o Q1_RaceCondition.out $(CFLAG)
	@$(COMPILER) Q1_sol1_removeLoopDependency.c -o Q1_sol1_removeLoopDependency.out $(CFLAG)
	@# @$(COMPILER) Q1_sol2_reduction.c -o Q1_sol2_reduction.out $(CFLAG)
	@# @$(COMPILER) Q1_sol3_critical.c -o Q1_sol3_critical.out $(CFLAG)

# "make run" command will compile ad run the code
# only one file will compile and run @#@ represented commented out lines and will not print to screen
run:
	@#@$(COMPILER) Q1_RaceCondition.c -o Q1_RaceCondition.out $(CFLAG)
	@$(COMPILER) Q1_sol1_removeLoopDependency.c -o Q1_sol1_removeLoopDependency.out $(CFLAG)
	@#@$(COMPILER) Q1_sol2_reduction.c -o Q1_sol2_reduction.out $(CFLAG)
	@#@$(COMPILER) Q1_sol3_critical.c -o Q1_sol3_critical.out $(CFLAG)
	@#@ ./Q1_RaceCondition.out
	@ ./Q1_sol1_removeLoopDependency.out
	@#@ ./Q1_sol2_reduction.out
	@#@ ./Q1_sol3_critical.out

# "make clean" will remove the executable. Redirection used to supress error on console	
clean:
	@rm -rf *.out

