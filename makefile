COMPILER=gcc
CFLAG= -lm
all:
	@$(COMPILER) MatrixVector_mul_l2norm.c -o MatrixVector_mul_l2norm $(CFLAG)

clean:
	@rm MatrixVector_mul_l2norm


