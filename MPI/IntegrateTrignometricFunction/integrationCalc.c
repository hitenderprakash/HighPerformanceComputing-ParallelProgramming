/* Program to integrate cos(x)*sin(x/2) from 0 to Pi/2
 * Hitender Prakash (hprakash@iu.edu)
 * Version 1: Sequential code
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //link with -lm compiler flag

#define PI 3.14159265


int main(){
	double c=0.0;
	double incr=0.0001;
	int reps=(PI/2)/incr;
	printf("\nReps will be : %f", reps);
	int count=0;
	double integral=0.0;
	while (c<=(PI/2)){
		integral+= cos(c)*sin(c/2);
		c+=incr;
		count++;
	}
	printf("\nThe val: %lf and done in %d times", integral*incr,count);
	printf("\nFinal val of c: %lf",c);
	return 0;
}
