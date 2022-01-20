#include <stdio.h>
#include "stdlib.h"
#include "Math.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>


void HinfinityInC(double gama,double dt,double duration,double Steadystate,
					gsl_matrix *pos,
					gsl_matrix *vel,
					gsl_matrix *poshat,
					gsl_matrix *velhat,
					gsl_matrix *poshatinf,
					gsl_matrix *velhatinf,
					gsl_matrix *hinfgains,
					gsl_matrix *kalmangains)

gsl_matrix * invert_a_matrix(gsl_matrix *matrix, size_t size)

void print_mat_contents(gsl_matrix *matrix, size_t size)

double rand_generator()

