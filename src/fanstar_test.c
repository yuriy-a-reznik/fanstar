/*!
 *  \file       fanstar_test.c
 *  \brief      Test program for closest-point functions for An* and few other lattices.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 */

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fanstar.h"

/*! Local constants */
#define NUM_SAMPLES  4000000      /*!< the number of random samples to test                   */
#define SEED         12345        /*<! random seed for repeatable results                     */

/************
 * 
 *  Helper functions:
 * 
 *   - generate random sample array
 *   - l2 distance
 *   - lift/reduce to An* ambient space
 */

/*!
 *  \brief Generate an array if input samples.
 *
 *  This function fills the provided array with random generated samples.
 *  Uniform distribution in [-r..r] range across each dimension is used.
 */
static void gen_samples(double x[][LQ_MAX_M], int samples, int n, int r)
{
    int i, j;
    double u;

    /* reset rand: */
    srand((unsigned)SEED);

    /* scan samples: */
    for (i = 0; i < samples; i++) {
        /* generate random sample: */
        for (j = 0; j < n; j++) {
            /* pull random number in [0..1] */
            u = ((double)rand() / (double)RAND_MAX);
            /* map it to [-r..r] range: */
            u = (u - 0.5) * (double)(2 * r);
            x[i][j] = u;
        }
    }
}

/*!
 *  \brief L2 distance:
 */
static double l2dist(double *x, double *y, int n)
{
    double err = 0., d;
    int i;
    for (i = 0; i < n; i++) {
        d = x[i] - y[i];
        err += d * d;
    }
    return sqrt(err);
}

/*!
 *  \brief  Isometric lift R^n -> zero-sum hyperplane in R^{n+1} (A_n ambient space).
 *
 *   Maps x (n coords) to y (n+1 coords) with sum(y)=0, using the orthonormal Helmert basis. 
 *   Distances are preserved exactly:  ||y - y'|| == ||x - x'||.
 *   Runs in O(n) via a running suffix sum.
 *
 *  \param[in]  x  n-vector.
 *  \param[out] y  (n+1)-vector on the hyperplane (caller-allocated).
 *  \param[in]  n  source dimension (n >= 1).
 */
static void lift(double* x, double* y, int n)
{
    double S = 0.0;             /* S_j = sum_{r=j}^{n} x[r-1]/sqrt(r(r+1)) */
    int    j;
    for (j = n + 1; j >= 1; --j) {
        double yj;
        if (j <= n) S += x[j - 1] / sqrt((double)j * (double)(j + 1));
        yj = S;
        if (j >= 2) yj -= x[j - 2] * sqrt((double)(j - 1) / (double)j);
        y[j - 1] = yj;
    }
}

/*!
 *  \brief Inverse map: zero-sum (n+1)-vector -> R^n. Exact inverse of lift().
 *
 *   If y is not exactly zero-sum, this returns the coordinates of y's orthogonal
 *   projection onto the hyperplane (the all-ones component is discarded), so
 *   lift(reduce(y)) == projection of y onto H.  Runs in O(n).
 *
 *  \param[in]  y  (n+1)-vector.
 *  \param[out] x  n-vector (caller-allocated).
 *  \param[in]  n  target dimension (n >= 1).
 */
static void reduce(double* y, double* x, int n)
{
    double P = 0.0;             /* P_r = sum_{j=1}^{r} y[j-1] */
    int    r;
    for (r = 1; r <= n; ++r) {
        P += y[r - 1];
        x[r - 1] = (P - (double)r * y[r]) / sqrt((double)r * (double)(r + 1));
    }
}

/*!
 *  \brief Lift all samles in a given array.
 */
static void lift_samples(double x[][LQ_MAX_M], double y[][LQ_MAX_M], int samples, int n)
{
    int i;
    for (i = 0; i < samples; i++)
        lift(x[i], y[i], n);
}

/*!
 *  \brief Reduce all samles in a given array.
 */
static void reduce_samples(double y[][LQ_MAX_M], double x[][LQ_MAX_M], int samples, int n)
{
    int i;
    for (i = 0; i < samples; i++)
        reduce(y[i], x[i], n);
}

/******************
 * 
 *  Data structures describing parameters/conditions of each test:
 * 
 */

/*! Parameters for each test: */
typedef struct {
	const char* lattice;                                /* lattice name */
	int n_min, n_max;                                   /* lattice dimension range */
    int lifted;                                         /* flag: 1 if we need to lift x[] to ambient space, 0 otherwise */
    const char* cp_method;                              /* closest-point method name */
    int (*cp_func)(const double*, double*, int);        /* closest point function pointer */
    int (*ref_cp_func)(const double*, double*, int);    /* reference closest point function pointer (NULL if N/A) */
} test_params_t;

/* list of test cases: */
static test_params_t tests[] = {
    /* tests of different implementations of closest-point functions for An* lattices: */
    { "An*", 2, 8, 1, "CS", closest_point_Anstar_CS, NULL},
    { "An*", 2, 8, 1, "VC", closest_point_Anstar_VC, closest_point_Anstar_CS },
    { "An*", 2, 8, 1, "MCQ", closest_point_Anstar_MCQ, closest_point_Anstar_CS },
    { "An*", 2, 8, 1, "MCSQ", closest_point_Anstar_MCSQ, closest_point_Anstar_CS },
    { "An*", 2, 8, 1, "MCQ_opt", closest_point_Anstar_MCQ_opt, closest_point_Anstar_CS },
    { "An*", 2, 8, 1, "MCSQ_opt", closest_point_Anstar_MCSQ_opt, closest_point_Anstar_CS },
    { "An*", 2, 8, 1, "Proposed", closest_point_Anstar_faster, closest_point_Anstar_CS },
    /* also test some alternatice lattices - for comparison: */
    { "Zn",  2, 8, 0, "CS", closest_point_Zn_CS, NULL },
    { "An",  2, 8, 1, "CS", closest_point_An_CS, NULL },
    { "Dn",  2, 8, 0, "CS", closest_point_Dn_CS, NULL},
    { "Dn*", 2, 8, 0, "CS", closest_point_Dnstar_CS, NULL},
    { "En",  8, 8, 0, "CS", closest_point_E8_CS, NULL}
};

/* number of tests: */
static int num_tests = sizeof(tests) / sizeof(tests[0]);

/*!
 *  \brief The main test program.
 */
int main(void)
{
    /* buffers: */
    static double x[NUM_SAMPLES][LQ_MAX_M];         /* random generated samples (n-dimensional values) */
    static double x_amb[NUM_SAMPLES][LQ_MAX_M];     /* samples converted to lattice ambient space */
    static double y[NUM_SAMPLES][LQ_MAX_M];         /* nearest points in n-dimensional space */
    static double y_amb[NUM_SAMPLES][LQ_MAX_M];     /* nearest points in ambient space */
    static double y_ref[NUM_SAMPLES][LQ_MAX_M];     /* nearest points found by reference algorithm */
    static double y_amb_ref[NUM_SAMPLES][LQ_MAX_M]; /* nearest points in ambient space found by reference algorithm */

    /* local variables: */
    clock_t start, end;
    double exec_time;                               /* execution time */
    double ave_err, max_err, ref_err, err;          /* average and wort case distances between x and y and y and y_ref */
	int n_mismatches;                               /* number of mismatches between y and y_ref */
    int r = 10;                                   /* range of per-coordinate values in x[] */
    int i, j, n, rc;

    /* scan test cases: */
    for (i = 0; i< num_tests; i++) {

        /* print test descptiption: */
        printf("\nLattice: %s, closest-point algorithm: %s\n\n", tests[i].lattice, tests[i].cp_method);

        /* print header for table with results: */
        printf("  %3s   %11s %11s %11s %11s\n",
            "n",                /* n                                               */
            "Avg L2 err",       /* Average L2 error                                */
            "Max L2 err",       /* Maximum L2 error                                */
            "Ex. Time",         /* Execution time [sec / 1M iterations]            */
            "Mismatches");      /* # of mismatches relative to reference algorithm */

        /* scan dimensions: */
        for (n = tests[i].n_min; n <= tests[i].n_max; n++) {

            /* zero everything: */
            memset(x, 0, sizeof(x));
            memset(x_amb, 0, sizeof(x_amb));
            memset(y, 0, sizeof(y));
            memset(y_amb, 0, sizeof(y_amb));
            memset(y_ref, 0, sizeof(y_ref));
            memset(y_amb_ref, 0, sizeof(y_amb_ref));

            /* generate an array if input samples: */
            gen_samples(x, NUM_SAMPLES, n, r);
 
            /* run algorithm under the test: */
            if (tests[i].lifted) {
                /* go to ambient space: */
                lift_samples(x, x_amb, NUM_SAMPLES, n);
                start = clock();
                /* run tests using lifted vectors: */
                for (j = 0; j < NUM_SAMPLES; j++) {
                    rc = tests[i].cp_func(x_amb[j], y_amb[j], n);
                    if (rc != LQ_SUCCESS) { fprintf(stderr, "\n => Closest-point error %d for n=%d, j=%d\n", rc, n, j); return 1; }
                }
                end = clock();
                /* reduce results back to R^n space: */
                reduce_samples(y_amb, y, NUM_SAMPLES, n);
            } else {
                /* run tests using original vectors: */
                start = clock();
                for (j = 0; j < NUM_SAMPLES; j++) {
                    rc = tests[i].cp_func(x[j], y[j], n);
                    if (rc != LQ_SUCCESS) { fprintf(stderr, "\n => Closest-point error %d for n=%d, j=%d\n", rc, n, j); return 1; }
                }
                end = clock();
            }

            /* compute execution time: */
            exec_time = (double)(end - start) / (double)CLOCKS_PER_SEC;

            /* compute average and worst-case errors for algorithm under test: */
            ave_err = 0.0; max_err = 0.0;
            for (j = 0; j < NUM_SAMPLES; j++) {
                err = l2dist(x[j], y[j], n);
                ave_err += err;
                if (err > max_err) max_err = err;
            }
            ave_err /= (double)NUM_SAMPLES;

            /* check if reference algorithm is provided: */
            n_mismatches = 0;
            if (tests[i].ref_cp_func != NULL) {

                /* run reference algorithm: */
                if (tests[i].lifted) {
                    lift_samples(x, x_amb, NUM_SAMPLES, n);
                    for (j = 0; j < NUM_SAMPLES; j++) {
                        rc = tests[i].ref_cp_func(x_amb[j], y_amb_ref[j], n);
                        if (rc != LQ_SUCCESS) { fprintf(stderr, "\n => Reference closest-point error %d for n=%d, j=%d\n", rc, n, j); return 1; }
                    }
                    reduce_samples(y_amb_ref, y_ref, NUM_SAMPLES, n);
                }
                else {
                    for (j = 0; j < NUM_SAMPLES; j++) {
                        rc = tests[i].ref_cp_func(x[j], y_ref[j], n);
                        if (rc != LQ_SUCCESS) { fprintf(stderr, "\n => Reference closest-point error %d for n=%d, j=%d\n", rc, n, j); return 1; }
                    }
                }

                /* check consistency with reference algorithm: */
                for (j = 0; j < NUM_SAMPLES; j++) {

                    /* compute distances between x and reconstructions: */
                    err = l2dist(x[j], y[j], n);
                    ref_err = l2dist(x[j], y_ref[j], n);
                    if (fabs(err - ref_err) > 1e-12) {
                        /* this is a problem. these distances must be the same! */
                        fprintf(stderr, "\n => Algorithm under test does not find a point at the same distance as reference algorithm."
                                        "\n => n=%d, j=%d, L2dist=%g, RefL2dist=%g\n", n, j, err, ref_err); return 1;
                    }

                    /* while being equidistant from x, the points returned by the reference algorithm and the algorithm under the test may not match. 
                     * these are not exactly errors, but we will still count such cases to see which algorithm is more consistent 
                     * to the reference: */
                    err = l2dist(y[j], y_ref[j], n);
                    if (err > 1e-12) 
                        n_mismatches++;
                }
            }

            /* print the results: */
            printf("  %3d   %11.4f %11.4f %11.4f %11d\n",
                n,                  /* n                                               */
                ave_err,            /* Average L2 error                                */
                max_err,            /* Maximum L2 error                                */
                exec_time,          /* Execution time [sec / 1M iterations]            */
                n_mismatches);      /* % of mismatches relative to reference algorithm */
        }
    }

    return 0;
}
 
/* fanstar_test.c - end of file */
