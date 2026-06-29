/*!
 *  \file       fanstar_MCSQ.c
 *  \brief      Closest-point algorithms for An* lattices -- using McKilliam, Clarkson, Smith, & Quinn method.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements closest-point finding algorithm for An* lattices by following McKilliam, Clarkson, Smith, & Quinn paper:
 *
 *   R. G. McKilliam, I. V. L. Clarkson, W. D. Smith, and B. G. Quinn, "A linear-time nearest point algorithm for the lattice A*_n",
 *   Int. Symp. Information Theory and its Applications (ISITA), Auckland, New Zealand, Dec. 2008.
 *
 *  The lattice A*_n is the projection of the cubic lattice Z^(n+1) onto the hyperplane orthogonal to 1:
 *
 *     A*_n = { Q k : k in Z^(n+1) },   Q = I - (1/(n+1)) 1 1^T.
 *
 *  Like the O(n log n) MCQ algorithm, this method tests the n+1 relevant candidates
 *
 *     u_i = f(i/(n+1)) = round_cs((x) + sum_{j in K_i} e_j,   i = 0..n,
 *
 *  (Lemma 3), and keeps the one minimising d_i = ||Q(x - u_i)||^2.  
 * 
 *  The novelty is that the index sets needed to build the candidates do NOT require a full
 *  comparison sort: each coordinate t is dropped into one of n+1 buckets by
 *
 *     bucket(t) = (n+1) - floor( (n+1) * (z_t + 1/2) ),   z = x - round_cs((x),
 *
 *  which is exactly an O(n) bucket sort (Algorithm 2 of the paper, implemented here with the bucket/link array pair).  
 *  Processing buckets in order i = 1..n+1 visits the coordinates in order of decreasing fractional part, so the recurrences
 *
 *     alpha_i = alpha_{i-1} - |S_i|,
 *     beta_i  = beta_{i-1} + |S_i| - 2 * sum_{t in S_i} z_t,
 *
 *  update d_i in O(1) per coordinate. The total cost is O(n) -- no sort.
 */

#include <math.h>
#include "fanstar.h"

/*!
 *  \brief Round x to the nearest integer, with ties broken toward zero (Conway & Sloan's f() function).
 *
 *  Examples: f(0.5) = 0, f(-0.5) = 0, f(1.5) = 1, f(-1.5) = -1.
 *
 *  \param[in] x  Real input.
 *  \return       The closest integer to x, returned as a double.
 */
static double round_cs(double x)
{
    double a = fabs(x), r = floor(a + 0.5);
    if (r - a >= 0.5) r -= 1.0;             /* exact tie -> toward zero */
    return (x >= 0.0) ? r : -r;
}

/*!
 *  \brief Closest point in An* finding by using MCSQ algorithm (McKilliam, Clarkson, Smith, & Quinn, 2008).
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest An* point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Anstar_MCSQ(const double* x, double* y, int n)
{
    double z[LQ_MAX_M];         /* centred fractional parts z = x - round_cs((x)        */
    double k[LQ_MAX_M];         /* k = round_cs((x), later incremented to u_m           */
    int    bucket[LQ_MAX_M];    /* bucket[b] = head index of bucket b, or -1        */
    int    link[LQ_MAX_M];      /* link[t]   = next index in t's bucket, or -1      */

    double alpha, beta;         /* z_i^T 1 and z_i^T z_i, updated recursively       */
    double d, D;                /* current candidate distance and running minimum   */
    double mu;                  /* mean, used by the projection Q k = k - mean(k)   */
    int    dim = n + 1;         /* ambient dimension n + 1                          */
    int    i, t, m, bi;

    /* check arguments: */
    if (x == NULL || y == NULL || n < 1) 
        return -1;

    /*
     * Compute: 
     *  z = x - round_cs((x), and k = round_cs((x).
     *  The residual lies in [-1/2, 1/2).
     */
    alpha = 0.0;
    beta = 0.0;
    for (i = 0; i < dim; i++) {
        k[i] = floor(x[i] + 0.5);
        z[i] = x[i] - k[i];
        alpha += z[i];
        beta += z[i] * z[i];
    }

    /*
     * Bucket sort (Algorithm 2, lines 4-8).
     *
     *  Bucket of coordinate t is
     *     b = dim - 1 - floor(dim * (z_t + 0.5))      [0-based]
     *  which places the largest fractional parts in bucket 0.
     *  The index is clamped into [0, dim-1] to absorb any floating-point edge effects.
     */
    for (bi = 0; bi < dim; bi++) {
        bucket[bi] = -1;
    }
    for (t = 0; t < dim; t++) {
        bi = dim - 1 - (int)floor((double)dim * (z[t] + 0.5));
        if (bi < 0) {
            bi = 0;
        }
        else if (bi > dim - 1) {
            bi = dim - 1;
        }
        link[t] = bucket[bi];
        bucket[bi] = t;
    }

    /*
     * Main loop (Algorithm 2, lines 9-19).
     *
     *  d_0 is the initial minimum;
     *  processing bucket bi accumulates S_{bi+1}, then d_{bi+1} is tested.
     */
    D = beta - alpha * alpha / (double)dim;
    m = 0;
    for (bi = 0; bi < dim; bi++) {
        for (t = bucket[bi]; t != -1; t = link[t]) {
            alpha = alpha - 1.0;
            beta = beta - 2.0 * z[t] + 1.0;
        }
        d = beta - alpha * alpha / (double)dim;
        if (d < D) {
            D = d;
            m = bi + 1;          /* paper's 1-based candidate index */
        }
    }

    /*
     * Build u_m: add 1 to every coordinate in buckets S_1..S_m, i.e. the
     * 0-based buckets 0..m-1 (Algorithm 2, lines 20-25).
     */
    for (bi = 0; bi < m; bi++) {
        for (t = bucket[bi]; t != -1; t = link[t]) {
            k[t] += 1.0;
        }
    }

    /* Output the lattice point itself: y = Q k = k - mean(k). */
    mu = 0.0;
    for (i = 0; i < dim; i++) 
        mu += k[i];
    mu /= (double)dim;
    for (i = 0; i < dim; i++) 
        y[i] = k[i] - mu;

    return LQ_SUCCESS;
}

/* fanstar_MCSQ.c -- end of file */
