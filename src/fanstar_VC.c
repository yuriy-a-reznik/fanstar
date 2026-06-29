/*!
 *  \file       fanstar_CS.c
 *  \brief      Closest-point search in An* lattices -- using Vaughan & Clarkson algorithm.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements closest-point finding algorithm for lattice An* by following Vaughan & Clarkson paper:
 *
 *   I. Vaughan L. Clarkson, "An Algorithm to Compute a Nearest Point in the
 *   Lattice A*_n", AAECC-13, LNCS 1719, pp. 104-120, 1999.
 *
 *  The lattice A*_n has rank n and lives in the hyperplane { u : sum(u_i) = 0 }
 *  of R^(n+1).  Its points are v = B z, with z in Z^(n+1) and the symmetric
 *  projection matrix B = I - (1/(n+1)) 1 1^T.
 *
 *  The nearest point is found by chaining the three algorithms of the paper:
 *
 *   Algorithm 1 : z_i = round_cs((x_i)              -> B z is alpha-close to B x
 *   Algorithm 2 : two while-loops on sorted d   -> B z is beta-close
 *   Algorithm 3 : prefix-sum minimisation       -> B z is gamma-close (nearest)
 *
 *  where d = B x - B z and "round_cs(" returns a nearest integer.
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
 *  \brief Sort the index array s (length len) so that key[s[0]] <= key[s[1]] <= ...
 *
 *  Insertion sort is used deliberately: the intended use is small n (n = 2..8,
 *  so len = n+1 = 3..9).  For such sizes insertion sort is faster in practice
 *  than any O(n log n) method (no recursion, no overhead, excellent cache and
 *  branch behaviour) while remaining simple and stable.
 * 
 *  \param[in] s    Array of indices to sort (length len).
 *  \param[in] key  Array of keys to sort by (length len).
 *  \param[in] len  Length of the arrays.
 */
static void sort_indices(int* s, const double* key, int len)
{
    int i, j, tmp;

    for (i = 0; i < len; i++)
        s[i] = i;

    for (i = 1; i < len; i++) {
        tmp = s[i];
        j = i - 1;
        while (j >= 0 && key[s[j]] > key[tmp]) {
            s[j + 1] = s[j];
            j--;
        }
        s[j + 1] = tmp;
    }
}

/*!
 *  \brief Closest point in An* finding by using VC algorithm (Vaughan & Clarkson, 1999).
 *
 *  The nearest point is found by chaining the three algorithms of the paper:
 *
 *   Algorithm 1 : z_i = round_cs((x_i)              -> B z is alpha-close to B x
 *   Algorithm 2 : two while-loops on sorted d   -> B z is beta-close
 *   Algorithm 3 : prefix-sum minimisation       -> B z is gamma-close (nearest)
 *
 *   where d = B x - B z and "round_cs(" returns a nearest integer.
 * 
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest An* point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Anstar_VC(const double* x, double* y, int n)
{
    double z[LQ_MAX_M];         /* integer lattice coordinates, stored as doubles  */
    double delta[LQ_MAX_M];     /* delta = B x - B z                               */
    int    s[LQ_MAX_M];         /* sorting permutation of the indices              */

    double mu;                  /* mean, used by the projection B u = u - mean(u)  */
    double t, tau;              /* running prefix sum and its minimum (Algorithm 3)*/
    int    dim = n + 1;         /* ambient dimension n + 1                         */
    int    i, k, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n < 2 || n > LQ_MAX_N) 
        return LQ_INVARG;

    /*
     * Algorithm 1: z_i = round_cs((x_i). 
     * 
     *  This yields a nearest integer with residual x_i - z_i in [-1/2, 1/2), 
     *  which is exactly the property used in the proof of Proposition 1.  
     *  The resulting B z is alpha-close to B x.
     */
    for (i = 0; i < dim; i++) 
        z[i] = round_cs(x[i]);            /* use C&S rounding, instead of floor(x[i] + 0.5); */

    /*
     * Algorithm 2: turn the alpha-close point into a beta-close point.
     *
     *  delta = project(x) - project(z) = B(x - z); element i equals
     *  (x_i - z_i) minus the mean of (x - z).
     */
    mu = 0.0;
    for (i = 0; i < dim; i++) {
        delta[i] = x[i] - z[i];
        mu += delta[i];
    }
    mu /= (double)dim;
    for (i = 0; i < dim; i++) 
        delta[i] -= mu;

    sort_indices(s, delta, dim);

    /* First while-loop (paper lines 7-10): fixes delta_{s_1} < -1/2.
     * 1-based condition  delta_{s_{m+1}} < m/(n+1) - 1/2  with body
     * m++, z_{s_m} -= 1, becomes (0-based, with loop variable m): */
    m = 0;
    while (m < dim &&
        delta[s[m]] < ((double)m) / ((double)dim) - 0.5) {
        z[s[m]] -= 1.0;
        m++;
    }

    /* Second while-loop (paper lines 12-15): fixes delta_{s_{n+1}} > 1/2.
     * 1-based condition  delta_{s_m} > m/(n+1) - 1/2  with body
     * z_{s_m} += 1, m--, starting from m = n+1, becomes (0-based): */
    m = dim;
    while (m > 0 &&
        delta[s[m - 1]] > ((double)m) / ((double)dim) - 0.5) {
        z[s[m - 1]] += 1.0;
        m--;
    }

    /*
     * Algorithm 3: turn the beta-close point into a gamma-close (nearest) point.  
     * 
     *  Recompute delta = B(x - z) for the updated z and re-sort.
     */
    mu = 0.0;
    for (i = 0; i < dim; i++) {
        delta[i] = x[i] - z[i];
        mu += delta[i];
    }
    mu /= (double)dim;
    for (i = 0; i < dim; i++) 
        delta[i] -= mu;

    sort_indices(s, delta, dim);

    /*
     * Find m = argmin_j  sum_{i=1..j} (delta_{s_i} - p_i),  where
     *
     *     p_i = (i-1)/(n+1) - n/(2(n+1)).
     *
     *  (Lemma 6 / Corollary 1 / Lemma 4.)  Note that the increment used here
     *  is delta_{s_i} - p_i.  The pseudocode printed in the paper accumulates
     *  delta_{s_i} + p_i, which is inconsistent with its own Lemma 4 and yields
     *  the wrong m (it would move an already-nearest point); the sign here is
     *  the one that makes all prefix sums non-negative after the rotation, as
     *  required by Lemma 4, and it is confirmed by brute-force testing.
     *
     *  With 0-based k = i-1:  p = k/dim - (dim-1)/(2*dim).
     */
    t = 0.0;
    tau = 0.0;
    m = 0;
    for (k = 0; k < dim - 1; k++) {     /* paper: i = 1..n */
        t += delta[s[k]]
            - ((double)k) / ((double)dim)
            + ((double)(dim - 1)) / (2.0 * (double)dim);
        if (t < tau) {
            tau = t;
            m = k + 1;                /* 1-based count of indices to adjust */
        }
    }

    for (k = 0; k < m; k++) 
        z[s[k]] -= 1.0;

    /* Output the lattice point itself: y = B z = z - mean(z): */
    mu = 0.0;
    for (i = 0; i < dim; i++) 
        mu += z[i];
    mu /= (double)dim;
    for (i = 0; i < dim; i++) 
        y[i] = z[i] - mu;

    return LQ_SUCCESS;
}

/* fanstar_VC.c -- end of file */
