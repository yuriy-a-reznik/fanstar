/*!
 *  \file       fanstar_MCQ.c
 *  \brief      Closest-point search in An* lattices -- using McKilliam, Clarkson, & Quinn algorithm.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements closest-point finding algorithm for An* lattices by following McKilliam, Clarkson, & Quinn paper:
 *
 *   R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An Algorithm to Compute the Nearest Point in the Lattice A*_n",
 *   IEEE Trans. Inform. Theory, vol. 54, no. 9, pp. 4378-4381, Sept. 2008.
 *
 *  The lattice A*_n is the projection of the cubic lattice Z^(n+1) onto hyperplane orthogonal to 1:
 *
 *     A*_n = { Q k : k in Z^(n+1) },   Q = I - (1/(n+1)) 1 1^T.
 *
 *  The idea of the algorithm (Section IV of the paper):
 *
 *  By Lemma 2, a nearest point Q k to x corresponds to a nearest integer point k to x + xi*1 for
 *  some scalar xi, and as xi sweeps [0,1) the rounded vector takes only n+1 distinct relevant values.
 *  Sorting the centred fractional parts {x_i} = x_i - round_cs(x_i) in DESCENDING order, these candidates are
 *
 *     u_0 = round_cs(x),   u_i = u_{i-1} + e_{s_i},   i = 1..n,
 *
 *  i.e. u_i increments the i coordinates with the largest fractional parts.
 *  The candidate that minimises
 *
 *     d_i = || Q (x - u_i) ||^2 = beta_i - alpha_i^2 / (n+1)
 *
 *  gives a nearest point x_hat = Q u_m.
 * 
 *  Variables z_i = x - u_i, alpha_i = z_i^T 1 and beta_i = z_i^T z_i are updated in O(1):
 *
 *     alpha_i = alpha_{i-1} - 1,
 *     beta_i  = beta_{i-1} - 2 {x_{s_i}} + 1.
 *
 *  The whole algorithm needs a single sort, hence its complexity is O(n log n); 
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
 *  \brief Sort the index array s (length len) so that key[s[0]] >= key[s[1]] >= ... (DESCENDING).
 *
 *  Insertion sort is used deliberately: the intended use is small n (n = 2..8, so len = n+1 = 3..9).
 *  For such sizes insertion sort is faster in practice than any O(n log n) method (no recursion, minimal overhead,
 *  good branch/cache behaviour) while staying simple and stable.
 * 
 *  \param[in] s    Array of indices to sort (length len).
 *  \param[in] key  Array of keys to sort by (length len).
 *  \param[in] len  Length of the arrays.
 */
static void sort_indices_desc(int* s, const double* key, int len)
{
    int i, j, tmp;

    for (i = 0; i < len; i++) 
        s[i] = i;

    for (i = 1; i < len; i++) {
        tmp = s[i];
        j = i - 1;
        while (j >= 0 && key[s[j]] < key[tmp]) {
            s[j + 1] = s[j];
            j--;
        }
        s[j + 1] = tmp;
    }
}

/*!
 *  \brief Closest point in An* finding by using MCQ algorithm (McKilliam, Clarkson, & Quinn, 2008).
 *
 *  Implementinng algorithm of Section IV of the paper, by following notations in formulae and the text.
 *  However, x and y are switched: we assume that x is input and y is output.
 * 
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest An* point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Anstar_MCQ(const double* x, double* y, int n)
{
    double u[LQ_MAX_M];     /* u_0 = round_cs(x), later incremented to u_m         */
    double z[LQ_MAX_M];     /* centred fractional parts z_i = x_i - round_cs(x_i)  */
    int    s[LQ_MAX_M];     /* permutation sorting z[] in descending order      */

    double alpha, beta;     /* z_i^T 1 and z_i^T z_i, updated recursively       */
    double d, D;            /* current candidate distance and running minimum   */
    double mu;              /* mean, used by the projection Q u = u - mean(u)   */
    int    dim = n + 1;     /* ambient dimension n + 1                          */
    int    i, m;

    /* check arguments: */
    if (x == NULL || y == NULL || n < 2 || n > LQ_MAX_N)
        return LQ_INVARG;

    /* 
     * Compute u_0 = round_cs(x) and centred fractional parts: z = x - u_0.
     * The residual lies in [-1/2, 1/2).
     */
    for (i = 0; i < dim; i++) {
        u[i] = round_cs(x[i]);         /* replacing floor(x[i] + 0.5); for consistency with C&S */
        z[i] = x[i] - u[i];
    }

    /* Sort so that z[s[0]] >= z[s[1]] >= ... >= z[s[n]] (descending, condition (3) of the paper). */
    sort_indices_desc(s, z, dim);

    /*
     * Compute d_0 = || Q z_0 ||^2 with z_0 = z (the i = 0 candidate u_0).
     */
    alpha = 0.0;
    beta = 0.0;
    for (i = 0; i < dim; i++) {
        alpha += z[i];
        beta += z[i] * z[i];
    }
    D = beta - alpha * alpha / (double)dim;
    m = 0;

    /*
     * The main loop: 
     *
     * Recursively evaluate d_i for i = 1..n and keep the minimising index.
     *   alpha_i = alpha_{i-1} - 1
     *   beta_i  = beta_{i-1} - 2 {x_{s_i}} + 1
     * (s[i-1] is the 0-based index of the i-th largest fractional part.)
     */
    for (i = 1; i <= n; i++) {
        alpha = alpha - 1.0;
        beta = beta - 2.0 * z[s[i - 1]] + 1.0;
        d = beta - alpha * alpha / (double)dim;
        if (d < D) {
            D = d;
            m = i;
        }
    }

    /*
     * Compute u_m = round_cs(x) + e_{s_1} + ... + e_{s_m}: add 1 to the m coordinates
     * with the largest fractional parts.
     */
    for (i = 0; i < m; i++) {
        u[s[i]] += 1.0;
    }

    /* Output the lattice point itself: y = Q u_m = u_m - mean(u_m). */
    mu = 0.0;
    for (i = 0; i < dim; i++) 
        mu += u[i];
    mu /= (double)dim;
    for (i = 0; i < dim; i++)
        y[i] = u[i] - mu;

    return LQ_SUCCESS;
}

/* fanstar_MCQ.c -- end of file */
