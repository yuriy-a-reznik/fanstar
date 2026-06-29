/*!
 *  \file       fanstar_CS.c
 *  \brief      Closest-point search in Zn, An, An*, Dn, Dn*, and E8 lattices -- using Conway & Sloan algorithms.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module provides reference implementation of closest-point finding algorithms based on Conway & Sloan's paper:
 *
 *     [CS]  J.H. Conway and N.J.A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers and codes,"
 *           IEEE Trans. Inform. Theory, vol. IT-28, no. 2, pp. 227-232, March 1982.
 */

#include <stdlib.h>
#include <math.h>
#include "fanstar.h"

/************
 * 
 *  Rounding functions (Conway & Sloane, [1, section III]):
 *
 *    f(x)  - round x to the nearest integer, with ties broken toward zero
 *    w(x)  - round x in the opposite direction from f ("wrong-way" rounding)
 */

/*!
 *  \brief Round x to the nearest integer, with ties broken toward zero.
 *
 *  Examples: f(0.5) = 0, f(-0.5) = 0, f(1.5) = 1, f(-1.5) = -1.
 *
 *  \param[in] x  Real input.
 *  \return       The closest integer to x, returned as a double.
 */
static double f(double x)
{
    double a = fabs(x), r = floor(a + 0.5);
    if (r - a >= 0.5) r -= 1.0;             /* exact tie -> toward zero */
    return (x >= 0.0) ? r : -r;
}

/*!
 *  \brief "Wrong-way" rounding - rounds in the opposite direction from f.
 *
 *  If delta(x) = x - f(x) > 0  then w(x) = f(x) + 1.
 *  If delta(x) < 0             then w(x) = f(x) - 1.
 *  If delta(x) = 0 (integer x) then w(x) = f(x) + 1, matching w(0) = 1
 * 
 *  \param[in] x  Real input.
 *  \return       The wrong-way rounded integer for x, returned as a double.
 */
static double w(double x)
{
    double fx = f(x), d = x - fx;
    if (d > 0.0) return fx + 1.0;
    if (d < 0.0) return fx - 1.0;
    return fx + 1.0;                 /* integer x: +1 (paper) */
}

/**********
 *
 *  Closest point algorithms:
 *
 *  closest_point_Zn_CS()      - closest point for lattice Zn
 *  closest_point_Dn_CS()      - closest point for lattice Dn
 *  closest_point_An_CS()      - closest point for lattice An
 *  closest_point_Dnstar_CS()  - closest point for lattice Dn*
 *  closest_point_Anstar_CS()  - closest point for lattice An*
 *  closest_point_E8_CS()      - closest point for lattice E8
*/

/*!
 *  \brief Closest point algorithm for lattice Z^n.
 *
 *  \param[in]  x    Input vector, n doubles.
 *  \param[out] out  Closest z_n point, n doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Zn_CS(const double* x, double* out, int n)
{
    int i;

    /* check arguments: */
    if (x == NULL || out == NULL)      return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)         return LQ_INVARG;

    /* round to Zn: */
    for (i = 0; i < n; i++) out[i] = f(x[i]);

    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point  algorithm for lattice Dn = { x in Z^n : sum x_i is even }.
 *
 *  Based on Conway & Sloan Alg. 2.
 *  Whichever of f(x) and g(x) has even sum is returned,
 *  where g(x) differs from f(x) in exactly one coordinate -- the one with the
 *  largest |delta(x_i)| -- which is replaced by w(x_k).
 *
 *  \param[in]  x    Input vector, n doubles.
 *  \param[out] out  Closest D_n point, n doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Dn_CS(const double* x, double* out, int n)
{
    int i, isum, k = 0;
    double sum = 0.0, best = -1.0, d;

    /* check arguments: */
    if (x == NULL || out == NULL)      return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)         return LQ_INVARG;

    /* find closest to Z^n, and track the worst component: */
    for (i = 0; i < n; i++) {
        out[i] = f(x[i]);
        sum += out[i];
        d = fabs(x[i] - out[i]);
        if (d > best) { best = d; k = i; }
    }
    /* if sum is odd, round the worst component the wrong way: */
    isum = (int)sum;
    if (isum & 1L) out[k] = w(x[k]);
    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for lattice An = { x in Z^{n+1} : sum x_i = 0 }.
 *
 *  Based on Conway & Sloan Alg. 5.
 *  Buffers x and out hold n+1 doubles; the lattice lives on the sum=0 hyperplane in R^{n+1}.
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest A_n point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_An_CS(const double* x, double* out, int n)
{
    int n1, idx[LQ_MAX_M];
    double xp[LQ_MAX_M], dx[LQ_MAX_M], s = 0.0;
    int i, j, t, Delta = 0;

    /* check arguments: */
    if (x == NULL || out == NULL)      return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)         return LQ_INVARG;

    /* compute s = mean(x) */
    n1 = n + 1;
    for (i = 0; i < n1; i++) s += x[i];
    s /= (double)n1;
 
    /* compute nearest-point quantized centered coordinates f(xp[]) and discrepancy Delta: */
    for (i = 0; i < n1; i++) {
        xp[i] = x[i] - s;
        out[i] = f(xp[i]);
        dx[i] = xp[i] - out[i];
        Delta += (int)out[i];
        idx[i] = i;
    }

    /* insertion sort idx by dx ascending; n+1 <= 9 elements */
    for (i = 1; i < n1; i++) {
        t = idx[i];
        for (j = i; j > 0 && dx[idx[j - 1]] > dx[t]; j--)
            idx[j] = idx[j - 1];
        idx[j] = t;
    }

    /* if Delta > 0, we need to round some components down; if Delta < 0, we need to round some components up. */
    if (Delta > 0)
        for (j = 0; j < Delta; j++)        out[idx[j]] -= 1.0;
    else if (Delta < 0)
        for (j = 0; j < -Delta; j++)       out[idx[n1 - 1 - j]] += 1.0;

    return LQ_SUCCESS;
}

/* Unified closest-point function type */
typedef int (*closest_point_t)(const double* x, double* y, int n);

/*!
 *  \brief Find closest point of a union of cosets (R[0] + L) U (R[1] + L) U ... U (R[d-1] + L)
 *         given a closest-point routine for the base lattice L.  (Conway/Sloane Algorithm 3.)
 *
 *  For each j we form xs = x - R[j], call cpf(xs, n, y), restore y += R[j],
 *  and keep the y with smallest squared Euclidean distance to x.
 *
 *  \param[in]  x              Input vector, m doubles.
 *  \param[in]  m              Ambient dimension of the vectors.
 *  \param[in]  d              Number of cosets.
 *  \param[in]  R              d * m doubles, row-major. Coset representatives.
 *  \param[in]  find_closest   Closest-point routine for the base lattice L.
 *  \param[in]  n              Base lattice-dimension argument forwarded to find_closest().
 *  \param[out] out            Closest lattice point, m doubles.
 *
 *  \return     LQ_SUCCESS or LQ_INVARG if x, R, or out is NULL or if m, d, or n are out of range.
 */
static int closest_point_coset(const double* x, int m, int d, double* R, closest_point_t find_closest, int n, double* out)
{
    double xs[LQ_MAX_M], y[LQ_MAX_M];
    double bd = -1.0, dd, t;
    int i, j;

    /* check arguments: */
    if (x == NULL || R == NULL || find_closest == NULL || out == NULL) return LQ_INVARG;
    if (m < 1 || m > LQ_MAX_N + 1) return LQ_INVARG;
    if (d < 1 || d > LQ_MAX_N + 1) return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)     return LQ_INVARG;

    /* scan cosets: */
    for (j = 0; j < d; j++) {
        /* form xs = x - R[j]: */
        for (i = 0; i < m; i++) xs[i] = x[i] - R[j * m + i];
        /* find closest point of L to xs: */
        find_closest((const double *)xs, y, n);
        /* restore y += R[j] and compute distance to x: */
        dd = 0.0;
        for (i = 0; i < m; i++) {
            y[i] += R[j * m + i];
            t = x[i] - y[i];
            dd += t * t;
        }
        /* keep best: */
        if (bd < 0.0 || dd < bd) {
            bd = dd;
            for (i = 0; i < m; i++) out[i] = y[i];
        }
    }

    return LQ_SUCCESS;
}

/*!
 *  \brief Closest point algorithm for lattice Dn*.
 *
 *  \param[in]  x    Input vector, n doubles.
 *  \param[out] out  Closest D_n* point, n doubles.
 *  \param[in]  n    Lattice dimension.
 * 
 *  \return LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Dnstar_CS(const double* x, double* out, int n)
{
    double R[2 * LQ_MAX_M];
    int i;

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)     return LQ_INVARG;

    /* form coset representatives: 0 and (0.5,...,0.5): */
    for (i = 0; i < n; i++) { R[i] = 0.0; R[n + i] = 0.5; }

    /* find closest point of the union of the two cosets: */
    return closest_point_coset(x, n, 2, R, closest_point_Zn_CS, n, out);
}

/*!
 *  \brief Closest point algorithm for lattice An*.
 *
 *  An* realized as the union of n+1 cosets of An with reps r_i (C&S, eq. (3)).
 *
 *  \param[in]  x    Input vector, n+1 doubles.
 *  \param[out] out  Closest An* point, n+1 doubles.
 *  \param[in]  n    Lattice dimension.
 *
 *  \return LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n is out of range.
 */
int closest_point_Anstar_CS(const double* x, double* out, int n)
{
    double R[LQ_MAX_M * LQ_MAX_M], a, b;
    int n1, i, j;

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n < 1 || n > LQ_MAX_N)     return LQ_INVARG;

    /* form coset representatives: r_i has i components -((n+1-i)/n+1) and n+1-i components (i/n+1): */
    n1 = n + 1;
    for (i = 0; i <= n; i++) {
        a = -((double)(n + 1 - i)) / (double)n1;
        b = ((double)i) / (double)n1;
        for (j = 0; j < i; j++) R[i * n1 + j] = a;
        for (j = i; j < n1; j++) R[i * n1 + j] = b;
    }

    /* find closest point of the union of the n+1 cosets: */
    return closest_point_coset(x, n1, n + 1, R, closest_point_An_CS, n, out);
}

/*!
 *  \brief Closest point algorithm for the Gosset lattice E8 = D8 U ((1/2)^8 + D8).
 *
 *  Based on Conway & Sloan, Alg. 4.
 *
 *  \param[in]  x    Input vector, 8 doubles.
 *  \param[out] out  Closest E_8 point, 8 doubles.
 *  \param[in]  n    Lattice dimension; must equal 8.
 *
 *  \return     LQ_SUCCESS, or LQ_INVARG if x or out is NULL or n != 8.
 */
int closest_point_E8_CS(const double* x, double* out, int n)
{
    static double R[2 * 8] = {
        0,0,0,0,0,0,0,0,
        0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5
    };

    /* check arguments: */
    if (x == NULL || out == NULL)  return LQ_INVARG;
    if (n != 8)                    return LQ_INVARG;
    
    /* find closest point of the union of the two cosets: */
    return closest_point_coset(x, 8, 2, R, closest_point_Dn_CS, 8, out);
}

/* fanstar_CS.c -- end of file */
