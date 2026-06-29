/*!
 *  \file       fanstar_MCQ_opt.c
 *  \brief      Closest-point search in An lattices -- unrolled and optimized version of McKilliam, Clarkson, & Quinn algorithm.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module contains optimized version of McKilliam, Clarkson, & Quinn algorithm for An* lattices for n=2..8.
 * 
 *  Implementation notes:
 *
 *  Pre-condition: x lies in the ambient (n+1)-space with x[0] + x[1] + ... + x[n] = 0.
 * 
 *  Per-n algorithm implementations:
 *    size-optimal branchless sorting networks + linear prefix scans.
 *    everything is unrolled, no loops.
 * 
 *  Zero-sum simplifications used throughout:
 *    (a) alpha = sum_j (x_j - k_j) = -sum_j k_j is an EXACT INTEGER. 
 *        We accumulate the integer round-sum sumk = sum_j k_j  during the rounding pass 
 *        and set alpha = -(double)sumk. This is cheaper than a floating-point residual 
 *        sum and free of accumulation error.
 *    (b) The output projection mean = (sum of final k)/N is obtained as
 *        (sumk + flips)/N, so no second summation of k is required.
 *    (c) x is already its own projection (Qy = y), so the input is never projected.
 *
 *  Common pipeline (per function):
 *    1. k_j   = round(x_j);  z_j = x_j - k_j in [-1/2,1/2);  sumk += k_j.
 *       alpha = -(double)sumk    (exact integer).
 *    2. Find the optimal nested-candidate prefix that minimises
 *          h = beta_var - (alpha - cnt)^2 / N
 *       where beta_var accumulates (1 - 2 z) over the flipped coordinates
 *       (the constant z'z offset of the papers is dropped).
 *       MCQ: sort residuals descending, scan i=1..n (cnt = i).
 *    3. Flip the chosen coordinates' k upward by 1; flips counts them.
 *    4. mean = (double)(sumk + flips)/N;   y_j = (double)k_j - mean.
 *
 *   Rounding control:
 *       For matching C&S implementation exactly, please #define LQ_ROUND_CS 1.
 *       If LQ_ROUND_CS=0, the rounding is implemented using half-way-up process, which is
 *       faster. Both methods produce equidistant points from x.
 */

#include <math.h>
#include "fanstar.h"

/*!
 *  Rounding control options:
 *    LQ_ROUND_CS = 1 -> round to the nearest integer with ties broken toward zero (Conway & Sloan 1982 convention).
 *    LQ_ROUND_CS = 0 -> round to the nearest integer with ties rounded up (faster execution).
 */
#define LQ_ROUND_CS 0

/*! 
 *  \def   ROUND
 *  \brief Rounding to the nearest integer macro.
 */
#if LQ_ROUND_CS
#define ROUND(x) ((x) >= 0.0 ? ceil((x)  - 0.5): floor((x) + 0.5))  /* ties broken toward zero */
#else
#define ROUND(x) floor((x) + 0.5)                                   /* ties rounded up */
#endif

/*!
 *  \def   RND_ACC
 *  \brief Round x[j] to k[j], store residual z[j], accumulate integer sumk.
 *         Exploits the zero-sum precondition: alpha = -(double)sumk later.
 */
#define RND_ACC(j, x, k, z, sumk)                                     \
    do {                                                              \
        k[j] = (long)ROUND(x[j]);                                     \
        z[j] = x[j] - (double)k[j];                                   \
        sumk += k[j];                                                 \
    } while (0)

/*!
 *  \def CMPXCHG
 *  \brief Descending compare-exchange on parallel value/index arrays.
 *         After execution v[p] >= v[q].  Compilers lower this to branchless
 *         min/max + blends, giving a fully predictable sorting network.
 */
#define CMPXCHG(v, ix, p, q)                                          \
    do {                                                              \
        if ((v)[q] > (v)[p]) {                                        \
            double _tv = (v)[p]; (v)[p] = (v)[q]; (v)[q] = _tv;       \
            int    _ti = (ix)[p]; (ix)[p] = (ix)[q]; (ix)[q] = _ti;   \
        }                                                             \
    } while (0)

/*!
 *  \def MCQ_SCAN
 *  \brief One MCQ prefix-scan step at index i (cnt == i implicitly).
 *         w is the i-th largest residual. alpha is the exact integer.
 *         Updates running minimum (best, m).
 */
#define MCQ_SCAN(i, w, INV, alpha, acc, best, m)                     \
    do {                                                             \
        double _a, _h;                                               \
        acc += 1.0 - 2.0*(w);                                        \
        _a = (alpha) - (double)(i);                                  \
        _h = acc - _a*_a*(INV);                                      \
        if (_h < (best)) { (best) = _h; (m) = (i); }                 \
    } while (0)


/****************
 * 
 *  Closest-point functions for lattices A_2^* ... A_6^*.
 * 
 ****/

/*!
 *  \brief Closest point of A2* (MCQ algorithm).
 * 
 *   A2*:   N = 3,  optimal network = 3 compare-exchanges
 * 
 *  \param[in]  x  3-vector, zero-sum.
 *  \param[out] y  3-vector, closest A2* point (zero-sum).
 */
static void closest_point_a2star(const double* x, double* y)
{
    long   k[3];  double v[3];  int ix[3];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 3.0;

    /* round, compute residuals, and sum */
    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    alpha = -(double)sumk;                      /* alpha = integer */

    /* sort residuals: */
    CMPXCHG(v, ix, 0, 2); CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 1, 2);

    /* find prefix that minimises h = beta_var - (alpha - cnt)^2 / N: */
    best = -(alpha * alpha) * INV;              /* h_0 (flip nothing) */
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);

    /* flip choosen vector coordinates upward by 1; flips counts them. */
    switch (m) {
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    /* produce lattice point: */
    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
}

/*!
 *  \brief Closest point of A3* (MCQ algorithm).
 * 
 *   A3*:   N = 4,  optimal network = 5 compare-exchanges
 *
 *  \param[in]  x  4-vector, zero-sum.
 *  \param[out] y  4-vector, closest A3* point (zero-sum).
 */
static void closest_point_A3star(const double* x, double* y)
{
    long   k[4];  double v[4];  int ix[4];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 4.0;

    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    alpha = -(double)sumk;

    CMPXCHG(v, ix, 0, 2); CMPXCHG(v, ix, 1, 3); CMPXCHG(v, ix, 0, 1);
    CMPXCHG(v, ix, 2, 3); CMPXCHG(v, ix, 1, 2);

    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);

    switch (m) {
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean;
    y[2] = (double)k[2] - mean; y[3] = (double)k[3] - mean;
}

/*!
 *  \brief Closest point of A4* (MCQ algorithm).
 *
 *   A4*:   N = 5,  optimal network = 9 compare-exchanges
 *
 *  \param[in]  x  5-vector, zero-sum.
 *  \param[out] y  5-vector, closest A4* point (zero-sum).
 */
static void closest_point_A4star(const double* x, double* y)
{
    long   k[5];  double v[5];  int ix[5];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 5.0;

    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    RND_ACC(4, x, k, v, sumk); ix[4] = 4;
    alpha = -(double)sumk;

    /* 9-op sorting network: (0,3) (1,4) (0,2) (1,3) (0,1) (2,4) (1,2) (3,4) (2,3) */
    CMPXCHG(v, ix, 0, 3); CMPXCHG(v, ix, 1, 4); CMPXCHG(v, ix, 0, 2);
    CMPXCHG(v, ix, 1, 3); CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 2, 4);
    CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 2, 3);

    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);
    MCQ_SCAN(4, v[3], INV, alpha, acc, best, m);

    switch (m) {
    case 4: k[ix[3]]++;
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean;
}

/*!
 *  \brief Closest point of A5* (MCQ algorithm).
 *
 *   A5*:   N = 6,  optimal network = 12 compare-exchanges
 *
 *  \param[in]  x  6-vector, zero-sum.
 *  \param[out] y  6-vector, closest A5* point (zero-sum).
 */
static void closest_point_A5star(const double* x, double* y)
{
    long   k[6];  double v[6];  int ix[6];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 6.0;

    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    RND_ACC(4, x, k, v, sumk); ix[4] = 4;
    RND_ACC(5, x, k, v, sumk); ix[5] = 5;
    alpha = -(double)sumk;

    /* 12-op sorting network: (0,5) (1,3) (2,4) (1,2) (3,4) (0,3) (2,5) (0,1) (2,3) (4,5) (1,2) (3,4) */
    CMPXCHG(v, ix, 0, 5); CMPXCHG(v, ix, 1, 3); CMPXCHG(v, ix, 2, 4);
    CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 0, 3);
    CMPXCHG(v, ix, 2, 5); CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 2, 3);
    CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4);

    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);
    MCQ_SCAN(4, v[3], INV, alpha, acc, best, m);
    MCQ_SCAN(5, v[4], INV, alpha, acc, best, m);

    switch (m) {
    case 5: k[ix[4]]++;
    case 4: k[ix[3]]++;
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
}

/*!
 *  \brief Closest point of A6* (MCQ algorithm).
 *
 *   A6*:   N = 7,  optimal network = 16 compare-exchanges
 *
 *  \param[in]  x  7-vector, zero-sum.
 *  \param[out] y  7-vector, closest A6* point (zero-sum).
 */
static void closest_point_A6star(const double* x, double* y)
{
    long   k[7];  double v[7];  int ix[7];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 7.0;

    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    RND_ACC(4, x, k, v, sumk); ix[4] = 4;
    RND_ACC(5, x, k, v, sumk); ix[5] = 5;
    RND_ACC(6, x, k, v, sumk); ix[6] = 6;
    alpha = -(double)sumk;

    /* 16-op sorting network: 
     *  (0,6) (2,3) (4,5) (0,2) (1,4) (3,6) (0,1) (2,5) (3,4) (1,2) (4,6) (2,3) (4,5) (1,2) (3,4) (5,6) */
    CMPXCHG(v, ix, 0, 6); CMPXCHG(v, ix, 2, 3); CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 0, 2); 
    CMPXCHG(v, ix, 1, 4); CMPXCHG(v, ix, 3, 6); CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 2, 5); 
    CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 4, 6); CMPXCHG(v, ix, 2, 3);
    CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 5, 6);

    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);
    MCQ_SCAN(4, v[3], INV, alpha, acc, best, m);
    MCQ_SCAN(5, v[4], INV, alpha, acc, best, m);
    MCQ_SCAN(6, v[5], INV, alpha, acc, best, m);

    switch (m) {
    case 6: k[ix[5]]++;
    case 5: k[ix[4]]++;
    case 4: k[ix[3]]++;
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
    y[6] = (double)k[6] - mean;
}

/*!
 *  \brief Closest point of A7* (MCQ algorithm).
 *
 *   A7*:   N = 8   19-comparator network
 *
 *  \param[in]  x  8-vector, zero-sum.
 *  \param[out] y  8-vector, closest A6* point (zero-sum).
 */
static void closest_point_A7star(const double* x, double* y)
{
    long   k[8];  double v[8];  int ix[8];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 8.0;

	/* round, compute residuals, and sum */
    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    RND_ACC(4, x, k, v, sumk); ix[4] = 4;
    RND_ACC(5, x, k, v, sumk); ix[5] = 5;
    RND_ACC(6, x, k, v, sumk); ix[6] = 6;
    RND_ACC(7, x, k, v, sumk); ix[7] = 7;
    alpha = -(double)sumk;

    /* 8-element descending sort (19 CE, Batcher):
     * (0,2)(1,3)(4,6)(5,7)(0,4)(1,5)(2,6)(3,7)(0,1)(2,3) (4,5)(6,7)(2,4)(3,5)(1,4)(3,6)(1,2)(3,4)(5,6). */
    CMPXCHG(v, ix, 0, 2); CMPXCHG(v, ix, 1, 3); CMPXCHG(v, ix, 4, 6); CMPXCHG(v, ix, 5, 7);
    CMPXCHG(v, ix, 0, 4); CMPXCHG(v, ix, 1, 5); CMPXCHG(v, ix, 2, 6); CMPXCHG(v, ix, 3, 7);
    CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 2, 3); CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 6, 7);
    CMPXCHG(v, ix, 2, 4); CMPXCHG(v, ix, 3, 5); CMPXCHG(v, ix, 1, 4); CMPXCHG(v, ix, 3, 6);
    CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 5, 6);

    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);
    MCQ_SCAN(4, v[3], INV, alpha, acc, best, m);
    MCQ_SCAN(5, v[4], INV, alpha, acc, best, m);
    MCQ_SCAN(6, v[5], INV, alpha, acc, best, m);
    MCQ_SCAN(7, v[6], INV, alpha, acc, best, m);

    switch (m) {
    case 7: k[ix[6]]++;
    case 6: k[ix[5]]++;
    case 5: k[ix[4]]++;
    case 4: k[ix[3]]++;
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    /* compute the output lattice point */
    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
    y[6] = (double)k[6] - mean; y[7] = (double)k[7] - mean;
}

/*!
 *  \brief Closest point of A8* (MCS algorithm).
 *
 *   A8*:   N = 9   25-comparator network
 *
 *  \param[in]  x  9-vector, zero-sum.
 *  \param[out] y  9-vector, closest A8* point (zero-sum).
 */
static void closest_point_A8star(const double* x, double* y)
{
    long   k[9];  double v[9];  int ix[9];
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int m = 0;
    const double INV = 1.0 / 9.0;

    /* round, compute residuals, and sum */
    RND_ACC(0, x, k, v, sumk); ix[0] = 0;
    RND_ACC(1, x, k, v, sumk); ix[1] = 1;
    RND_ACC(2, x, k, v, sumk); ix[2] = 2;
    RND_ACC(3, x, k, v, sumk); ix[3] = 3;
    RND_ACC(4, x, k, v, sumk); ix[4] = 4;
    RND_ACC(5, x, k, v, sumk); ix[5] = 5;
    RND_ACC(6, x, k, v, sumk); ix[6] = 6;
    RND_ACC(7, x, k, v, sumk); ix[7] = 7;
    RND_ACC(8, x, k, v, sumk); ix[8] = 8;
    alpha = -(double)sumk;

    /* 9-element descending sort (25 CE, size-optimal):
     *  (0,3)(1,7)(2,5)(4,8)(0,7)(2,4)(3,8)(5,6)(0,2)(1,3)(4,5)(7,8) (1,4)(3,6)(5,7)(0,1)(2,4)(3,5)(6,8)(2,3)(4,5)(6,7)(1,2)(3,4)(5,6) */
    CMPXCHG(v, ix, 0, 3); CMPXCHG(v, ix, 1, 7); CMPXCHG(v, ix, 2, 5); CMPXCHG(v, ix, 4, 8); CMPXCHG(v, ix, 0, 7);
    CMPXCHG(v, ix, 2, 4); CMPXCHG(v, ix, 3, 8); CMPXCHG(v, ix, 5, 6); CMPXCHG(v, ix, 0, 2); CMPXCHG(v, ix, 1, 3);
    CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 7, 8); CMPXCHG(v, ix, 1, 4); CMPXCHG(v, ix, 3, 6); CMPXCHG(v, ix, 5, 7);
    CMPXCHG(v, ix, 0, 1); CMPXCHG(v, ix, 2, 4); CMPXCHG(v, ix, 3, 5); CMPXCHG(v, ix, 6, 8); CMPXCHG(v, ix, 2, 3);
    CMPXCHG(v, ix, 4, 5); CMPXCHG(v, ix, 6, 7); CMPXCHG(v, ix, 1, 2); CMPXCHG(v, ix, 3, 4); CMPXCHG(v, ix, 5, 6);


    best = -(alpha * alpha) * INV;
    MCQ_SCAN(1, v[0], INV, alpha, acc, best, m);
    MCQ_SCAN(2, v[1], INV, alpha, acc, best, m);
    MCQ_SCAN(3, v[2], INV, alpha, acc, best, m);
    MCQ_SCAN(4, v[3], INV, alpha, acc, best, m);
    MCQ_SCAN(5, v[4], INV, alpha, acc, best, m);
    MCQ_SCAN(6, v[5], INV, alpha, acc, best, m);
    MCQ_SCAN(7, v[6], INV, alpha, acc, best, m);
    MCQ_SCAN(8, v[7], INV, alpha, acc, best, m);

    switch (m) {
    case 8: k[ix[7]]++;
    case 7: k[ix[6]]++;
    case 6: k[ix[5]]++;
    case 5: k[ix[4]]++;
    case 4: k[ix[3]]++;
    case 3: k[ix[2]]++;
    case 2: k[ix[1]]++;
    case 1: k[ix[0]]++;
    default: break;
    }

    /* compute the output lattice point */
    mean = (double)(sumk + m) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
    y[6] = (double)k[6] - mean; y[7] = (double)k[7] - mean; y[8] = (double)k[8] - mean;
}

/* Unified closest-point function type */
typedef void (*closest_point_Anstar_t)(const double* x, double* y);

/* Dispatch table: */
static const closest_point_Anstar_t closest_point_Anstar[] = {
    NULL, NULL, closest_point_a2star, closest_point_A3star, closest_point_A4star,
    closest_point_A5star, closest_point_A6star, closest_point_A7star, closest_point_A8star
};

/*!
 *  \brief Closest point of An* (faster method for n=2..8).
 */
int closest_point_Anstar_MCQ_opt(const double* x, double* y, int n)
{
    /* check parameters: */
	if (n < 2 || n > 8) return LQ_INVARG;

    /* call the appropriate size-optimized function: */
	closest_point_Anstar[n](x, y);

    return LQ_SUCCESS;
}

/* fanstar_MCQ_opt.c -- end of file */
