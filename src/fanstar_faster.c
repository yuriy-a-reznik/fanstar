/*!
 *  \file       fanstar_faster.c
 *  \brief      Closest-point search in An lattices -- proposed mixed algorithms method, achieving best speed at each n.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 *
 *  This module implements closest-point finding algorithms for An* lattices for n=2..8.
 *  Based on:
 *   [MCQ]     R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in 
 *             the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378–4381, Sep. 2008.
 * 
 *   [MCSQ]    R. G. McKilliam, I. V. L. Clarkson, W. D. Smith, and B. G. Quinn, "A linear-time nearest
 *             point algorithm for the lattice An*," in Proc. Int. Symp. Information Theory and its Applications
 *             (ISITA), Auckland, New Zealand, Dec. 2008, pp. 1–5.
 * 
 *   [FANSTAR] Y. Reznik, "Faster closest point algorithms for the An* lattices" -- in preparation, 2026
 *
 *  Implementation notes:
 * 
 *  Per-n method selection (best speed on x64):
 *     n = 2,3 :      modified MCQ -> size-optimal branchless sorting network + linear prefix scan.
 *     n = 4,5,6,7,8: modified MCSQ -> single-step O(n) bucketing (no comparisons).
 *
 *  Pre-condition: x lies in the ambient (n+1)-space with x[0] + x[1] + ... + x[n] = 0.
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
 *       MCSQ: bucket residuals, scan buckets 1..N (cnt tracked explicitly, no bucket lists).
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

/*!
 *  \def MCSQ_BUCKET
 *  \brief Bucket element j into b[j] = clamp(ceil(N*(1/2 - z_j)), 1, N) and
 *         scatter (1 - 2 z_j) and a unit count.  O(1), no comparisons.
 */
#define MCSQ_BUCKET(j, N, z, b, bd, bc)                              \
    do {                                                             \
        double _t = (double)(N) * (0.5 - z[j]);                      \
        int _bi = (int)ceil(_t);                                     \
        if (_bi < 1) _bi = 1; else if (_bi > (N)) _bi = (N);         \
        b[j] = _bi;                                                  \
        bd[_bi] += 1.0 - 2.0 * z[j];                                 \
        bc[_bi] += 1;                                                \
    } while (0)

/*!
 *  \def MCSQ_SCAN
 *  \brief One MCSQ scan step over bucket i. cnt tracks the true number of
 *         flipped coordinates (a bucket may hold 0, 1 or several elements).
 *         Captures fbest = cnt at the winning bucket for the projection sum.
 */
#define MCSQ_SCAN(i, bd, bc, INV, alpha, acc, cnt, best, m, fbest)   \
    do {                                                             \
        double _a, _h;                                               \
        cnt += bc[i];                                                \
        acc += bd[i];                                                \
        _a = (alpha) - (double)cnt;                                  \
        _h = acc - _a*_a*(INV);                                      \
        if (_h < (best)) { (best) = _h; (m) = (i); (fbest) = cnt; }  \
    } while (0)

/*!
 *  \def MCSQ_RECON
 *  \brief Reconstruction test: flip k[j] up by 1 iff its bucket <= m.
 */
#define MCSQ_RECON(j, b, m, k)                                       \
    do {                                                             \
        if (b[j] <= (m)) k[j]++;                                     \
    } while (0)


/****************
 * 
 *  Closest-point functions for lattices A_2^* ... A_6^*.
 * 
 ****/

/*!
 *  \brief Closest point of A2* (modified MCQ algorithm).
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
 *  \brief Closest point of A3* (modified MCQ algorithm).
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
 *  \brief Closest point of A4* (modified MCSQ algorithm).
 *
 *   A4*:   N = 5,  single-pass bucketing beats 9-comparator network
 *
 *  \param[in]  x  5-vector, zero-sum.
 *  \param[out] y  5-vector, closest A4* point (zero-sum).
 */
static void closest_point_A4star(const double* x, double* y)
{
    long   k[5];  double z[5];  int b[5];
    double bd[6] = { 0.0 };  int bc[6] = { 0 };
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int cnt = 0, m = 0, fbest = 0;
    const double INV = 1.0 / 5.0;

    /* round, compute residuals, and sum */
    RND_ACC(0, x, k, z, sumk); RND_ACC(1, x, k, z, sumk);
    RND_ACC(2, x, k, z, sumk); RND_ACC(3, x, k, z, sumk);
    RND_ACC(4, x, k, z, sumk);
    alpha = -(double)sumk;

    /* bucket residuals into b[j] = clamp(ceil(N*(1/2 - z_j)), 1, N) and scatter (1 - 2 z_j) and a unit count */
    MCSQ_BUCKET(0, 5, z, b, bd, bc); MCSQ_BUCKET(1, 5, z, b, bd, bc);
    MCSQ_BUCKET(2, 5, z, b, bd, bc); MCSQ_BUCKET(3, 5, z, b, bd, bc);
    MCSQ_BUCKET(4, 5, z, b, bd, bc);

    /* find prefix that minimises h = beta_var - (alpha - cnt)^2 / N: */
    best = -(alpha * alpha) * INV;
    MCSQ_SCAN(1, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(2, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(3, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(4, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(5, bd, bc, INV, alpha, acc, cnt, best, m, fbest);

    /* flip the chosen coordinates' k upward by 1 */
    MCSQ_RECON(0, b, m, k); MCSQ_RECON(1, b, m, k); MCSQ_RECON(2, b, m, k);
    MCSQ_RECON(3, b, m, k); MCSQ_RECON(4, b, m, k);

    /* compute the output lattice point */
    mean = (double)(sumk + fbest) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean;
}

/*!
 *  \brief Closest point of A5* (modified MCSQ algorithm).
 *
 *   A5*:   N = 6,  single-pass bucketing beats 12-comparator network
 *
 *  \param[in]  x  6-vector, zero-sum.
 *  \param[out] y  6-vector, closest A5* point (zero-sum).
 */
static void closest_point_A5star(const double* x, double* y)
{
    long   k[6];  double z[6];  int b[6];
    double bd[7] = { 0.0 };  int bc[7] = { 0 };
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int cnt = 0, m = 0, fbest = 0;
    const double INV = 1.0 / 6.0;

    /* round, compute residuals, and sum */
    RND_ACC(0, x, k, z, sumk); RND_ACC(1, x, k, z, sumk);
    RND_ACC(2, x, k, z, sumk); RND_ACC(3, x, k, z, sumk);
    RND_ACC(4, x, k, z, sumk); RND_ACC(5, x, k, z, sumk);
    alpha = -(double)sumk;

    /* bucket residuals into b[j] = clamp(ceil(N*(1/2 - z_j)), 1, N) and scatter (1 - 2 z_j) and a unit count */
    MCSQ_BUCKET(0, 6, z, b, bd, bc); MCSQ_BUCKET(1, 6, z, b, bd, bc);
    MCSQ_BUCKET(2, 6, z, b, bd, bc); MCSQ_BUCKET(3, 6, z, b, bd, bc);
    MCSQ_BUCKET(4, 6, z, b, bd, bc); MCSQ_BUCKET(5, 6, z, b, bd, bc);

    /* find prefix that minimises h = beta_var - (alpha - cnt)^2 / N: */
    best = -(alpha * alpha) * INV;
    MCSQ_SCAN(1, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(2, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(3, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(4, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(5, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(6, bd, bc, INV, alpha, acc, cnt, best, m, fbest);

    /* flip the chosen coordinates' k upward by 1 */
    MCSQ_RECON(0, b, m, k); MCSQ_RECON(1, b, m, k); MCSQ_RECON(2, b, m, k);
    MCSQ_RECON(3, b, m, k); MCSQ_RECON(4, b, m, k); MCSQ_RECON(5, b, m, k);

    /* compute the output lattice point */
    mean = (double)(sumk + fbest) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
}

/*!
 *  \brief Closest point of A6* (modified MCSQ algorithm).
 *
 *   A6*:   N = 7,  single-pass bucketing beats 16-comparator network
 *
 *  \param[in]  x  7-vector, zero-sum.
 *  \param[out] y  7-vector, closest A6* point (zero-sum).
 */
static void closest_point_A6star(const double* x, double* y)
{
    long   k[7];  double z[7];  int b[7];
    double bd[8] = { 0.0 };  int bc[8] = { 0 };
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int cnt = 0, m = 0, fbest = 0;
    const double INV = 1.0 / 7.0;

    /* round, compute residuals, and sum */
    RND_ACC(0, x, k, z, sumk); RND_ACC(1, x, k, z, sumk);
    RND_ACC(2, x, k, z, sumk); RND_ACC(3, x, k, z, sumk);
    RND_ACC(4, x, k, z, sumk); RND_ACC(5, x, k, z, sumk);
    RND_ACC(6, x, k, z, sumk);
    alpha = -(double)sumk;

    /* bucket residuals into b[j] = clamp(ceil(N*(1/2 - z_j)), 1, N) and scatter (1 - 2 z_j) and a unit count */
    MCSQ_BUCKET(0, 7, z, b, bd, bc); MCSQ_BUCKET(1, 7, z, b, bd, bc);
    MCSQ_BUCKET(2, 7, z, b, bd, bc); MCSQ_BUCKET(3, 7, z, b, bd, bc);
    MCSQ_BUCKET(4, 7, z, b, bd, bc); MCSQ_BUCKET(5, 7, z, b, bd, bc);
    MCSQ_BUCKET(6, 7, z, b, bd, bc);

    /* find prefix that minimises h = beta_var - (alpha - cnt)^2 / N: */
    best = -(alpha * alpha) * INV;
    MCSQ_SCAN(1, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(2, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(3, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(4, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(5, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(6, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(7, bd, bc, INV, alpha, acc, cnt, best, m, fbest);

    /* flip the chosen coordinates' k upward by 1 */
    MCSQ_RECON(0, b, m, k); MCSQ_RECON(1, b, m, k); MCSQ_RECON(2, b, m, k);
    MCSQ_RECON(3, b, m, k); MCSQ_RECON(4, b, m, k); MCSQ_RECON(5, b, m, k);
    MCSQ_RECON(6, b, m, k);

    /* compute the output lattice point */
    mean = (double)(sumk + fbest) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
    y[6] = (double)k[6] - mean;
}

/*!
 *  \brief Closest point of A7* (modified MCSQ algorithm).
 *
 *   A7*:   N = 8   single-pass bucketing beats 19-comparator network
 *
 *  \param[in]  x  8-vector, zero-sum.
 *  \param[out] y  8-vector, closest A6* point (zero-sum).
 */
static void closest_point_A7star(const double* x, double* y)
{
    long   k[8];  double z[8];  int b[8];
    double bd[9] = { 0.0 };  int bc[9] = { 0 };
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int cnt = 0, m = 0, fbest = 0;
    const double INV = 1.0 / 8.0;

	/* round, compute residuals, and sum */
    RND_ACC(0, x, k, z, sumk); RND_ACC(1, x, k, z, sumk);
    RND_ACC(2, x, k, z, sumk); RND_ACC(3, x, k, z, sumk);
    RND_ACC(4, x, k, z, sumk); RND_ACC(5, x, k, z, sumk);
    RND_ACC(6, x, k, z, sumk); RND_ACC(7, x, k, z, sumk);
    alpha = -(double)sumk;

	/* bucket residuals into b[j] = clamp(ceil(N*(1/2 - z_j)), 1, N) and scatter (1 - 2 z_j) and a unit count */
    MCSQ_BUCKET(0, 8, z, b, bd, bc); MCSQ_BUCKET(1, 8, z, b, bd, bc);
    MCSQ_BUCKET(2, 8, z, b, bd, bc); MCSQ_BUCKET(3, 8, z, b, bd, bc);
    MCSQ_BUCKET(4, 8, z, b, bd, bc); MCSQ_BUCKET(5, 8, z, b, bd, bc);
    MCSQ_BUCKET(6, 8, z, b, bd, bc); MCSQ_BUCKET(7, 8, z, b, bd, bc);

	/* find prefix that minimises h = beta_var - (alpha - cnt)^2 / N: */
    best = -(alpha * alpha) * INV;
    MCSQ_SCAN(1, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(2, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(3, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(4, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(5, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(6, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(7, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(8, bd, bc, INV, alpha, acc, cnt, best, m, fbest);

	/* flip the chosen coordinates' k upward by 1 */
    MCSQ_RECON(0, b, m, k); MCSQ_RECON(1, b, m, k); MCSQ_RECON(2, b, m, k);
    MCSQ_RECON(3, b, m, k); MCSQ_RECON(4, b, m, k); MCSQ_RECON(5, b, m, k);
    MCSQ_RECON(6, b, m, k); MCSQ_RECON(7, b, m, k);

    /* compute the output lattice point */
    mean = (double)(sumk + fbest) * INV;
    y[0] = (double)k[0] - mean; y[1] = (double)k[1] - mean; y[2] = (double)k[2] - mean;
    y[3] = (double)k[3] - mean; y[4] = (double)k[4] - mean; y[5] = (double)k[5] - mean;
    y[6] = (double)k[6] - mean; y[7] = (double)k[7] - mean;
}

/*!
 *  \brief Closest point of A8* (modified MCSQ algorithm).
 *
 *   A8*:   N = 9   single-pass bucketing beats 25-comparator network
 *
 *  \param[in]  x  9-vector, zero-sum.
 *  \param[out] y  9-vector, closest A8* point (zero-sum).
 */
static void closest_point_A8star(const double* x, double* y)
{
    long   k[9];  double z[9];  int b[9];
    double bd[10] = { 0.0 };  int bc[10] = { 0 };
    double alpha, acc = 0.0, best, mean;
    long   sumk = 0;  int cnt = 0, m = 0, fbest = 0;
    const double INV = 1.0 / 9.0;

    RND_ACC(0, x, k, z, sumk); RND_ACC(1, x, k, z, sumk);
    RND_ACC(2, x, k, z, sumk); RND_ACC(3, x, k, z, sumk);
    RND_ACC(4, x, k, z, sumk); RND_ACC(5, x, k, z, sumk);
    RND_ACC(6, x, k, z, sumk); RND_ACC(7, x, k, z, sumk);
    RND_ACC(8, x, k, z, sumk);
    alpha = -(double)sumk;

    MCSQ_BUCKET(0, 9, z, b, bd, bc); MCSQ_BUCKET(1, 9, z, b, bd, bc);
    MCSQ_BUCKET(2, 9, z, b, bd, bc); MCSQ_BUCKET(3, 9, z, b, bd, bc);
    MCSQ_BUCKET(4, 9, z, b, bd, bc); MCSQ_BUCKET(5, 9, z, b, bd, bc);
    MCSQ_BUCKET(6, 9, z, b, bd, bc); MCSQ_BUCKET(7, 9, z, b, bd, bc);
    MCSQ_BUCKET(8, 9, z, b, bd, bc);

    best = -(alpha * alpha) * INV;
    MCSQ_SCAN(1, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(2, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(3, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(4, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(5, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(6, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(7, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(8, bd, bc, INV, alpha, acc, cnt, best, m, fbest);
    MCSQ_SCAN(9, bd, bc, INV, alpha, acc, cnt, best, m, fbest);

    MCSQ_RECON(0, b, m, k); MCSQ_RECON(1, b, m, k); MCSQ_RECON(2, b, m, k);
    MCSQ_RECON(3, b, m, k); MCSQ_RECON(4, b, m, k); MCSQ_RECON(5, b, m, k);
    MCSQ_RECON(6, b, m, k); MCSQ_RECON(7, b, m, k); MCSQ_RECON(8, b, m, k);

    mean = (double)(sumk + fbest) * INV;
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
 *  \brief Closest point of An* (fastest methods for n=2..8).
 */
int closest_point_Anstar_faster(const double* x, double* y, int n)
{
    /* check parameters: */
	if (n < 2 || n > 8) return LQ_INVARG;

    /* call the appropriate size-optimized function: */
	closest_point_Anstar[n](x, y);

    return LQ_SUCCESS;
}

/* fanstar_faster.c -- end of file */
