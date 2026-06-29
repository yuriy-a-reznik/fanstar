/*!
 *  \file       fanstar.h
 *  \brief      FANSTAR:  -- internal header file for the project.
 *
 *  \author     Yuriy A. Reznik <yreznik@mit.edu>
 *  \version    1.05
 *  \date       May 25, 2026
 *
 *  \copyright  Copyright (C) 2026 Yuriy A. Reznik
 *  \license    MIT License: https://opensource.org/licenses/MIT
 * 
 *  This project implements and compares several closest-point algorithms for the An* lattices based on the following papers:
 * 
 *   [CS]      J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers
 *             and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227–232, Mar. 1982.
 * 
 *   [VC]      I. Vaughan and L. Clarkson, "An algorithm to compute a nearest point in the lattice A∗n," in Applied Algebra, 
 *             Algebraic Algorithms and Error-Correcting Codes (Lecture Notes in Computer Science), vol. 1719, 
 *             M. Fossorier, H. Imai, S. Lin, and A. Poli, Eds. Berlin, Germany: Springer-Verlag, 1999, pp. 104–120.
 * 
 *   [MCQ]     R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in 
 *             the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378–4381, Sep. 2008.
 * 
 *   [MCSQ]    R. G. McKilliam, I. V. L. Clarkson, W. D. Smith, and B. G. Quinn, "A linear-time nearest
 *             point algorithm for the lattice An*," in Proc. Int. Symp. Information Theory and its Applications
 *             (ISITA), Auckland, New Zealand, Dec. 2008, pp. 1–5.
 * 
 *   [FANSTAR] Y. Reznik, "Faster closest point algorithms for the An* lattices" -- in preparation, 2026
 *             This is our proposed faster method. 
 *
 *  For reference, it also implements the closest-point algorithms for the Zn, An, Dn, Dn*, and E8 lattices by following [CS].
 *  The test program compares the speed and accuracy of all these methods for n = 2..8.
 */

#ifndef _FANSTAR_H
#define _FANSTAR_H

#ifdef __cplusplus
extern "C" {
#endif

/*! Maximum supported dimensions: */
#define LQ_MAX_N    8                   /*!< Maximum lattice dimension n                */
#define LQ_MAX_M    (LQ_MAX_N+1)        /*!< Maximum lattice ambient dimension m        */

/*! Error codes: */
enum {
    LQ_SUCCESS = 0,         /*!< Success                                                */
    LQ_INVARG = 1,          /*!< Invalid argument(s)                                    */
    LQ_NOTSUP = 2,          /*!< Lattice type/dimension combination is not supported    */
};

/*! 
 *  \brief Closest point algorithms for An* lattices.
 *
 *  \param[in]  x        Input vector (in lattice augmented dimension corresponding to n).
 *  \param[out] y        Output vector for the closest point.
 *  \param[in]  n        Lattice dimension parameter.
 */
extern int closest_point_Anstar_CS(const double* x, double* y, int n);          /* Conway & Sloan's closest point algorithm for An* lattices */
extern int closest_point_Anstar_VC(const double* x, double* y, int n);          /* Vaughan & Clarkson's closest point algorithm for An* lattices */
extern int closest_point_Anstar_MCQ(const double* x, double* y, int n);         /* McKilliam, Clarkson & Quinn's closest point algorithm for An* lattices */
extern int closest_point_Anstar_MCSQ(const double* x, double* y, int n);        /* McKilliam, Clarkson, Smith & Quinn's closest point algorithm for An* lattices */
extern int closest_point_Anstar_MCQ_opt(const double* x, double* y, int n);     /* Optimized McKilliam, Clarkson, & Quinn's closest point algorithm for An* lattices */
extern int closest_point_Anstar_MCSQ_opt(const double* x, double* y, int n);    /* Optimized McKilliam, Clarkson, Smith & Quinn's closest point algorithm for An* lattices */
extern int closest_point_Anstar_faster(const double* x, double* y, int n);     /* Proposed method, rolling in all optimizations */

/* Closest point algorithms for some other lattices */
extern int closest_point_Zn_CS(const double* x, double* y, int n);              /* Conway & Sloan's closest point algorithm for Zn lattices  */
extern int closest_point_An_CS(const double* x, double* y, int n);              /* Conway & Sloan's closest point algorithm for An lattices  */
extern int closest_point_Dn_CS(const double* x, double* y, int n);              /* Conway & Sloan's closest point algorithm for Dn lattices  */
extern int closest_point_Dnstar_CS(const double* x, double* y, int n);          /* Conway & Sloan's closest point algorithm for Dn* lattices */
extern int closest_point_E8_CS(const double* x, double* y, int n);              /* Conway & Sloan's closest point algorithm for E8 lattice   */

#ifdef __cplusplus
}
#endif

#endif /* _FANSTAR_H */
