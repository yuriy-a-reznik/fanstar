# fanstar
Faster algorithms for the An* lattices.

This project implements closest-point algorithms for the An* lattices based on the following papers:

[CS]      J. H. Conway and N. J. A. Sloane, "Fast quantizing and decoding algorithms for lattice quantizers
          and codes," IEEE Trans. Inf. Theory, vol. 28, no. 2, pp. 227–232, Mar. 1982.
          
[VC]      I. Vaughan and L. Clarkson, "An algorithm to compute a nearest point in the lattice An*," in Applied Algebra, 
          Algebraic Algorithms and Error-Correcting Codes (Lecture Notes in Computer Science), vol. 1719, 
          Springer, 1999, pp. 104–120.
          
[MCQ]     R. G. McKilliam, I. V. L. Clarkson, and B. G. Quinn, "An algorithm to compute the nearest point in 
          the lattice An*," IEEE Trans. Inf. Theory, vol. 54, no. 9, pp. 4378–4381, Sep. 2008.
          
[MCSQ]    R. G. McKilliam, I. V. L. Clarkson, W. D. Smith, and B. G. Quinn, "A linear-time nearest
          point algorithm for the lattice An*," in Proc. Int. Symp. Information Theory and its Applications
          (ISITA), Dec. 2008, pp. 1–5.
          
[FANSTAR] Y. Reznik, "Faster closest-point algorithms for the An* lattices" -- in preparation, 2026

The [FANSTAR] algorithm is ours, rolling in several improvements over [MCQ] and [MCSQ] methods.  

For reference, the project also implements the closest-point algorithms for the Zn, An, Dn, Dn* and E8 lattices based on the [CS] reference. The test program compares the speed and accuracy of all these methods for n = 2..8.
