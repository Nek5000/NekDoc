.. _intro_compute_approach:

Computational Approach
======================

The spatial discretization is based on the spectral element method (SEM) [Patera1984]_, which is a
high-order weighted residual technique similar to the finite element method.   In the SEM, the
solution and data are represented in terms of \(N\)th-order tensor-product polynomials within each
of :math:`E` deformable hexahedral (brick) elements. Typical discretizations involve
:math:`E`\=100--10,000 elements of order :math:`N`\=8--16 (corresponding to 512--4096 points per
element).  Vectorization and cache efficiency derive from the local lexicographical ordering within
each macro-element and from the fact that the action of discrete operators, which nominally have
:math:`O(EN^6)` nonzeros, can be evaluated in only :math:`O(EN^4)` work and :math:`O(EN^3)` storage
through the use of tensor-product-sum factorization (sao80).   The SEM exhibits very little
numerical dispersion and dissipation, which can be important, for example, in stability
calculations, for long time integrations, and for high Reynolds number flows. We refer to
[Denville2002]_ for more details.
