.. todo::
    Add hyperlinks or bibliography refs to VisIt
    Proofread "design principles"

============
Introduction
============

----------------------
Computational Approach
----------------------

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

The code Nek5000 is based on the following design principles:

- Accessible both to beginners and experts alike.
- Accessible interface via Fortran to include user-defined modules.
- The code intrinsics can be accessed and modified via the user files for more experienced developers
- portability.
- Minimal use of external libraries to assure fast compile times.
- Fast matrix free operator evaluation with minimal storage.
- Matrix operations are implemented in assembler code :math:`M \times M` routines to speed up computations.
- The parallelism is "under the hood", demanding from the user only care in handling local versus global operations and array sizes.
- By testing at the beginning of each run which one of the three readily implemented parallel-algorithms behaves optimally, the parallelism of Nek5000 can be automatically tuned to each machine.
- Direct access to parameters at runtime.
- Geometry and boundary conditions are exposed to the user via the .rea file.
- Handling complex geometries that can be imported from external codes .

--------------------------------------
Incompressible Navier-Stokes Equations
--------------------------------------

-----------------------------
Non-Dimensional Navier-Stokes
-----------------------------

---------------
Energy Equation
---------------

------------------------------------------------
Non-Dimensional Energy / Passive Scalar Equation
------------------------------------------------

---------------
Passive Scalars
---------------

---------------
Unsteady Stokes
---------------

-------------
Steady Stokes
-------------

--------------------
Linearized Equations
--------------------

-----------------
Steady Conduction
-----------------

----------------------
Low-Mach Navier-Stokes
----------------------

----------------------------
Incompressible MHD Equations
----------------------------

----------------------------------
Adpative Lagrangian-Eulerian (ALE)
----------------------------------
