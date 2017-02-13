Nek5000
=======

Nek5000 is designed to simulate laminar, transitional, and turbulent incompressible or low
Mach-number flows with heat transfer and species transport. It is also suitable for incompressible
magnetohydrodynamics (MHD). Nek5000 is written in Fortran 77 and C. It uses MPI for message passing
(but can be compiled without MPI for serial applications) and some LAPACK routines for eigenvalue
computations (depending on the particular solver employed).  Nek5000 output formats can be read by
either ``postx`` (distributed with Nek5000) or the parallel visualization package VisIt developed
by Hank Childs and colleagues at LLNL/LBNL.  VisIt is mandatory for large problems (e.g., more than
100,000 spectral elements).

.. toctree::
   :maxdepth: 2

   intro
   quickstart
   user_files
   geometry
   large_scale
   routines
   appendix
   bibliography


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
