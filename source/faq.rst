.. _faq:

==============
FAQ
==============

--------------
General
--------------

**Where can I download the code and what version of the code should I use?**
   
   We strongly recommend to always use the latest `release <https://github.com/Nek5000/Nek5000/releases>`__  available.
   We do **not** recommend using the master branch on `GitHub <https://github.com/Nek5000/Nek5000>`__
   in a production environment!

**How can I properly reference Nek5000?**

   If you use our software, please cite the following:

::

  NEK5000 Version X.Y. Release date. Argonne National Laboratory, Illinois. 
  Available: https://nek5000.mcs.anl.gov.

**What is the license for Nek5000?**

   Nek5000 is licensed under BSD.  
   For more information see the ``LICENSE`` file included with the distribution in the root level directory.

**Where can I get help?**

   If you have a question, first check the `Google group <https://groups.google.com/forum/#!forum/nek5000>`__ to see if your question is already answered somewhere. 
   The google group also serves as a primary support channel with the Nek5000 user community. 
   Please subscribe to the google group by clicking the button "Apply to join group" right above the "Welcome to the new Nek5000 User Group!" sign.
   
   The past `mailing list archive <https://lists.mcs.anl.gov/pipermail/nek5000-users>`__ can also be checked for potential answers.

**How can I report a bug / feature request?**

  Nek5000 is hosted on GitHub and all bugs are reported and tracked through the `Issues <https://github.com/Nek5000/Nek5000/issues>`__ feature on GitHub. 
  However, GitHub Issues should not be used for common troubleshooting purposes. 
  If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User’s Group mailing list. 
  If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

**How can I contribute to the Nek5000 project?**

  Our project is hosted on `GitHub <https://github.com/Nek5000>`__. Here are the most important things you need to know:
  
  - follow the usual “fork-and-pull” Git workflow
  - all development happens on the master branch
  - anything in master is always deployable
  - upcoming releases get their own tags out of master
  
  If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.

----------------------------------
Installing, Compiling, and Running
----------------------------------

**Which platforms are supported by Nek5000?**

   All posix compliant operations system. 
   Compiling and running on Microsoft Windows 10 has been demonstrated using CMake and MinGW, although this is **not** supported.

**Can Nek5000 run on GPUs?**

   Currently there is no support for GPUs. We are actively working on a high-performance GPU version but this is
   a topic of active research. 

**Which compilers work?**

   Currently Intel, PGI, IBM and GNU compilers are supported.

**How much memory is required?**

   The memory footprint of a run depends on many factors and is printed to
   screen whenever Nek5000 exits. What follows is a first rough guess:

.. code-block:: none

   lx1*ly1*lz1*lelt * 3000byte + lelg * 12byte + MPI + optional libraries (e.g. CVODE)
..

   where ``lelt`` (the maximum number of local elements) is computed as lelg/lpmin.
   The memory allocated by MPI will depend heavily on the total number of ranks and the considered MPI implementation. 
   For large rank counts (say > 100'000) it's easily 50-100MB.

   Note, the output of GNU`s SIZE utility is inaccurate as it does not take into account the dynamic memory alloation of MPI, gslib, CVODE, etc. 

**Why does the compiler issue relocation errors?**

   ``relocation truncated to fit: R_X86_64 against symbol 'foo_' defined in COMMON section in obj/bar.o``

   This happens when the resultant executable requires more than 2GB of static data.  
   The best way to avoid this is to increase the minimum number of MPI ranks, ``lpmin`` in ``SIZE``.  
   Alternatively, you can add ``-mcmodel=medium`` to the ``FFLAGS/CFLAGS`` variable in ``makenek``.

**Why does the code use static memory allocation?**

   It is easy to program, produces fast code and is less error prone. The drawback of recompilation is not relevant 
   as full rebuild takes less than 1 minute. 

**How do I launch a parallel run?**
  
  Assuming your machine has MPI installed correctly, parallel runs can be launced with the included ``nekmpi`` and ``nekbmpi`` scripts in ``Nek5000/bin``. 
  If you are running on a machine with a queuing system, consult your sysadmin how to submit a job.

**How do I run the examples?**

  The examples are included by default in the release tarball (see example directory). There is nothing special you need
  to do as they are ready to run.  

---------------------------
Pre-Processing and Numerics
---------------------------

**How can I generate a mesh for use with Nek5000?**

   Please see quickstart section on :ref:`qstart_meshing`.

**What element types are supported?**

   Conformal curved quadrilateral/hexahedral elements.

**How do I import/convert a mesh to Nek5000?**

   We currently support conversion from the exodusII with the ``exo2nek`` converter. This enables the import from popular mesh generators like ANSYS ICEM and CUBIT.

**Why is it important to non-dimensionalize my case?**

  Nek5000 can be run with dimensions, but we STRONGLY recommend that the case has been non-dimensionalized properly.
  An advantage of the nondimensional form is that physical simulation times, tolerances, etc. tend to
  be easy to set based on prior experience with other simulations.

**How do I choose solver tolerances?**

  Depends on how accurate you need your simulation to be.  
  Typical values (for engineering type of problems) are :math:`10^{-7}` for velocity and scalars.
  In Pn/Pn-2 the pressure tolerance is equal to desired error in divergence. This is in contrast to Pn/Pn where the divergence
  error is mainly a function of spatial resolution and a tolerance of :math:`10^{-4}` is typically good enough.   
  Note the tolerances are related to the residual in the linear solve and do not represent the accuracy of the solution. 

**What formulation Pn/Pn or Pn/Pn-2 should I use?**

   There is no simple answer but we typically recommend to use the Pn/Pn formulation altough not all features are supported (at least for now). 

**What polynomial order should I use?**

  The code supports a large range of polynomial orders, e.g. :math:`N=1` through :math:`N=32`.
  You can effectively realize the same number of grid points
  by using relatively few high-order elements or more low-order elements.
  For example, a 3D grid with resolution of 64x64x64 could be implemented
  as a 16x16x16 array of elements of order :math:`N=3` or as a
  8x8x8 array of elements of order :math:`N=7`.  In Nek5000, the 
  latter is preferred. The solution will be more accurate and the code
  is optimized for this range of :math:`N`.

  The sweet spot is typically :math:`N=7` (``lx1=8``). 

.. Unless you have a very good reason to change it do not deviate from this best practice. 

.. Note, do never use :math:`N<5` as this results in a very poor performance. 

**How do I specify/change the polynomial order?**

   Change ``lx1`` in the SIZE file. Note, the polynomial order is :math:`N=lx1-1`. 

**How do I specify/change the solver runtime parameters?**

   See the section on the :ref:`case_files_par` file.

**Why is ``userbc`` only called for certain element faces?**

   ``userbc`` is ONLY called for element boundary conditions specified with a lower-case letter, e.g. 'v', 't', or 'o' but NOT 'W', 'E', or 'O'.  Note that this implies it is not necesarily called on all MPI ranks.

**How do I solve for a scalar?**

   Nek5000 supports solving up to 99 additional scalars.  
   To solve an additional scalar equation, increase ``ldimt`` in the ``SIZE`` file to accomodate the additional scalar and specify the appropriate parameter in the :ref:`case_files_par` file. See ``shear4`` example for more details. 

**Are there any books/papers that describe the numerics of Nek5000?**

  There are a number of descriptions of the various numerical methods used in Nek5000 available.
  Probably the best starting point is the book *High-Order Methods for Incompressible Fluid Flow* by Deville et al. (2002).
  There are also other, perhaps shorter, expositions of the material.
  Two papers that we found particularly useful (there are of course many more) are:

  - Fischer. An Overlapping Schwarz Method for Spectral Element Solution of the Incompressible Navier–Stokes Equations. *J. Comput. Phys.* 133, 84–101 (1997)

  - Fischer et al. Simulation of high-Reynolds number vascular flows. *Comput. Methods Appl. Mech. Engrg.* 196 (2007) 3049–3060

  and also the `lecture notes <http://www.mcs.anl.gov/~fischer/kth/kth_crs_2016s.pdf>`_ by Paul Fischer (given at KTH in 2016).

**Why can I see sometimes the imprint of elements in the solution?**

  Nek5000 is based on the spectral-element method, which relies on an expansion of the solution in terms of element-local basis functions.
  These basis functions are the Lagrange interpolants to the Legendre polynomials of a specific order.
  If using PnPn-2, the velocity is on the Gauss-Lobatto-Legendre mesh (i.e. including the boundary points), and the pressure is on the Gauss-Legendre mesh (without boundary points).
  These functions are defined within each element, and the continuity between elements is C0, i.e. only the function value is the same.
  The ansatz functions are polynomials, so you can differentiate them inside each element; however, derivatives are not continuous over element boundaries (even though this difference reduces spectrally fast). 
  Note that for the PnPn-2 method, the pressure is non-continuous.

  This means that when visualising e.g. derivatives, one might see discontinuities, which then appear as imprints of the elements.
  This is due to the mentioned properties of the discretisation, and as such not a sign of a wrong solution.
  With increasing resolution (either p or h-type) these jumps will most certainly get smaller.

---------------------------
Physical Models
---------------------------

**What turbulence models are available in Nek5000?**

   For LES we provide an explicit filtering approach or a relaxation term model. 
   RANS turbulence models (k-ω, k-ω SST, etc.) are not an integral part of the code but available through examples.

-------------------
Computational Speed
-------------------

**Are there any compiler specific flags I should use?**

  Compile with vector instructions like AVX, AVX2 using FFLAGS and CFLAGS 
  in makenek.   

**How many elements should I have per process?**

  The upper limit is given by the available memory. The lower limit is (technically) 1 but you may want to have more
  elements (work) to get a reasonable (whatever that means for you) parallel efficiency. 
  On most machines you need more than 10 elements per MPI rank to get a parallel efficiency of 0.5 (assuming N=7).  

**Should I use residual projection?**

  Typically projection is used for pressure but not velocity, however
  this is highly case specific and a simple experiment will show if it pays off or not.  
  Projection will speed up the solution to a scalar, but takes time to compute itself.
  A scalar solve requiring ~40 iterations or greater is a good candidate for use.

**What other things can I do to get best performance?**

  - Design your mesh for a polynomial order N=7
  - Tune your solver tolerances
  - Increase time step size by switching to 2nd order BDF and OIFS extrapolation (target Courant number 2-5)
  - Use AMG instead of XXT as coarse grid solver
  - Avoid unnecessary time consuming operations in ``usrchk/userbc``
  - Use binary input files e.g. ``.re2`` and ``.ma2`` to minimize solver initialization time
  - Use a high performance MXM implementation for your platform (see ``makenek`` options)

---------------------------
Troubleshooting
---------------------------

**My simulation diverges. What should I do?**

  There are many potential root causes but here are some things you can experiment with:

  * lower the time step (in particular during initial transients) 
  * reduce time integration order (e.g. use 2 instead of 3)
  * increase spatial resolution
  * provide a better initial condition
  * check that your boundary conditions are meaningful and correctly implemented 
  * visualize the solution and look for anomalies

---------------
Post-Processing
---------------

**What options are available**

   * For data analysis you use Nek5000's internal machinery through the usr file
   * Solution files can be read by VisIt and Paraview (for more information see :ref:`qstart_vis`)
   * Various user contributions in `NekBazaar <https://github.com/Nek5000/NekBazaar/>`_ 

**The local coordinate axes of my elements are not aligned with the global coordinate system, is this normal?**

   Yes, there is no guarantee that the elements are generated with any particular orientation (except if you use genbox).

**Where are my solution files?**

   By default Nek5000 outputs solution files in binary ``<casename>.f%05d``.  

**I have calculated additional fields from my solution, how do I visualize them?**

   Using the ``.par`` file, define an additional scalar and include the ``solver=none`` option.
   For example:

.. code-block:: none

   [SCALAR01] # lambda2 vortex criterion
   solver = none

..

   Then store the calculated field in ``t(1,1,1,1,iscal+1)`` where ``iscal`` is your passive scalar index (in this example 1).
   The scalar will then be output by default with the solution files.

**How do I obtain values of variables at a specific point?**


  The simplest way is through the use of history points. See the section on the :ref:`features_his` file.
  You can also use the spectral interpolation tool (see examples for more details).
