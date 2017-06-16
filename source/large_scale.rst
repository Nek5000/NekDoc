=============================================
Performing Large-Scale Simulations in Nek5000
=============================================

-----------------------
Large-Scale Simulations
-----------------------

The largest simulations performed so far with Nek5000 are in the range of :math:`5\times 10^6`.
Performance aspects to keep in mind

- design your SEM-mesh for a polynomial order N=7 or N=9 (lx1 in SIZE)
- turn on dealiasing only if needed and try to minimize the polynomial order used for the fine grid (lxd in SIZE)
- ensure that you have at least 50 elements per MPI process (the more the better)
- explore the pressure solver performance using AMG solver (sweet spot depends on number of processors and elements)
- make sure usrchk() does not contain time consuming operations getting called in the time loop
- enable tuned MxM implementation for your platform (see makenek options)
- turn on residual projection scheme (see .rea file parameters)
- tune your pressure/velocity tolerances (e.g. use 5e-5 for pressure and 1e-8 for velocity solver) having in mind your overall accuracy
- try to maximize the timestep, if needed turn on OIFS scheme with a target Courant number of around 2.0
- use the parallel I/O (MPI-IO or custom kernel) to write checkpoints to disk
- understand where you spend most of the time - turn on solver and MPI timings to monitor solver performance (see makenek options)

In such context it is recommendable to save the .rea file as a binary re2.

Also for input/output it may be necessary to use MPI I/O. In this case the code has to be compiled with MPI I/O, i.e. the line PPLIST="MPIIO" in the makefile should not be commented. For output we may use iofiles, the parameter 65 in the .rea file specifies the number of directories and separate files that have to be created as specified by the user.

For large scale simulations the AMG solver is a better, faster choice for solving the Poisson problem (the default solver is XXT).

..........
AMG solver
..........

The code should be compiled once with the settings ``AMG=true``, ``AMG_DUMP=true``. In the tools folder of Nek5000 we can find the AMG solver, a Matlab version for the moment which is subject to further integration in the main code. The user should run the script run which will read the AMG dump files and create new ones. The new files are now to be used in the code and with ``AMG_DUMP`` commented out the user should recompile and run his Nek5000 version.

The AMG solver is a 3 stage process.

The first step will generate the files needed for the matlab code. Next matlab must run the setup code and generate a set of .dat files. Then nek can run with the ``.dat`` files and use the AMG pressure solver.

AMG dump stage

    Make sure ``IFAMG`` and ``IFAMG_DUMP`` in makenek are uncommented and set to true

    Run ``makenek clean`` then ``makenek <casename>``

    Run Nek (this will produce a set of ``*.dat`` files)

MATLAB AMG stage

    Move the ``amgdmp_*.dat`` files to ``nek5_svn/trunk/tools/amg_matlab``:

    ``mv amgdmp*.dat ../../trunk/tools/amg_matlab``

    ``cd ../../trunk/tools/amg_matlab``

    Run the script: ``tools/amg_matlab/run`` (this may take several hours and will produce set of files)

AMG run stage

    Move all ``*.dat`` files produced back to your case directory:

    ``mv *.dat /path/to/case/dir``

    Comment ``IFAMG_DUMP`` in ``makenek`` (``IFAMG`` should still be set to ``TRUE``)

    Run ``makenek clean``, then run ``makenek <casename>``

    Run Nek (the generated AMG files will be read during the AMG setup phase)

Notes on improving AMG results:

    To help speed up the matlab process, try running the 1st stage, the AMG dump stage, with ``lx1=3`` in the SIZE file. Using a lower lx1 number will create a sparser matrix and thus a speedier matlab resolution. lx1 can be increased when ready to run the 2nd stage, the AMG run stage, after the .dat files are produced.

    To increase accuracy in the AMG results, try tightening the tolerances in the run script, in ``trunk/tools/amg_matlab``. Specifically, the first tolerance (default set to 0.5). Lowering this (say, to 0.1), will increase the time the matlab code stage takes, but the result will be a faster convergence in the pressure solves of the AMG run stage.

...................
Size related issues
...................

Large simulations require a high number of nodes and thus a high number of variables. So far Nek5000 supports by default 4 byte integers, however the node numbering for big meshes may exceed this size and 8 byte integers may be needed. How is this done?

If more than 9 passive scalars are needed

Exiting Nek5000 while a batch job in being executed should be done not using ``"qdel"`` but ``echo -1 > ioinfo``.

.......
MAKENEK
.......

The shell script makenek is designed to assist the compilation process of NEK5000. The script will create a makefile based on the user settings section in makenek. The GNU gmake utility is used to build NEK5000.
Available configurations options:

.. _tab:bdms:

.. table:: Compiler options

   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | name           | values     | default       | description                                                                              |
   +================+============+===============+==========================================================================================+
   | PPLIST         | string     |               | list of pre-processor symbols (BG,MOAB,BLAS_MXM, MPIIO)                                  |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | IFMPI          | true,false | true          | use MPI (needed for a multiprocessor computation)                                        |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | IFAMG_DUMP     | true,false | false         | dump AMG pre-processing files                                                            |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | IFAMG          | true,false | false         | use AMG as coarse grid solver for pressure preconditioner else XXT                       |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | F77            | string     | mandatory     | Fortran compiler (e.g. MPI: mpif77)                                                      |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | CC             | string     | mandatory     | C compiler (e.g. MPI: mpicc)                                                             |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | G              | string     | optional      | optional compilation flags                                                               |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | OPT_FLAGS_STD7 | string     | optional      | optimization flags for L1,L2,L3                                                          |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | OPT_FLAGS_MAG  | string     | optional      | optimization flags for L4 (highest opt level)                                            |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | SOURCE_ROOT    | string     | mandatory     | path of nek5000 source                                                                   |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | USR            | string     | optional      | object list of additional files to compile (make intructions (makefile_usr.inc required) |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | USR_LFLAGS     | string     | optional      | optional linking flags                                                                   |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | MOAB_DIR       | string     | NEK with MOAB | Path to MOAB directories                                                                 |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | IFVISIT        | true,false | false         | Toggles Visit in situ. See Visit_in_situ for details                                     |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | VISIT_INSTALL  | string     | VISIT in situ | Path to VISIT install path. See Visit_in_situ for details.                               |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+
   | VISIT_STOP     | true,false | false         | When running VISIT in situ, simulation stops after step 1 to connect VISIT.              |
   +----------------+------------+---------------+------------------------------------------------------------------------------------------+

...............
Binary geometry
...............

Reatore2
Jump to: navigation, search

The NEK5000 tool, reatore2 allows users to split an ASCII .rea file to an ASCII .rea and a binary .re2 file. The .re2 file contains the mesh and boundary condition data that is normally written in ASCII in the .rea file. For large simulations, this information can be substantial, so storing it in binary lowers the memory footprint for the simulation.
Running reatore2

Be sure that your nekton tools are up-to-date and compiled.
At the command prompt type: reatore2

NOTE-If the executables for the tools were not placed in the bin directory(default),
include the path to the reatore2 executable

    User is prompted for name of .rea file

    -Enter the name to the .rea file, excluding the .rea extenstion

    User is prompted for the new files name

    -Enter the name for your new files

----------------------
Parallelism in Nek5000
----------------------

The parallelism of Nek5000 is accomplished via domain decomposition methods and a suitable gather-scatter code. All this is implemented in such way that the user does not have to be concerned with the parallelism and only focus on the actual solvers while keepin in mind a few simple rules and routines that switch from local to global and back.

- Locally, the SEM is structured.
- Globally, the SEM is unstructured.
- Vectorization and serial performance derive from the structured aspects of the computation.
- Parallelism and geometric flexibility derive from the unstructured, element-by-element, operator evaluation.
- Elements, or groups of elements are distributed across processors, but an element is never subdivided.

For the most part, the global element numbering is not relevant since Nek5000 assigns it randomly but following certain rules.

There are two types of array sizes, starting with lx1, lelv, etc. which give an upper bound of the arrays. And nx1, nelv, etc. which give the actual number of elements/grid points per processors. For the example in :numref:`fig:procsplit` we have

- on proc 0, ``nelt=2``  (nelt = no elements in temperature domain)
- on proc 1, ``nelt=3``  (nelv = no elements in fluid domain, usually = nelt)

.. _fig:procsplit:

.. figure:: figs/serial_parallel.png
    :align: center
    :figclass: align-center
    :alt: element-splitting

    A simple SEM row of elements and a potential splitting

Arrays ``lglel`` that distinguish which processor has which elements,

- on proc 0, ``nelt=2, lglel=(2,5)``, local element ``1->2`` and ``2->5``
- on proc 1, ``nelt=3, lglel=(1,3,4)``, local element ``1->1``, ``2->3`` and ``4->3``


Now for global to local we have two common arrays (scaling as ``nelgt``, but only two such arrays)

- ``gllel=(1,1,2,3,2)``, assigns a global element to its local correspondent, i.e. global element ``1->1``, ``2->1`` and ``3->2`` etc.
- ``gllnid=(1,0,1,1,0)``, assigns a global element to its processor, i.e. ``1->1``, ``2->0`` and ``3->1`` etc.

All data contiguously packed (and quad-aligned) ``real  u(lx1,ly1,lz1,lelt)`` indicates that ``u`` is a collection of elements, ``e=1,\(\ldots\),Nelt =< lelt``, each of size :math:`(N+1)d, d=2 or 3`.

**Example: Summation**

Serial version

.. code-block:: fortran

   s = 0
   do e=1,nelv
   do iz=1,nz1
   do iy=1,ny1
   do ix=1,nx1
   s=s+u(ix,iy,iz,e)
   end do,...,end do

Second approach, serial version (works in parallel in Nek)

.. code-block:: fortran

   n=nx1*ny1*nz1*nelv
   s=0
   do i=1,n
   s=s+u(i,1,1,1)
   end do

Nek Parallel Version

.. code-block:: fortran

   s=glsum(s,n)

If you want a local max ``s=vlmax(u,n)``, or a global max ``s=glmax(u,n)``.

