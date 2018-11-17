==========
Appendices
==========

-----------------
Build Options
-----------------

The shell script ``makenek`` is designed to assist the compilation process of Nek5000. The script will create a ``makefile`` based on the user settings section in ``makenek``. The GNU gmake utility is used to build Nek5000.
Available configurations options:

.. _tab:bdms:

.. csv-table:: Compiler options
   :header: name,values,default,description
   :widths: 12,7,12,20

   PPLIST, string, , "list of pre-processor symbols (CVODE, ...)"                                     
   MPI, "1, 0", 1, use MPI (needed for a multiprocessor computation)                                           

   FC, string, optional, Fortran compiler (mpif77)                                                         
   CC, string, optional, C compiler (mpicc)                                                               
   FCLAGS, string, optional, optional Fortan compilation flags        
   CCLAGS, string, optional, optional C compilation flags                                                                  
   SOURCE_ROOT, string, optional, path of Nek5000 source                                                                      
   USR, string, optional, object list of additional files to compile make intructions (``makefile_usr.inc`` required) 
   USR_LFLAGS, string, optional, optional linking flags                                                                      
   PROFILING, "1, 0", 1, enable internal timers for performance statistics                                       
   VISIT, "1, 0", 0, Toggles Visit in situ. See Visit_in_situ for details                                        
   VISIT_INSTALL, string, VISIT in situ, Path to VISIT install path. See Visit_in_situ for details.                                 
   VISIT_STOP, "true, false", false, "When running VISIT in situ, simulation stops after step 1 to connect VISIT."                 


The ``PPLIST`` field can be used to activate several features at compilation time. 
A list of possible options is below:

.. _tab:PPLIST

.. csv-table:: PPLIST options
   :header: Symbol, Description

   NOMPIIO, deactivate MPI-IO support
   BGQ, use Blue Gene Q optimized mxm
   XSMM, use libxsmm for mxm
   CVODE, compile with CVODE support for scalars
   VENDOR_BLAS, use VENDOR BLAS/LAPACK
   EXTBAR, add underscore to exit call (for BGQ)
   NEKNEK, activate overlapping mesh solver (experimental)
   CMTNEK, activate discontinuous Galerkin compressible-flow solver (experimental)

In addition to these preprocessor items, the user can add compilation and linking flags. 
``FFLAGS`` allows the user to add Fortran compilation flags while ``CCFAGS`` allows the user to 
add C compilation flags. 
These will be compiler dependent and the user is encouraged to consult the manual of the compiler if specific options are needed/desired. 
A commonly used flag is ``-mcmodel`` which allows for arrays of size larger than 2GB. 
This option  tells the compiler to use a specific memory model to generate code and store data. 
It can affect code size and performance. 
If your program has global and static data with a total size smaller than 2GB, ``-mcmodel=small`` is sufficient. 
Global and static data larger than 2GB requires ``-mcmodel=medium`` or ``-mcmodel=large``.

.. Another useful flag is related to implicit typesetting. 
.. Nek5000 relies often on implicit typesetting as default in the example cases. 
.. This means in practice that if the user defines a new variable in the user file and forgets to define its type explicitly then variable beginning with a character from I to N, its type is ``INTEGER``. 
.. Otherwise, it is ``REAL``.  
.. To avoid confusion the user not accustomed to implicit typesetting may use the warning flag ``-Wimplicit``. 
.. This flag warns whenever a variable, array, or function is implicitly declared and has an effect similar to using the ``IMPLICIT NONE`` statement in every program unit.

----------------------------------
Internal Input Parameters/Switches
----------------------------------

....................
Parameters
....................

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

| **P001**  density for the case of constant properties
|
| **P002**  dynamic viscosity 
|
| **P007**  heat capacity for the case of constant properties 
|
| **P008**  conductivity for the case of constant properties
|
| **P010**  simulation end time
| 
| **P011**  number of time steps
| 
| **P012**  time step size
| 
| **P014**  time fequency to dump fld files
| 
| **P015**  step frequency to dump fld files
| 
| **P021**  pressure solver tolernace
| 
| **P022**  velocity solver tolerance
| 
| **P023**  number of passive scalars
| 
| **P024**  relative tolerance for Helmholtz solver
| 
| **P025**  absolute tolerance for Helmholtz solver
| 
| **P026**  target Courant number (determines the number of RK4 substeps for OIFS)
| 
| **P027**  temporal discretization order
| 
| **P028**  temporal discretization order for mesh solver
| 
| **P029**  magnetic viscosity
| 
| **P030**  material properties (0: constant, 1: uservp) 
| 
| **P031**  number of perturbation modes in linearized N-S.
| 
| **P032**  number of boundary conditions in .re2 file
| 
| **P033**  first field index in .re2
| 
| **P040**  pressure coarse grid solver (0: XXT, 1: AMG)
| 
| **P041** 1 :math:`\rightarrow` multiplicative SEMG
| 
| **P042** linear solver for the pressure equation (0: GMRES, 1: CG)
|
| **P043** 0: additive multilevel scheme - 1: original two level scheme.
| 
| **P044** 0=E-based additive Schwarz for PnPn-2; 1=A-based.
| 
| **P045** Free-surface stability control (defaults to 1.0)
| 
| **P046** if :math:`>0`, do not set Initial Condition (no call to subroutine ``SETICS``).
| 
| **P047** Poisson ratio for mesh elasticity solve (default 0.4)
| 
| **P054** direction of fixed flowrate (1: x, 2: y, 3: z), negative means fixed bulk
| 
| **P055** volumetric flowrate or bulk velocity (see p054) for periodic case
| 
| **P059** deformed element switch
| 
| **P060** initialize velocity to 1e-10 (for steady Stokes problem).
| 
| **P062** byte swap for output
| 
| **P063** output precision (4: SP, 8: DP)
| 
| **P064** restart perturbation solution
| 
| **P065** number of I/O nodes (if :math:`< 0` write in separate subdirectories).
| 
| **P066** output format (0: ASCII, 4: legacy binary, 6: binary)
| 
| **P067** read format 
| 
| **P068** averaging frequency in ``avg_all`` (0: every timestep).
| 
| **P084** custom inital time step
| 
| **P086** use skew-symmetric instead of convective form.
| 
| **P093** number of previous solutions to use for residual projection.
| 
| **P094** number of steps starting residual projection for velocity and passive scalars
| 
| **P095** number of steps starting residual projection for pressure 
| 
| **P099** dealiasing mode (:math:`<0`: disabled, 3: old dealiasing, 4: new dealiasing)
| 
| **P100** :red:`RESERVED!` pressure preconditioner when using CG solver (0: Jacobi, :math:`>0`: two-level Schwarz) :red:`or viseversa?`
| 
| **P101** number of additional modes to filter
| 
| **P103** filter weight for last mode
| 
| **P107** if :math:`\neq0`, add it to ``h2`` in ``sethlm``
| 
| **P116 NELX** number of elements in :math:`x` for FTP
| 
| **P117 NELY** number of elements in :math:`y` for FTP
| 
| **P118 NELZ** number of elements in :math:`z` for FTP
| 


.. _sec:switches:

................
Logical switches
................

**IFFLOW** solve for fluid (velocity, pressure)

**IFHEAT** solve for heat (temperature and/or scalars)

**IFTRAN** solve transient equations (otherwise, solve the steady Stokes flow)

**IFADVC** specify the fields with convection

**IFTMSH** specify the field(s) defined on T mesh  (first field is the ALE mesh)

**IFAXIS** axisymmetric formulation

**IFSTRS** use stress formulation

**IFLOMACH** use low Mach number formulation

**IFMGRID** moving grid

**IFMVBD** moving boundary (for free surface flow)

**IFCHAR** use characteristics for convection operator

**IFSYNC** use upfront synchronization

**IFUSERVP** user-defined properties

.....................
Other Inut Varaibles
.....................

**filterType** 0: explicit, 1: HPT-RT

**restol(:)** field solver tolerance 

------------------------------
Commonly used Variables
------------------------------

..................
Solution Variables
..................

.. table::

  +---------------+------------------------------------+---------+------------------------------------------+
  | Variable Name | Size                               | Type    | Short Description                        |
  +===============+====================================+=========+==========================================+
  | ``vx``        | (lx1,ly1,lz1,lelv)                 | real    | x-velocity (u)                           |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vy``        | (lx1,ly1,lz1,lelv)                 | real    | y-velocity (v)                           |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vz``        | (lx1,ly1,lz1,lelv)                 | real    | z-velocity (w)                           |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``pr``        | (lx2,ly2,lz2,lelv)                 | real    | pressure (pr)                            |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``t``         | (lx1,ly1,lz1,lelt,ldimt)           | real    | temperature (t) and passive scalars (ps) |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vtrans``    | (lx1,ly1,lz1,lelt,ldimt+1)         | real    | convective coefficient                   |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vdiff``     | (lx1,ly1,lz1,lelt,ldimt+1)         | real    | diffusion coefficient                    |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vxlag``     | (lx1,ly1,lz1,lelv,2)               | real    | x-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vylag``     | (lx1,ly1,lz1,lelv,2)               | real    | y-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vzlag``     | (lx1,ly1,lz1,lelv,2)               | real    | z-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``prlag``     | (lx2,ly2,lz2,lelv,lorder2)         | real    | pressure at previous time steps          |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``tlag``      | (lx1,ly1,lz1,lelv,lorder-1,ldimt+1)| real    | t and ps at previous time steps          |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``time``      | --                                 | real    | physical time                            |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``dt``        | --                                 | real    | time step size                           |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``dtlag``     | ( 10 )                             | real    | previous time step sizes                 |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``istep``     | --                                 | integer | time step number                         |
  +---------------+------------------------------------+---------+------------------------------------------+

..................
Geometry Variables
..................

.. table::

  +---------------+---------------------------+-------------+-------------------------------------------+
  | Variable Name | Size                      | Type        | Short Description                         |
  +===============+===========================+=============+===========================================+
  | ``xm1``       | (lx1,ly1,lz1,lelt)        | real        | x-coordinates for velocity mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``ym1``       | (lx1,ly1,lz1,lelt)        | real        | y-coordinates for velocity mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``zm1``       | (lx1,ly1,lz1,lelt)        | real        | z-coordinates for velocity mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``bm1``       | (lx1,ly1,lz1,lelt)        | real        | mass matrix for velocity mesh             |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``binvm1``    | (lx1,ly1,lz1,lelv)        | real        | inverse mass matrix for velocity mesh     |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``bintm1``    | (lx1,ly1,lz1,lelt)        | real        | inverse mass matrix for t mesh            |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``volvm1``    | --                        | real        | total volume for velocity mesh            |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``voltm1``    | --                        | real        | total volume for t mesh                   |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``xm2``       | (lx2,ly2,lz2,lelv)        | real        | x-coordinates for pressure mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``ym2``       | (lx2,ly2,lz2,lelv)        | real        | y-coordinates for pressure mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``zm2``       | (lx2,ly2,lz2,lelv)        | real        | z-coordinates for pressure mesh           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``unx``       | (lx1,ly1,6,lelt)          | real        | x-component of face unit normal           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``uny``       | (lx1,ly1,6,lelt)          | real        | y-component of face unit normal           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``unz``       | (lx1,ly1,6,lelt)          | real        | z-component of face unit normal           |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``area``      | (lx1,ly1,6,lelt)          | real        | face area (surface integral weights)      |
  +---------------+---------------------------+-------------+-------------------------------------------+

.......................
Problem Setup Variables
.......................

.. table::

  +---------------+---------------------------+-------------+-------------------------------------------+
  | Variable Name | Size                      | Type        | Short Description                         |
  +===============+===========================+=============+===========================================+
  | ``nid``       | --                        | integer     | MPI rank id (lowest rank is always 0)     |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``nio``       | --                        | integer     | I/O node id                               |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``nelv``      | --                        | integer     | number of elements in velocity mesh       |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``nelt``      | --                        | integer     | number of elements in t mesh              |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``ndim``      | --                        | integer     | dimensionality of problem (i.e. 2 or 3)   |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``nsteps``    | --                        | integer     | number of time steps to run               |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``iostep``    | --                        | integer     | time steps between data output            |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``cbc``       | (6,lelt,ldimt+1)          | character*3 | boundary condition                        |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``lglel``     | (lelt)                    | integer     | local to global element number map        |
  +---------------+---------------------------+-------------+-------------------------------------------+
  | ``gllel``     | (lelg)                    | integer     | global to local element number map        |
  +---------------+---------------------------+-------------+-------------------------------------------+

...................
Averaging Variables
...................

Arrays associated with the ``avg_all`` subroutine

.. table::

  +---------------+---------------------------+---------+-----------------------------------------------+
  | Variable Name | Size                      | Type    | Short Description                             |
  +===============+===========================+=========+===============================================+
  | ``uavg``      | (ax1,ay1,az1,lelt)        | real    | time averaged x-velocity                      |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``vavg``      | (ax1,ay1,az1,lelt)        | real    | time averaged y-velocity                      |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``wavg``      | (ax1,ay1,az1,lelt)        | real    | time averaged z-velocity                      |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``pavg``      | (ax2,ay2,az2,lelt)        | real    | time averaged pressure                        |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``tavg``      | (ax1,ay1,az1,lelt,ldimt)  | real    | time averaged temperature and passive scalars |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``urms``      | (ax1,ay1,az1,lelt)        | real    | time averaged u^2                             |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``vrms``      | (ax1,ay1,az1,lelt)        | real    | time averaged v^2                             |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``wrms``      | (ax1,ay1,az1,lelt)        | real    | time averaged w^2                             |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``prms``      | (ax1,ay1,az1,lelt)        | real    | time averaged pr^2                            |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``trms``      | (ax1,ay1,az1,lelt,ldimt)  | real    | time averaged t^2 and ps^2                    |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``uvms``      | (ax1,ay1,az1,lelt)        | real    | time averaged uv                              |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``vwms``      | (ax1,ay1,az1,lelt)        | real    | time averaged vw                              |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``wums``      | (ax1,ay1,az1,lelt)        | real    | time averaged wu                              |
  +---------------+---------------------------+---------+-----------------------------------------------+
  | ``iastep``    | --                        | integer | time steps between averaged data output       |
  +---------------+---------------------------+---------+-----------------------------------------------+

.. include:: routines.rst

.. include:: geometry.rst

