==========
Appendices
==========

-------------------------------
List of Parameters in .rea File
-------------------------------

..........
Parameters
..........

This section tells Nek5000

- If the input file reflects a 2D or 3D job (it should match the ``ldim`` parameter in the SIZE file).
- The combination of heat transfer, Stokes, Navier-Stokes, steady or unsteady to be run.
- The relevant physical parameters.
- The solution algorithm within Nek5000 to use.
- The timestep size or Courant number to use, or whether to run variable DT (:math:`dt`), etc.

A ``.rea`` file starts with the following three parameters:

**NEKTON VERSION** the version of Nek5000

**DIMENSIONAL RUN** number of spatial dimensions (``NDIM`` =2,3 - has to match the setting in the SIZE file).

**PARAMETERS FOLLOW** the number of parameters which are going to be followed in the ``.rea`` file.(``NPARAM``)

The latter specifies how many lines of ``.rea`` file, starting from the next line, are the parameters and have to be read by the program.

....................
Available Parameters
....................

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

| **P001  DENSITY** density for the case of constant properties (for variable density see parameter ``P030``).
|
| **P002  VISCOS**  kinematic viscosity (if :math:`<0 \rightarrow Re`, otherwise :math:`1/Re`).
|
| **P003  BETAG** if :math:`>0`, natural convection is turned on (Boussinesq approximation). :red:`NOT IN USE!`
|
| **P004  GTHETA** model parameter for Boussinesq approximation (see parameter P003). :red:`NOT IN USE!`
|
| **P005  PGRADX** :red:`NOT IN USE!`
|
| **P006** :red:`NOT IN USE!`
|
| **P007  RHOCP** ``navier5.f:      param(7) = param(1)  ! rhoCP   = rho`` :red:`NOT IN USE!`
|
| **P008  CONDUCT** conductivity for the case of constant properties (if :math:`<0`, it defines the Peclet number, see parameter P030).
|    ``connect2.f:      if(param(8) .lt.0.0) param(8)  = -1.0/param(8)``
|    ``navier5.f:      param(8) = param(2)  ! conduct = dyn. visc``
|
| **P009** :red:`NOT IN USE!` (passed to ``CPFLD(2,3)``!)
|    ``connect2.f:      CPFLD(2,3)=PARAM(9)``
|
| **P010  FINTIME** if :math:`>0`, specifies simulation end time. Otherwise, use ``NSTEP`` (P011).
|    ``drive2.f:      FINTIM = PARAM(10)``
| 
| **P011  NSTEP** number of time steps.
|     ``connect2.f:            param(11) = 1.0``
|     ``drive2.f:      NSTEPS = PARAM(11)``
| 
| **P012  DT** upper bound on time step :math:`dt`   (if :math:`<0`, then :math:`dt=|P012|` constant)
|     ``connect2.f:            param(12) = 1.0``
|     ``drive2.f:      DT     = abs(PARAM(12))``
| 
| **P013  IOCOMM** frequency of iteration histories
|     ``drive2.f:      IOCOMM = PARAM(13)``
| 
| **P014  IOTIME** if :math:`>0`, time interval to dump the fld file. Otherwise, use ``IOSTEP`` (P015).
|     ``drive2.f:      TIMEIO = PARAM(14)``
| 
| **P015  IOSTEP** dump frequency, number of time steps between dumps.
|     ``drive2.f:      IOSTEP = PARAM(15)``
|     ``navier5.f:      if  (iastep.eq.0) iastep=param(15)   ! same as iostep``
| 
| **P016  PSSOLVER** heat/passive scalar solver:
|    1. Helmholz
|    2. ``CVODE``
|    3. ``CVODE`` with user-supplied Jacobian
|    Note: a negative number will set source terms to 0.
| 
| **P017  AXIS**  :red:`NOT IN USE!`
| 
| **P018  GRID** :red:`NOT IN USE!`
| 
| **P019  INTYPE** :red:`NOT IN USE!`
|     ``connect2.f:            param(19) = 0.0``
| 
| **P020  NORDER**  :red:`NOT IN USE!`
| 
| **P021  DIVERGENCE** tolerance for the pressure solver.
|     ``drive2.f:      TOLPDF = abs(PARAM(21))``
|     ``hmholtz.f:      if (name.eq.'PRES'.and.param(21).ne.0) tol=abs(param(21))``
| 
| **P022  HELMHOLTZ** tolerance for the velocity solver.
|     ``drive2.f:      TOLHDF = abs(PARAM(22))``
|     ``hmholtz.f:      if (param(22).ne.0) tol=abs(param(22))``
|     ``hmholtz.f:         if (param(22).lt.0) tol=abs(param(22))*rbn0``
|     ``navier4.f:      if (param(22).ne.0) tol = abs(param(22))``
| 
| **P023  NPSCAL** number of passive scalars.
|     ``connect2.f:      NPSCAL=INT(PARAM(23))``
| 
| **P024  TOLREL** relative tolerance for the passive scalar solver (``CVODE``).
|     ``drive2.f:      TOLREL = abs(PARAM(24))``
| 
| **P025  TOLABS** absolute tolerance for the passive scalar solver (``CVODE``).
|     ``drive2.f:      TOLABS = abs(PARAM(25))``
| 
| **P026  COURANT** maximum Courant number (number of RK4 substeps if OIFS is used).
|     ``drive2.f:      CTARG  = PARAM(26)``
| 
| **P027  TORDER** temporal discretization order (2 or 3).
|     ``drive2.f:      NBDINP = PARAM(27)``
| 
| **P028  NABMSH** Order of temporal integration for mesh velocity. If 1, 2, or 3 use Adams-Bashforth of corresponding order. Otherwise, extrapolation of order ``TORDER`` (P027).
| 
| **P029  MHD_VISCOS** if :math:`>0 \rightarrow` magnetic viscosity, if :math:`<0 \rightarrow` magnetic Reynolds number.
|     ``connect2.f:      if(param(29).lt.0.0) param(29) = -1.0/param(29)``
|     ``connect2.f:      if (param(29).ne.0.) ifmhd  = .true.``
|     ``connect2.f:         cpfld(ifldmhd,1) = param(29)  ! magnetic viscosity``
| 
| **P030  USERVP** if
|    0. constant properties
|    1. user-defined properties via ``uservp`` subroutine (each scalar separately)
|    2. user-defined properties via ``uservp`` subroutine (all scalars at once)
| 
| **P031  NPERT**  if :math:`\neq 0`, number of perturbation modes in linearized N-S.
|     ``connect2.f:      if (param(31).ne.0.) ifpert = .true.``
|     ``connect2.f:      if (param(31).lt.0.) ifbase = .false.   ! don't time adv base flow``
|     ``connect2.f:      npert = abs(param(31))``
| 
| **P032  NBCRE2** if :math:`>0`, number of BCs in ``.re2`` file, 0: all.
|     ``connect2.f:      if (param(32).gt.0) nfldt = ibc + param(32)-1``
| 
| **P033** :red:`NOT IN USE!`
| 
| **P034** :red:`NOT IN USE!`
| 
| **P035** :red:`NOT IN USE!`
| 
| **P036 XMAGNET** :red:`NOT IN USE!`
| 
| **P037 NGRIDES** :red:`NOT IN USE!`
| 
| **P038 NORDER2** :red:`NOT IN USE!`
| 
| **P039 NORDER3** :red:`NOT IN USE!`
| 
| **P040** :red:`NOT IN USE!`
| 
| **P041** 1 :math:`\rightarrow` multiplicative SEMG
|     ``hsmg.f:c     if (param(41).eq.1) ifhybrid = .true.`` :math:`\leftarrow` :red:`NOT IN USE!`
| 
| **P042** linear solver for the pressure equation, 0 :math:`\rightarrow` GMRES or 1 :math:`\rightarrow` PCG
| 
| **P043** 0: additive multilevel scheme - 1: original two level scheme.
|     ``navier6.f:      if (lx1.eq.2) param(43)=1.``
|     ``navier6.f:            if (param(43).eq.0) call hsmg_setup``
| 
| **P044** 0=E-based additive Schwarz for PnPn-2; 1=A-based.
| 
| **P045** Free-surface stability control (defaults to 1.0)
|     ``subs1.f:      FACTOR = PARAM(45)``
| 
| **P046** if :math:`>0`, do not set Initial Condition (no call to subroutine ``SETICS``).
|    ``drive2.f:      irst = param(46)``
|    ``ic.f:      irst = param(46)        ! for lee's restart (rarely used)``
|    ``subs1.f:      irst = param(46)``
| 
| **P047** parameter for moving mesh (Poisson ratio for mesh elasticity solve (default 0.4)).
|     ``mvmesh.f:      VNU    = param(47)``
| 
| **P048** :red:`NOT IN USE!`
| 
| **P049** if :math:`<0`, mixing length factor :red:`NOT IN USE!`.
|     ``drive2.f:c     IF (PARAM(49) .LE. 0.0) PARAM(49) = TLFAC``
|     ``turb.f:      TLFAC = PARAM(49)``
| 
| **P050** :red:`NOT IN USE!`
| 
| **P051** :red:`NOT IN USE!`
| 
| **P052  HISTEP** if :math:`>1`, history points dump frequency (in number of steps).
|     ``prepost.f:      if (param(52).ge.1) iohis=param(52)``
| 
| **P053** :red:`NOT IN USE!`
| 
| **P054** direction of fixed mass flowrate (1: :math:`x`-, 2: :math:`y`-, 3: :math:`z`-direction). If 0: :math:`x`-direction.
|     ``drive2.f:      if (param(54).ne.0) icvflow = abs(param(54))``
|     ``drive2.f:      if (param(54).lt.0) iavflow = 1 ! mean velocity``
| 
| **P055** volumetric flow rate for periodic case;  if p54:math:`<0`, then p55:=mean velocity.
|     ``drive2.f:      flowrate = param(55)``
| 
| **P056** :red:`NOT IN USE!`
| 
| **P057** :red:`NOT IN USE!`
| 
| **P058** :red:`NOT IN USE!`
| 
| **P059** if :math:`\neq0`, deformed elements (only relevant for FDM). !=0 :math:`\rightarrow` full Jac. eval. for each el.
| 
| **P060** if :math:`\neq0`, initialize velocity to 1e-10 (for steady Stokes problem).
| 
| **P061** :red:`NOT IN USE!`
| 
| **P062** if :math:`>0`, swap bytes for output.
| 
| **P063  WDSIZO** real output wordsize (8: 8-byte reals, else 4-byte).
|     ``prepost.f:      if (param(63).gt.0) wdsizo = 8         ! 64-bit .fld file``
| 
| **P064** if :math:`=1`, restart perturbation solution
|     ``pertsupport.f:      if(param(64).ne.1) then !fresh start, param(64) is restart flag``
| 
| **P065** number of I/O nodes (if :math:`< 0` write in separate subdirectories).
| 
| **P066** Output format: (only ``postx`` uses ``.rea`` value; other nondefault should be set in ``usrdat``) (if :math:`\geq 0` binary else ASCII).
|     ``connect2.f:         param(66) = 6        ! binary is default``
|     ``connect2.f:         param(66) = 0        ! ASCII``
| 
| **P067** read format (if :math:`\geq 0` binary else ASCII).
| 
| **P068** averaging frequency in ``avg_all`` (0: every timestep).
| 
| **P069** :red:`NOT IN USE!`
| 
| **P070** :red:`NOT IN USE!`
| 
| **P071** :red:`NOT IN USE!`
| 
| **P072** :red:`NOT IN USE!`
| 
| **P073** :red:`NOT IN USE!`
| 
| **P074** if :math:`> 0` print Helmholtz solver iterations.
|     ``hmholtz.f:         if (nid.eq.0.and.ifprint.and.param(74).ne.0) ifprinthmh=.true.``
| 
| **P075** :red:`NOT IN USE!`
| 
| **P076** :red:`NOT IN USE!`
| 
| **P077** :red:`NOT IN USE!`
| 
| **P078** :red:`NOT IN USE!`
| 
| **P079** :red:`NOT IN USE!`
| 
| **P080** :red:`NOT IN USE!`
| 
| **P081** :red:`NOT IN USE!`
| 
| **P082** coarse-grid dimension (2: linear). :red:`NOT IN USE!`
| 
| **P083** :red:`NOT IN USE!`
| 
| **P084** if :math:`<0`, force initial time step to this value.
| 
| **P085** set :math:`dt` in *setdt*.
|     ``subs1.f:            dt=dtopf*param(85)i``
| 
| **P086** :red:`RESERVED!` if :math:`\neq0`, use skew-symmetric form, else convective form.
|     ``drive2.f:      PARAM(86) = 0 ! No skew-symm. convection for now``
|     ``navier1.f:      if (param(86).ne.0.0) then  ! skew-symmetric form``
| 
| **P087** :red:`NOT IN USE!`
| 
| **P088** :red:`NOT IN USE!`
| 
| **P089** :red:`RESERVED!`
| 
| **P090** :red:`NOT IN USE!`
| 
| **P091** :red:`NOT IN USE!`
| 
| **P092** :red:`NOT IN USE!`
| 
| **P093**  number of previous solutions to use for residual projection.
|    (adjust ``MXPREF`` in ``SIZEu`` accordingly)
| 
| **P094** if :math:`>0`, start projecting velocity and passive scalars after P094 steps
| 
| **P095** if :math:`>0`, start projecting pressure after P095 steps
| 
| **P096** :red:`NOT IN USE!`
| 
| **P097** :red:`NOT IN USE!`
| 
| **P098** :red:`NOT IN USE!`
| 
| **P099** dealiasing:
|    :math:`<0`:  disable
|    3:  old dealiasing
|    4:  new dealiasing
| 
| **P100** :red:`RESERVED!` pressure preconditioner when using CG solver (0: Jacobi, :math:`>0`: two-level Schwarz) :red:`or viseversa?`
| 
| **P101** number of additional modes to filter (0: only last mode)
|     ``navier5.f:         ncut = param(101)+1``
| 
| **P102** :red:`NOT IN USE!`
| 
| **P103** filter weight for last mode (:math:`<0`: disabled)
| 
| **P107** if :math:`\neq0`, add it to ``h2`` in ``sethlm``
| 
| **P116 NELX** number of elements in :math:`x` for Fast Tensor Product FTP solver (0: do not use FTP).
|    NOTE: box geometries, constant properties only!
| 
| **P117  NELY** number of elements in :math:`y` for FTP
| 
| **P118  NELZ** number of elements in :math:`z` for FTP

..........................
Available Logical Switches
..........................

This part of ``.rea`` file starts with such a line::

   n   LOGICAL SWITCHES FOLLOW

where ``n`` is the number of logical switches which is set in the following lines.

.. _sec:switches:

................
Logical switches
................

Note that by default all logical switches are set to false.

**IFFLOW** solve for fluid (velocity, pressure).

**IFHEAT** solve for heat (temperature and/or scalars).

**IFTRAN** solve transient equations (otherwise, solve the steady Stokes flow).

**IFADVC** specify the fields with convection.

**IFTMSH** specify the field(s) defined on T mesh  (first field is the ALE mesh).

**IFAXIS** axisymmetric formulation.

**IFSTRS** use stress formulation in the incompressible case.

**IFLOMACH** use low Mach number compressible flow.

**IFMGRID** moving grid (for free surface flow).

**IFMVBD** moving boundary (for free surface flow).

**IFCHAR** use characteristics for convection operator.

**IFSYNC** use mpi barriers to provide better timing information.

**IFUSERVP** user-defined properties (e.g., :math:`\mu`, :math:`\rho` varying with space and time.

-------------------------------
List of Parameters in SIZE File
-------------------------------

| **ldim**: number of spatial dimensions (2 or 3). 
| 
| **lx1**: number of (GLL) points in the :math:`x` -direction within each element of mesh1 (velocity) which is equal to the (polynomial order :math:`+1`) by definition. 
| 
| (``lx1`` recomeneded odd for better performance)
| 
| **lx2**: number of (GLL) points in the :math:`x` -directions within each element of mesh2 (pressure). Use ``lx2=lx1`` for PN/PN formulation or ``lx2=lx1-2`` for PN/PN-2 formulation.
| 
| **lxd**: number of points for over integration (dealiasing), use three half rule e.g. for ``lx1=8`` use ``lxd=12``.
| 
| **lelx, lely, lelz**: maximum number of elements per rank for global FDM (Fast Diagonalization Method) solver.
| 
| **ldimt**:  maximum number of T-array fields (temperature + additional scalars).
| 
| **lpmax**: maximum number of ranks.
|
| **lpmin**: minimum number of ranks. 
|
| **lelg**: maximum (global) number of elements (it is usually set more than the # of elements existing in the mesh, for making maximum use of memory is can be set to the exact number of mesh elements).
| 
| **lelt**: maximum number of local elements for T-mesh (per rank, ``lelt`` :math:`\geq` ``lelg/np +1``).
| 
| **lpelt**: Number of elements of the perturbation field, number of perturbation fields
| 
| **lbelt**: Total Number of elements of the B-field (MHD)
| 
| **lx1m**: when the mesh is a moving type ``lx1m=lx1``, otherwise it is set to 1.
| 
| **lorder**: maximum time integration order (2 or 3).
| 
| **maxobj**: maximum number of objects. :red:`zero if not using objects?`
| 
| **maxmbr**: maximum number of members in an object.
| 
| **lhis**: maximum number of history points a single rank will read in (``NP*LHIS`` :math:`<` number of points in ``hpts.in``).
| 
| **mxprev**: maximum number of history entries for residual projection (recommended value: 20).
| 
| **lgmres**: dimension of Krylov subspace in GMRES (recommended value: 40).

...........................
Parameters in SIZE.inc File
...........................

The following parameters appeared in the SIZE file in previous versions, and are now moved to the internal SIZE.inc file. They are automatically set based on SIZE parameters.


| **ly1, lz1**: number of (GLL) points in the :math:`y` and :math:`z`-directions, respectively, within each element of mesh1 (velocity) which is equal to the (polynomial order :math:`+1`) by definition. ``ly1`` is usually the same as ``lx1`` and for 2D cases ``lz1=1``.
| (is ``lx1`` :math:`\neq` ``ly1`` supported?)
| 

| **ly2, lz2**: number of (GLL) points in the :math:`y` and :math:`z` directions, respectively, within each element of mesh2 (pressure).
| 
| **lx3, ly3, lz3**: number of (GLL) points in the :math:`x`, :math:`y` and :math:`z` directions, respectively, within each element of mesh3. These are set to the same number of (GLL) point on mesh1.
| (mesh3 is rarely used)
| 
| **lyd, lzd**: number of points for over integration (dealiasing). ``lyd = lxd``, and ``lzd = lxd`` for 3D and ``lzd = 1`` for 2D.
| 
| **lp**: ``lp = lpmax``
| 
| **lelv**: maximum number of local elements for V-mesh (``lelv = lelt``).
| 
| **lpelv**: Number of elements of the perturbation field, number of perturbation fields. ``lpelv = lpelt``.
| 
| **lpx1, lpy1, lpz1**: Number of point in :math:`x`, :math:`y`, :math:`z` direction of perturbation field within each element of mesh1. ``lpx1 = lx1``, ``lpy1 = lpx1``, and ``lpz1 = lpx1`` for 3D and ``lpz1 = 1`` for 2D.
| 
| **lbelv**: Total Number of elements of the B-field (MHD). ``lbelv = lbelt``.
| 
| **lbx1, lby1, lbz1**: Number of point in :math:`x`, :math:`y`, :math:`z` direction of B-field within each element of mesh1. ``lbx1 = lx1``, ``lby1 = lbx1``, and ``lbz1 = lbx1`` for 3D and ``lbz1 = 1`` for 2D.
| 
| **lbx2, lby2, lbz2**: Number of point in :math:`x`, :math:`y`, :math:`z` direction of B-field within each element of mesh2. ``lbx2 = lx2``, ``lby2 = lbx2``, and ``lbz2 = lbx2`` for 3D and ``lbz2 = 1`` for 2D.
| 
| **ly1m, lz1m**: when the mesh is a moving type ``lx1m=lx1``, otherwise it is set to 1. ``ly1m = lx1m``, and ``lz1m = lx1m`` for 3D and ``lz1m = 1`` for 2D.
| 
| **lxz**: lxz = lx1*lz1
|      ``connect1.f:      common /scruz/  snx(lxz) , sny(lxz) , snz(lxz) ,  efc(lxz)``
| 
| **lctmp0**: ``lctmp0 = 2*lx1*ly1*lz1*lelt``
|      ``drive1.f:c      COMMON /CTMP0/ DUMMY0(LCTMP0)``
| 
| **lctmp1**: ``lctmp1 = 4*lx1*ly1*lz1*lelt``
|      ``drive1.f:c      COMMON /CTMP1/ DUMMY1(LCTMP1)``
|      ``drive2.f:      COMMON /SCRNS/ WORK(LCTMP1)``
| 
| **maxmor**: ``=lelt``
| 
| **lzl**: for 2D cases ``lzl=1`` and for 3D cases ``lzl=3`` (computed automatically).
 
.....................
Deprecated Parameters
.....................

The following parameters are deprecated and were subsequently removed in newer versions.

| **lpert**: Number of elements of the perturbation field, number of perturbation fields
|
| **lstore**: :red:`NOT IN USE!`
|
| **lsvec**: :red:`NOT IN USE!`
|
| **lmvec**: :red:`NOT IN USE !`
|
| **lvec**: :red:`NOT IN USE!`
|
| **lelgec**: ``lelgec = 1``
|
| **lxyz2**: ``lxyz2 = 1``
|
| **lxz21**: ``lxz21 = 1``
|
| **lmaxv**: ``lmaxv = lx1*ly1*lz1*lelv``
|
| **lmaxt**: ``lmaxt = lx1*ly1*lz1*lelt``
|
| **lmaxp**: ``lmaxp = lx1*ly1*lz1*lelv``

--------------------
The .fld File Format
--------------------

The ``fld`` file format is used to write and read data both in serial and parallel
in Nek5000. This section describes the format and should allow third party tool
developers to implement pre and postprocessing tools.

The file is composed of:

  - the *header* in ASCII format,
  - mesh data, including the geometry, saved unrolled as scalar vector
    fields

We will go through each of these categories and give a description of its
composition.

......
Header
......

The header provides structural information about the stored data that is needed
to parse it correctly. The header is composed of 11 values in ASCII format. It
has a fixed size of 132 bytes and starts with the string ``#std``. All
header entries are padded to the right. After the header with 132 bytes, 4 bytes
follow that determine the endianess of the binary file.  It is the binary
representation of the number 6.54321 either in little or big endian.

.. table::

   +-------+---------+-------------+-----------------------------------------------+
   | Entry | Padding |  Name       | Short Description                             |
   +=======+=========+=============+===============================================+
   | 1     | 2       | ``wdsizo``  | sets the precision to 4 (float) or 8 (double) |
   +-------+---------+-------------+-----------------------------------------------+
   | 2     | 3       | ``nx``      | number of coordinates in x direction          |
   +-------+---------+-------------+-----------------------------------------------+
   | 2     | 3       | ``ny``      | number of coordinates in y direction          |
   +-------+---------+-------------+-----------------------------------------------+
   | 2     | 3       | ``nz``      | number of coordinates in z direction          |
   +-------+---------+-------------+-----------------------------------------------+
   | 5     | 11      | ``nelo``    | number of elements                            |
   +-------+---------+-------------+-----------------------------------------------+
   | 6     | 11      | ``nelgt``   | :red:`----`                                   |
   +-------+---------+-------------+-----------------------------------------------+
   | 7     | 21      | ``time``    | time stamp                                    |
   +-------+---------+-------------+-----------------------------------------------+
   | 8     | 10      | ``iostep``  | time step                                     |
   +-------+---------+-------------+-----------------------------------------------+
   | 9     | 7       | ``fid0``    | :red:`field id`                               |
   +-------+---------+-------------+-----------------------------------------------+
   | 10    | 7       | ``nfileoo`` | :red:`number of files`                        |
   +-------+---------+-------------+-----------------------------------------------+
   | 11    | 4       | ``rdcode``  | Fields written                                |
   +-------+---------+-------------+-----------------------------------------------+

Example of a header:::

    #std 4  6  6  1         36         36  0.1000000000000E+03     10000     0      1 XUP                                              úaÑ@

``wdsize`` sets the precision of the floating point numbers in the file. This
is either 4 bytes for floats or 8 bytes for double precision.

``nx``, ``ny`` and ``nz`` set the number of coordinates in  :math:`x`, :math:`y` and :math:`z`
direction for each element (polynomial order), respectively. ``nelo`` sets
the number of total elements on the mesh.

``time`` is the simulation time while ``iostep`` is the time step when the file was written.

``rdcode`` determines which fields are contained in the file:

  - X: Geometry
  - U: Velocity
  - P: Pressure
  - T: Temperature
  - S: Passive scalar

....
Data
....

The data field begins after the first 136 bytes of the file. The values are
stored unrolled for each element and for each direction.
Example code for reading the geometry field in python:

.. code-block:: python

    for iel in range(nelo):
        x=ifilebuf.read(nxyzo8*wdsizo)
        xup=numpy.array(struct.unpack(nxyzo8*c,x),dtype=c)
        xfield[iel,:]=xup
        y=ifilebuf.read(nxyzo8*wdsizo)
        yup=numpy.array(struct.unpack(nxyzo8*c,y),dtype=c)
        yfield[iel,:]=yup
        if if3d:
            z=ifilebuf.read(nxyzo8*wdsizo)
            zup=numpy.array(struct.unpack(nxyzo8*c,z),dtype=c)
            zfield[iel,:]=zup
