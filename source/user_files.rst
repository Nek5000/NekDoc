.. _user_files:

==========
Case Files
==========

Each simulation is defined by six case files: 

- SESSION.NAME
- par
- re2
- usr
- SIZE
- map/ma2
- fld/f%05d

.. _user_files_session:

------------
SESSION.NAME
------------

To run Nek5000, each simulation must have a SESSION.NAME file. This file is read in by the code and
gives the path to the relevant files describing the structure and parameters of the simulation. The
SESSION.NAME file is a file that contains the name of the simulation and the full path to
supporting files. For example, to run the eddy example from the repository, the SESSION.NAME file
would look like:

.. code-block:: none

  eddy_uv
  /home/user_name/Nek5000/short_tests/eddy/ 

.. _user_files_par:

-----------------------------------
par
-----------------------------------

The simulation paramaters are defined in the ``.par`` file.
The keys are grouped in different sections and a specific value is assigned to each key.
The ``.par`` file follows the structure exemplified below.

.. code-block:: none

   #
   # nek parameter file
   #

   [SECTION]
   key = value
   ...

   [SECTION]
   key = value
   ...

The sections are:

* ``GENERAL`` (mandatory)
* ``PROBLEMTYPE``
* ``MESH``
* ``VELOCITY``
* ``PRESSURE`` (required for velocity)
* ``TEMPERATURE`` 
* ``SCALAR%%`` 
* ``CVODE``

When scalars are used, the keys of each scalar are defined under the section ``SCALAR%%`` varying 
between ``SCALAR01`` and ``SCALAR99``. The descripton of the keys of each section is given in the 
following tables (all keys/values are case insensitive). The value assigned to each key can be a 
user input (e.g. a <real> value) or one of the avaliable options listed in the tables below.
Values in parentheses denote the default value.


.. _tab:generalparams:

.. table:: ``GENERAL`` keys in the ``.par`` file

   +-------------------------+-----------------+----------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                |
   +=========================+=================+==============================================+
   | ``startFrom``           | | ``<string>``  | | Absolute/relative path of the field file   |
   |                         |                 | | to restart the simulation from             |
   +-------------------------+-----------------+----------------------------------------------+
   | ``stopAt``              | | ``(numSteps)``| | Stop mode                                  |
   |                         | | ``runTime``   |                                              |
   +-------------------------+-----------------+----------------------------------------------+
   | ``endTime``             | | ``<real>``    | | Final physical time at which we want to    |
   |                         |                 | | our simulation to stop                     |
   +-------------------------+-----------------+----------------------------------------------+
   | ``numSteps``            | | ``<real>``    | | Number of time steps instead of specifying |
   |                         |                 | | final physical time                        |
   +-------------------------+-----------------+----------------------------------------------+
   | ``dt``                  | | ``<real>``    | | Specifies the step size or in case of a    |
   |                         |                 | | a variable time step the maximum step size | 
   +-------------------------+-----------------+----------------------------------------------+
   | ``variableDT``          | | ``(no)``      | | Controls if the step size will be adjusted |
   |                         | | ``yes``       | | to match the targetCFL                     |
   +-------------------------+-----------------+----------------------------------------------+
   | ``targetCFL``           | | ``<real>``    | | Sets stability/target CFL number for       |
   |                         |                 | | OIFS or variable time steps                |
   |                         |                 | | (fixed to 0.5 for standard extrapolation   | 
   +-------------------------+-----------------+----------------------------------------------+
   | ``writeControl``        | | ``(timeStep)``| | Specifies whether checkpointing is based   |
   |                         | | ``runTime``   | | on number of time steps or physical time   |
   +-------------------------+-----------------+----------------------------------------------+
   | ``writeInterval``       | | ``<real>``    | | Checkpoint frequency in time steps or      | 
   |                         |                 | | physical time                              | 
   +-------------------------+-----------------+----------------------------------------------+
   | ``filtering``           | | ``(none)``    | | Specifies the filtering method             | 
   |                         | | ``explicit``  |                                              | 
   |                         | | ``hpfrt``     |                                              | 
   +-------------------------+-----------------+----------------------------------------------+
   | ``filterCutoffRatio``   | | ``<real>``    | | Ratio of modeal modes not affected         |
   |                         |                 | | Use i.e. for stabilization or LES 0.9/0.65 |  
   +-------------------------+-----------------+----------------------------------------------+
   | ``filterWeight``        | | ``<real>``    | | Sets the filter strength of transfer       |
   |                         |                 | | function of the last mode (explicit) or the|
   |                         |                 | | relaxation parameter in case of hpfrt      |  
   +-------------------------+-----------------+----------------------------------------------+
   | ``writeDoublePrecision``| | ``no``        | | Sets the precision of the field files      |
   |                         | | ``(yes)``     |                                              |
   +-------------------------+-----------------+----------------------------------------------+
   | ``writeNFiles``         | | ``(1)``       | | Sets the number of output files            | 
   |                         |                 | | By default a parallel shared file is used  |
   +-------------------------+-----------------+----------------------------------------------+
   | ``dealiasing``          | | ``no``        | | Enable/diasble over-integration            |
   |                         | | ``(yes)``     |                                              |
   +-------------------------+-----------------+----------------------------------------------+
   | ``timeStepper``         | | ``BDF1``      | | Time integration order                     |
   |                         | | ``(BDF2)``    |                                              |
   |                         | | ``BDF3``      |                                              |
   +-------------------------+-----------------+----------------------------------------------+
   | ``extrapolation``       | | ``(standard)``| | Extrapolation method                       |
   |                         | | ``OIFS``      |                                              |
   +-------------------------+-----------------+----------------------------------------------+
   | ``optLevel``            | | ``(2)``       | | Optimization level                         |
   +-------------------------+-----------------+----------------------------------------------+
   | ``logLevel``            | | ``(2)``       | | Verbosity level                            |
   +-------------------------+-----------------+----------------------------------------------+
   | ``userParam%%``         | | ``<real>``    | | User parameter (can be accessed through    |
   |                         |                 | | uparam(%) array in ``.usr``                |
   +-------------------------+-----------------+----------------------------------------------+



.. _tab:probtypeparams:

.. table:: ``PROBLEMTYPE`` keys in the ``.par`` file

   +---------------------------+---------------------+--------------------------------------------------+
   |   Key                     | | Value(s)          | | Description                                    |
   +===========================+=====================+==================================================+
   | ``equation``              | | ``(incompNS)``    | | Specifies equation type                        |
   |                           | | ``lowMachNS``     |                                                  |
   |                           | | ``steadyStokes``  |                                                  |
   |                           | | ``incompLinNS``   |                                                  |
   |                           | | ``incompLinAdjNS``|                                                  |
   |                           | | ``incompMHD``     |                                                  |
   |                           | | ``compNS``        |                                                  |
   |                           |                     |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``axiSymmetry``           | | ``(no)``          | | Axisymmetric problem                           |
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``swirl``                 | | ``(no)``          | | Enable axisymmetric azimuthal velocity         |
   |                           | | ``yes``           | | component (stored in temperature field         |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``cyclicBoundaries``      | | ``(no)``          | | Sets cyclic periodic boundaries                | 
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``numberOfPerturbations`` | | ``(1)``           | | Number of perturbations for linearized NS      |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``solveBaseFlow``         | | ``(no)``          | | Solve for base flow in case of linearized NS   |
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``variableProperties``    | | ``(no)``          | | Enable variable transport properties           |
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``stressFormulation``     | | ``(no)``          | | Enable stress formulation                      |
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``dp0dt``                 | | ``(no)``          | | Enable time-varying thermodynamic pressure     |
   |                           | | ``yes``           |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+

.. _tab:fieldparams:

.. table:: ``COMMON`` keys for all field variables in the ``.par`` file

   +-------------------------+-----------------+-------------------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                         |
   +=========================+=================+=======================================================+
   | ``residualTol``         | | ``<real>``    | | Residual tolerance used by solver (not for CVODE)   | 
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``residualProj``        | | ``(no)``      | | Controls the residual projection                    |
   |                         | | ``yes``       |                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``writeToFieldFile``    | | ``no``        | | Controls if fields will be written on output        |
   |                         | | ``(yes)``     |                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+

.. _tab:meshparams:

.. table:: ``MESH`` keys in the ``.par`` file

   +-------------------------+-----------------+-------------------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                         |
   +=========================+=================+=======================================================+
   | ``motion``              | | ``(none)``    | | Mesh motion solver                                  |
   |                         | | ``user``      |                                                       |
   |                         | | ``elasticity``|                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``viscosity``           | | ``(0.4)``     | | Diffusivity for elasticity solver                   |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``numberOfBCFields``    | | ``(nfields)`` | | Number of field variables which have a boundary     |
   |                         |                 | |  condition in ``.re2`` file                         |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``firstBCFieldIndex``   | | ``(1 or 2)``  | | Field index of the first BC specified in ``.re2``   |
   |                         |                 | | file                                                |
   +-------------------------+-----------------+-------------------------------------------------------+

.. _tab:velocityparams:

.. table:: ``VELOCITY`` keys in the ``.par`` file

   +-------------------------+--------------+------------------------------------------------+
   |   Key                   | | Value(s)   | | Description                                  |
   +=========================+==============+================================================+
   | ``viscosity``           | | ``<real>`` | | Dynamic viscosity                            |
   |                         |              | | A negative value sets the Reynolds number    |
   +-------------------------+--------------+------------------------------------------------+
   | ``density``             | | ``<real>`` | | Density                                      |
   +-------------------------+--------------+------------------------------------------------+

.. _tab:pressureparams:

.. table:: ``PRESSURE`` keys in the ``.par`` file

   +-------------------------+------------------+-----------------------------------------------+
   |   Key                   | | Value(s)       | | Description                                 |
   +=========================+==================+===============================================+
   | ``preconditioner``      | | ``(semg_xxt)`` | | Preconditioning method                      |
   |                         | | ``semg_amg``   | | First time usage of AMG will write three    |
   |                         |                  | | dump files to disc. Subsequently please run |
   |                         |                  | | the amg_hypre tool to create the setup files|
   |                         |                  | | required for the AMG solver initialization  |
   +-------------------------+------------------+-----------------------------------------------+

.. _tab:fieldparams:

.. table:: ``COMMON`` keys for temperature and scalar fields in the ``.par`` file

   +-------------------------+--------------+--------------------------------------------+
   |   Key                   | | Value(s)   | | Description                              |
   +=========================+==============+============================================+
   | ``solver``              | | ``(helm)`` | | Solver for scalar                        | 
   |                         | | ``cvode``  |                                            |  
   |                         | | ``none``   |                                            |
   +-------------------------+--------------+--------------------------------------------+
   | ``advection``           | | ``no``     | | Controls if advection is present         |
   |                         | | ``(yes)``  |                                            |
   +-------------------------+--------------+--------------------------------------------+
   | ``absoluteTol``         | | ``<real>`` | | Absolute tolerance used by CVODE         |
   +-------------------------+--------------+--------------------------------------------+

.. _tab:temperatureparams:

.. table:: ``TEMPERATURE`` keys in the ``.par`` file

   +--------------------------+--------------+----------------------------------------------+
   |   Key                    | | Value(s)   | | Description                                |
   +==========================+==============+==============================================+
   |``ConjugatedHeatTransfer``| | ``(no)``   | | Controls conjugate heat transfer           |
   |                          | | ``yes``    |                                              |
   +--------------------------+--------------+----------------------------------------------+
   | ``conductivity``         | | ``<real>`` | | Thermal conductivity                       |
   +--------------------------+--------------+----------------------------------------------+
   | ``rhoCp``                | | ``<real>`` | | Product of density and heat capacity       |
   +--------------------------+--------------+----------------------------------------------+

.. _tab:scalarparams:

.. table:: ``SCALAR%%`` keys in the ``.par`` file

   +--------------------------+----------------+--------------------------------------------+
   |   Key                    | | Value(s)     | | Description                              |
   +==========================+================+============================================+
   | ``density``              | | ``<real>``   | | Density                                  |
   +--------------------------+----------------+--------------------------------------------+
   | ``diffusivity``          | | ``<real>``   | | Diffusivity                              | 
   +--------------------------+----------------+--------------------------------------------+

.. _tab:cvodeparams:

.. table:: ``CVODE`` keys in the ``.par`` file

   +--------------------------+----------------+----------------------------------------------+
   |   Key                    | | Value(s)     | | Description                                |
   +==========================+================+==============================================+
   | ``relativeTol``          | | ``<real>``   | | Relative tolerance (applies to all scalars)|
   +--------------------------+----------------+----------------------------------------------+
   | ``stiff``                | | ``no``       | | Controls if BDF or Adams Moulton is used   |
   |                          | | ``(yes)``    |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``preconditioner``       | | ``(none)``   | | Preconditioner method                      |
   |                          | | ``user``     |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``dtMax``                | | ``<real>``   | | Maximum internal step size                 |
   |                          |                | | Controls splitting error of velocity       |
   |                          |                | | scalar coupling (e.g. set to 1-4 dt)       |
   +--------------------------+----------------+----------------------------------------------+


.. _user_files_re2:
-----------------------------------
re2
-----------------------------------

Stores the mesh and boundary condition. 
TODO: Add more details


.. _user_files_usr:

----------------------
usr
----------------------

This file implements the the user interface to Nek5000. What follows is a brief description of the available
subroutines. 

...................
uservp()
...................

This function can be used  to specify customized or solution dependent material
properties.  

Example:

.. code-block:: fortran

      if (ifield.eq.1) then
         udiff  = a * exp(-b*temp) ! dynamic viscosity
         utrans = 1.0              ! density
      else if (ifield.eq.2) then
         udiff  = 1.0              ! conductivity 
         utrans = 1.0              ! rho*cp 
      endif

...................
userf()
...................

This functions sets the source term (which will be subsequently be multiplied by 
the density) for the momentum equation. 

Example:

.. code-block:: fortran

      parameter(g = 9.81)

      ffx = 0.0 
      ffy = 0.0
      ffz = -g ! gravitational acceleration 
 
...................
userf()
...................

This functions sets the source term for the scalar equation. 

...................
userbc()
...................

This functions sets boundary conditions. Note, this function is only called
for special boundary condition types and only for points on the boundary surface.   

...................
userchk()
...................

This is a general purpose function that gets executed before the time stepper and after every time
step.

...................
userqtl()
...................

This function can be used  to specify a cutomzized thermal diveregence for the low Mach solver.
step.

...................
usrdat()
...................

This function can be used to modify the element vertices.

...................
usrdat2()
...................

This function can be used to modify the spectral element mesh.

...................
usrdat3()
...................

This function can be used to initialize case/user specific data.


------------------------
SIZE
------------------------

SIZE file defines the problem size, i.e. spatial points at which the solution is to be evaluated within each element, number of elements per processor etc.
The SIZE file governs the memory allocation for most of the arrays
in Nek5000, with the exception of those required by the C utilities.
The *basic* parameters of interest in SIZE are:

* **ldim** = 2 or 3.  This must be set to 2 for two-dimensional or axisymmetric simulations  (the latter only partially supported) or to 3 for three-dimensional simulations.
* **lx1** controls the polynomial order of the solution, :math:`N = {\tt lx1-1}`.
* **lxd** controls the polynomial order of the (over-)integration/dealiasing. Strictly speaking :math:`{\tt lxd=3 * lx1/2}` is required but often smaller values are good enough.
* **lx2** = ``lx1`` or ``lx1-2`` and is an approximation order for pressure that determines the formulation for the Navier-Stokes  solver (i.e., the choice between the :math:`\mathbb{P}_N - \mathbb{P}_N` or :math:`\mathbb{P}_N - \mathbb{P}_{N-2}` spectral-element methods). 
* **lelg**, an upper bound on the total number of elements in your mesh. 
* **lpmax**, a maximum number of processors that can be used
* **lpmin**, a minimum number of processors that can be used (see also  **Memory Requirements**).
* **ldimt**, an upper bound on a number of auxilary fields to solve (temperature + other scalars, minimum is 1).

The *optional*
upper bounds on parameters in SIZE are (minimum being 1 unless otherwise noted):

* **lhis**, a maximum history (i.e. monitoring) points.
* **maxobj**, a maximum number of objects.
* **lpert**, a maximum perturbations.
* **toteq**, a maximum number of conserved scalars in CMT (minimum could be 0).
* **nsessmax**, a maximum number of (ensemble-average) sessions.
* **lxo**, a maximum number of points per element for field file output (:math:`{\tt lxo \geq lx1}`).
* **lelx**, **lely**, **lelz**, a maximum number of element in each direction for global tensor product solver and/or dimentions.
* **mxprev**, a maximum dimension of projection space (e.g. 20).
* **lgmres**, a maximum dimension of Krylov space (e.g. 30).
* **lorder**, a maximum order of temporal discretization (minimum is2 see also characteristic/OIFS method).
* **lelt** determines the maximum number of elements *per processor* (should be not smaller than nelgt/lpmin, e.g. lelg/lpmin+1).
* **lx1m**, a polynomial order for mesh solver that should be equal to lx1 in case of ALE and in case of stress-formulation (=1 otherwise).
* **lbelt** determines the maximum number of elements per processor for MHD solver that should be equalt to lelt (=1 otherwise).
* **lpelt** determines the maximum number of elements per processor for linear stability solver that should be equalt to lelt (=1 otherwise).
* **lcvelt** determines the maximum number of elements per processor for CVODE solver that should be equalt to lelt (=1 otherwise).
* **lfdm** equals to 1 for global tensor product solver (that uses fast diagonalization method) being 0 otherwise.

Note that one also need to include the following line to SIZE file:

.. code-block:: fortran

      include 'SIZE.inc'

that defines addional internal parameters.


.. _user_files_ma2:
-----------------------------------
map/ma2
-----------------------------------

TODO: Add more details


.. _user_files_fld:
-----------------------------------
fld/f%05d
-----------------------------------

TODO: Add fld details

The binary ``.f%05d`` file format is used to write and read data both in serial and parallel
in Nek5000.

The file is composed of:

  - header
  - mesh data
  - field data
  - bounding box data

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
   | 1     | 2       | ``wdsizo``  | sets the precision to 4 or 8                  |
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

    #std 4  6  6  1         36         36  0.1000000000000E+03     10000     0      1 XUP                                          

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


TODO: Add more details
