.. _user_files:

==========
Case Files
==========

Each simulation is defined by six case files: 

- par
- re2
- usr
- SIZE
- ma2
- SESSION.NAME

.. _user_files_session:

------------
SESSION File
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
par File
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
* ``PRESSURE``
* ``TEMPERATURE`` 
* ``SCALAR%%`` 
* ``CVODE``

When scalars are used, the keys of each scalar are defined under the section ``SCALAR%%`` varying 
between ``SCALAR01`` and ``SCALAR99``. The descripton of the keys of each section is given in the 
following tables (all keys/values are case insensitive). The value assigned to each key can be a 
user input (e.g. a <real> value) or one of the avaliable options listed in the tables below.


.. _tab:generalparams:

.. table:: ``GENERAL`` keys in the ``.par`` file

   +-------------------------+---------------+----------------------------------------------+
   |   Key                   | | Value(s)    | | Description                                |
   +=========================+===============+==============================================+
   | ``startFrom``           | | ``<string>``| | Absolute/relative path of the field file   |
   |                         |               | | to restart the simulation from             |
   +-------------------------+---------------+----------------------------------------------+
   | ``stopAt``              | | ``endTime`` | | Stop mode                                  |
   |                         | | ``numStep`` |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``endTime``             | | ``<real>``  | | Final physical time at which we want to    |
   |                         |               | |  our simulation to stop                    |
   +-------------------------+---------------+----------------------------------------------+
   | ``numStep``             | | ``<real>``  | | Number of time steps instead of specifying |
   |                         |               | |  final physical time                       |
   +-------------------------+---------------+----------------------------------------------+
   | ``dt``                  | | ``<real>``  | | Specifies the step size or in case of a    |
   |                         |               |   a variable time step the maximum step size | 
   +-------------------------+---------------+----------------------------------------------+
   | ``variableDT``          | | ``yes``     | | Controls if the step size will be adjusted |
   |                         | | ``no``      | | to match the targetCFL                     |
   +-------------------------+---------------+----------------------------------------------+
   | ``targetCFL``           | | ``<real>``  | | Sets target CFL number                     |
   +-------------------------+---------------+----------------------------------------------+
   | ``writeControl``        | | ``runTime`` | | Specifies whether checkpointing is based   |
   |                         | | ``timeStep``| | on number of time steps or physical time   |
   +-------------------------+---------------+----------------------------------------------+
   | ``writeInterval``       | | ``<real>``  | | Checkpoint frequency in time steps or      | 
   |                         |               | | physical time                              | 
   +-------------------------+---------------+----------------------------------------------+
   | ``filtering``           | | ``explicit``| | Specifies the filtering method             | 
   |                         | | ``hpfrt``   |                                              | 
   +-------------------------+---------------+----------------------------------------------+
   | ``filterCutoffRatio``   | | ``<real>``  | | Ratio of modeal modes not affected         |
   |                         |               | | Use i.e. for stabilization or LES 0.9/0.65 |  
   +-------------------------+---------------+----------------------------------------------+
   | ``filterWeight``        | | ``<real>``  | | Sets the filter strength of transfer       |
   |                         |               | | function of the last mode (explicit) or the|
   |                         |               | | relaxation parameter in case of hpfrt      |  
   +-------------------------+---------------+----------------------------------------------+
   | ``writeDoublePrecision``| | ``yes``     | | Sets the precision of the field files      |
   |                         | | ``no``      |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``writeNFiles``         | | ``<real>``  | | Sets the number of output files            |  
   +-------------------------+---------------+----------------------------------------------+
   | ``dealiasing``          | | ``yes``     | | Enable/diasble over-integration            |
   |                         | | ``no``      |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``timeStepper``         | | ``BDF2``    | | Time integration order                     |
   |                         | | ``BDF3``    |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``extrapolation``       | | ``standard``| | Extrapolation method                       |
   |                         | | ``OIFS``    |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``optLevel``            | | ``<real>``  | | Optimization level (1-3,default is 1)      |
   +-------------------------+---------------+----------------------------------------------+
   | ``logLevel``            | | ``<real>``  | | Logfile verbosity level (1-2, default is 1)|
   +-------------------------+---------------+----------------------------------------------+
   | ``userParam%%``         | | ``<real>``  | | User parameter (can be accessed through    |
   |                         |               | | uparam(%) array in ``.usr``                |
   +-------------------------+---------------+----------------------------------------------+



.. _tab:probtypeparams:

.. table:: ``PROBLEMTYPE`` keys in the ``.par`` file

   +---------------------------+---------------------+--------------------------------------------------+
   |   Key                     | | Value(s)          | | Description                                    |
   +===========================+=====================+==================================================+
   | ``equation``              | | ``incompNS``      | | Specifies equation type                        |
   |                           | | ``lowMachNS``     |                                                  |
   |                           | | ``steadyStokes``  |                                                  |
   |                           | | ``incompLinNS``   |                                                  |
   |                           | | ``incompLinAdjNS``|                                                  |
   |                           | | ``incompMHD``     |                                                  |
   |                           | | ``compNS``        |                                                  |
   |                           |                     |                                                  |
   |                           |                     |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``axiSymmetry``           | | ``yes``           | | Axisymmetric problem                           |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``swirl``                 | | ``yes``           | | Enable axisymmetric azimuthal velocity         |
   |                           | | ``no``            | | component (stored in temperature field         |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``cyclicBoundaries``      | | ``yes``           | | Sets cyclic periodic boundaries                | 
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``numberOfPerturbations`` | | ``<real>``        | | Number of perturbations for linearized NS      |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``solveBaseFlow``         | | ``yes``           | | Solve for base flow in case of linearized NS   |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``variableProperties``    | | ``yes``           | | Enable variable transport properties           |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``stressFormulation``     | | ``yes``           | | Enable stress formulation                      |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``dp0dt``                 | | ``yes``           | | Enable time-varying thermodynamic pressure     |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+

.. _tab:fieldparams:

.. table:: ``COMMON`` keys for all field varaibles in the ``.par`` file

   +-------------------------+-----------------+-------------------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                         |
   +=========================+=================+=======================================================+
   | ``residualTol``         | | ``<real>``    | | Residual tolerance used by solver (not for CVODE)   | 
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``residualProj``        | | ``yes``       | | Controls the residual projection                    |
   |                         | | ``no``        |                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``writeToFieldFile``    | | ``yes``       | | Controls if fields will be written on output        |
   |                         | | ``no``        |                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+

.. _tab:meshparams:

.. table:: ``MESH`` keys in the ``.par`` file

   +-------------------------+-----------------+-------------------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                         |
   +=========================+=================+=======================================================+
   | ``motion``              | | ``none``      | | Mesh motion solver                                  |
   |                         | | ``user``      |                                                       |
   |                         | | ``elasticity``|                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``viscosity``           | | ``<real>``    | | Diffusivity for elasticity solver                   |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``numberOfBCFields``    | | ``<real>``    | | Number of field variables which have a boundary     |
   |                         |                 | |  condition in ``.re2`` file                         |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``firstBCFieldIndex``   | | ``<real>``    | | Field index of the first BC specified in ``.re2``   |
   |                         |                 | |  file                                               |
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

   +-------------------------+----------------+-----------------------------------------------+
   |   Key                   | | Value(s)     | | Description                                 |
   +=========================+================+===============================================+
   | ``preconditioner``      | | ``semg_xxt`` | | Preconditioning method                      |
   |                         | | ``semg_amg`` | | First time usage of AMG will write three    |
   |                         |                | | dump files to disc. Subsequently please run |
   |                         |                | | the amg_hypre tool to create the setup files|
   |                         |                | | required for the AMG solver initialization  |
   +-------------------------+----------------+-----------------------------------------------+

.. _tab:fieldparams:

.. table:: ``COMMON`` keys for temperature and scalar fields in the ``.par`` file

   +-------------------------+--------------+--------------------------------------------+
   |   Key                   | | Value(s)   | | Description                              |
   +=========================+==============+============================================+
   | ``solver``              | | ``helm``   | | Scalar solver                            | 
   |                         | | ``cvode``  |                                            |  
   |                         | | ``none``   |                                            |
   +-------------------------+--------------+--------------------------------------------+
   | ``advection``           | | ``yes``    | | Controls if advection is present         |
   |                         | | ``no``     |                                            |
   +-------------------------+--------------+--------------------------------------------+
   | ``absoluteTol``         | | ``<real>`` | | Absolute tolerance used by CVODE         |
   +-------------------------+--------------+--------------------------------------------+

.. _tab:temperatureparams:

.. table:: ``TEMPERATURE`` keys in the ``.par`` file

   +--------------------------+--------------+----------------------------------------------+
   |   Key                    | | Value(s)   | | Description                                |
   +==========================+==============+==============================================+
   |``ConjugatedHeatTransfer``| | ``yes``    | | Controls conjugate heat transfer           |
   |                          | | ``no``     |                                              |
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
   | ``stiff``                | | ``yes``      | | Controls if BDF or Adams Moulton is used   |
   |                          | | ``no``       |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``preconditioner``       | | ``user``     | | Preconditioner method                      |
   |                          | | ``none``     |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``dtMax``                | | ``<real>``   | | Maximum internal step size                 |
   |                          |                | | Controls splitting error of velocity       |
   |                          |                | | scalar coupling (e.g. set to 1-4 dt)       |
   +--------------------------+----------------+----------------------------------------------+


.. _user_files_re2:
-----------------------------------
re2 File
-----------------------------------

Stores the mesh and boundary condition. 
TODO: Add more details


.. _user_files_usr:

----------------------
usr File
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
SIZE File
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

-----------------------------------
ma2 File
-----------------------------------

TODO: Add more details

...................
Memory Requirements
...................

The memory footprint of a Nek5000 run depends on many factors such as 

- Local problem size per MPI rank (defined by ``lx1`` and ``lelt``)
- Total number of elements ``lelg``
- Number of MPI ranks
- Size of residual projection space ``mxprev`` 
- Dimension of Krylov subspace ``lgmres``.   
- Solvers (XXT, AMG, Helmholz, CVODE, ...)
- Number of scalars ``ldimt``
- Allocation policy (static, first touch, etc.) 

An a-priori estimate is hard to find and the output of the GNU SIZE utility does not
take into account the dynamic memory alloation of MPI, gslib, CVODE, etc. 
However, the true maximum memory usage is printed whenever Nek5000 exits.


