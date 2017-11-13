.. _user_files:

==========
User Files
==========

Each simulation is defined by three files: the .rea file, the ``.usr`` file, and the SIZE file.  In
addition, there is a derived ``.map`` file that is generated from the ``.rea`` file by running ``genmap``,
which will determine how the elements will be split across processors in the case of a parallel
run.  SIZE controls (at compile time) the polynomial degree used in the simulation, as well as the
space dimension :math:`d=2` or :math:`3`.

The SESSION.NAME file provides the name and path of the .rea file and the path to it.  It does not
however need to correspond to a ``.usr`` file of an identical name. This allows for different test
cases (``.usr`` files) that use the same geometry and boundary conditions (``.rea`` files).

This chapter provides an introduction to the basic files required to set up a Nek5000 simulation.

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

.. _user_files_usr:


-----------------------------------
Parameters File (.par)
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
   key = value

   [SECTION]
   key = value
   key = value

The sections are:

* ``GENERAL``
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
   | ``StartFrom``           | | ``<string>``| | Absolute/relative path of a field file     |
   |                         |               | | to restart the simulation from             |
   +-------------------------+---------------+----------------------------------------------+
   | ``StopAt``              | | ``endTime`` | | Type of total runtime assignment           |
   |                         | | ``numStep`` |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``endTime``             | | ``<real>``  | | Total runtime if ``StopAt=endTime``        |
   +-------------------------+---------------+----------------------------------------------+
   | ``numStep``             | | ``<real>``  | | Total number of steps if ``StopAt=numStep``|
   +-------------------------+---------------+----------------------------------------------+
   | ``VariableDT``          | | ``yes``     | | Enable variable time step (default is no)  |
   |                         | | ``no``      |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``DT``                  | | ``<real>``  | | - If ``VariableDT=yes`` denotes upper limit|
   |                         |               |   of time step                               | 
   |                         |               | | - If ``VariableDT=no`` denotes time step   |
   +-------------------------+---------------+----------------------------------------------+
   | ``WriteControl``        | | ``runTime`` | | Specify whether checkpointing is based on  |
   |                         | | ``timeStep``| | number of time steps or time intervals     |
   +-------------------------+---------------+----------------------------------------------+
   | ``WriteInterval``       | | ``<real>``  | | - If ``WriteControl=timeStep`` denotes the | 
   |                         |               | |   checkpoint frequency in number of time   | 
   |                         |               |     steps                                    |
   |                         |               | | - If ``WriteControl=runTime`` denotes the  |
   |                         |               | |   checkpoint frequency in runtime          |   
   +-------------------------+---------------+----------------------------------------------+
   | ``TargetCFL``           | | ``<real>``  | | If ``VariableDT=yes`` denotes the target   |
   |                         |               | | CFL number                                 |  
   +-------------------------+---------------+----------------------------------------------+
   | ``Filtering``           | | ``explicit``| | Specify filtering method                   | 
   |                         | | ``hpfrt``   |                                              | 
   +-------------------------+---------------+----------------------------------------------+
   | ``filterCutoffRatio``   | | ``<real>``  | | Fraction of unfiltered modes (0 to 1)      |
   +-------------------------+---------------+----------------------------------------------+
   | ``filterWeight``        | | ``<real>``  | | Filter weight                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``writeDoublePrecision``| | ``yes``     | | Output files in double precision           |
   |                         | | ``no``      |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``writeNFiles``         | | ``<real>``  | | Number of output files (default is 1)      |  
   +-------------------------+---------------+----------------------------------------------+
   | ``Dealiasing``          | | ``yes``     | | Enable/diasble over-integration            |
   |                         | | ``no``      |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``TimeStepper``         | | ``BDF2``    | | Time integration order                     |
   |                         | | ``BDF3``    |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``extrapolation``       | | ``standard``| | Extrapolation method                       |
   |                         | | ``OIFS``    |                                              |
   +-------------------------+---------------+----------------------------------------------+
   | ``OptLevel``            | | ``<real>``  | | Optimization level (1-3,default is 1)      |
   +-------------------------+---------------+----------------------------------------------+
   | ``LogLevel``            | | ``<real>``  | | Logfile verbosity level (1-2, default is 1)|
   +-------------------------+---------------+----------------------------------------------+
   | ``UserParam%%``         | | ``<real>``  | | User parameter (can be accessed through    |
   |                         |               | | uparam(%) array in ``.usr``                |
   +-------------------------+---------------+----------------------------------------------+



.. _tab:probtypeparams:

.. table:: ``PROBLEMTYPE`` keys in the ``.par`` file

   +---------------------------+---------------------+--------------------------------------------------+
   |   Key                     | | Value(s)          | | Description                                    |
   +===========================+=====================+==================================================+
   | ``equation``              | | ``incompNS``      | | Specify equation to solve:                     |
   |                           | | ``lowMachNS``     | | ``incompNS`` incompressible NS                 |
   |                           | | ``steadyStokes``  | | ``lowMachNS`` low-Mach NS                      |
   |                           | | ``incompLinNS``   | | ``steadyStokes`` steady stokes                 |
   |                           | | ``incompLinAdjNS``| | ``incompLinNS`` incompressible linearized NS   |
   |                           | | ``incompMHD``     | | ``incompLinAdjNS`` incompressible linearized   |
   |                           | | ``compNS``        |    adjoint NS                                    |
   |                           |                     | | ``incompMHD`` incompressible MHD               |
   |                           |                     | | ``compNS``  compressible NS                    |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``Axisymmetry``           | | ``yes``           | | Axisymmetric problem                           |
   |                           | | ``no``            |   (sets ``IFAXIS=.true.``)                       |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``Swirl``                 | | ``yes``           | | Enable axisymmetric azimuthal velocity         |
   |                           | | ``no``            | | component (sets ``IFAZIV=.true.``, stored      |
   |                           |                     | | in ``t(,,,,1)``)                               |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``CyclicBoundaries``      | | ``yes``           | | Cyclic periodic problem                        | 
   |                           | | ``no``            |   (sets ``IFCYCLIC=.true.``)                     |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``numberOfPerturbations`` | | ``<real>``        | | Number of perturbations if                     |
   |                           |                     | | ``equation=incompLinNS`` or                    |
   |                           |                     | | ``equation=incompLinAdjNS``                    |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``SolveBaseFlow``         | | ``yes``           | | Solve for base flow if                         |
   |                           | | ``no``            | | ``equation=incompLinNS`` or                    |
   |                           |                     | | ``equation=incompLinAdjNS``                    |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``VariableProperties``    | | ``yes``           | | Enable variable transport properties           |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``StressFormulation``     | | ``yes``           | | Enable stress formulation                      |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+
   | ``dp0dt``                 | | ``yes``           | | Enable time-varying thermodynamic pressure     |
   |                           | | ``no``            |                                                  |
   +---------------------------+---------------------+--------------------------------------------------+


.. _tab:meshparams:

.. table:: ``MESH`` keys in the ``.par`` file

   +-------------------------+-----------------+-------------------------------------------------------+
   |   Key                   | | Value(s)      | | Description                                         |
   +=========================+=================+=======================================================+
   | ``Motion``              | | ``none``      | | Enable mesh motion.                                 |
   |                         | | ``user``      | | ``user``: user-specified mesh velocity              |
   |                         | | ``elasticity``| | ``elasticity``: elasticity solver                   |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``Viscosity``           | | ``<real>``    | | Mesh solver diffusivity if ``Motion=elasticity``    |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``residualTol``         | | ``<real>``    | | Mesh solver residual tolerance if                   |
   |                         |                 |   ``Motion=elasticity``                               |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``residualProj``        | | ``yes``       | | Enable mesh solver residual projection if           |
   |                         | | ``no``        |   ``Motion=elasticity``                               |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``WriteToFieldFile``    | | ``yes``       | | Write mesh in field file                            |
   |                         | | ``no``        |                                                       |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``NumberOfBCFields``    | | ``<real>``    | | Number of BC fields in ``.re2`` file                |
   +-------------------------+-----------------+-------------------------------------------------------+
   | ``firstBCFieldIndex``   | | ``<real>``    | | Field index of the first BC specified in ``.re2``   |
   |                         |                 |   file                                                |
   +-------------------------+-----------------+-------------------------------------------------------+




.. _tab:velocityparams:

.. table:: ``VELOCITY`` keys in the ``.par`` file

   +-------------------------+--------------+------------------------------------------------+
   |   Key                   | | Value(s)   | | Description                                  |
   +=========================+==============+================================================+
   | ``ResidualTol``         | | ``<real>`` | | Residual tolerance                           | 
   +-------------------------+--------------+------------------------------------------------+
   | ``ResidualProj``        | | ``yes``    | | Enable residual projection                   |
   |                         | | ``no``     |                                                |
   +-------------------------+--------------+------------------------------------------------+
   | ``WriteToFieldFile``    | | ``yes``    | | Write to field file                          |
   |                         | | ``no``     |                                                |
   +-------------------------+--------------+------------------------------------------------+
   | ``Advection``           | | ``yes``    | | Enable advection                             |
   |                         | | ``no``     |                                                |
   +-------------------------+--------------+------------------------------------------------+
   | ``Viscosity``           | | ``<real>`` | | Positive value denotes dynamic viscosity,    |
   |                         |              | | negative value denotes Reynolds number       |
   |                         |              | | (required only if ``VariableProperties=no``  |
   +-------------------------+--------------+------------------------------------------------+
   | ``Density``             | | ``<real>`` | | Density                                      |
   |                         |              | | (required only if ``VariableProperties=no``) |
   +-------------------------+--------------+------------------------------------------------+



.. _tab:pressureparams:

.. table:: ``PRESSURE`` keys in the ``.par`` file

   +-------------------------+----------------+----------------------------------------------+
   |   Key                   | | Value(s)     | | Description                                |
   +=========================+================+==============================================+
   | ``Preconditioner``      | | ``semg_amg`` | | Preconditioner                             |
   |                         | | ``semg_xxt`` |                                              |
   +-------------------------+----------------+----------------------------------------------+
   | ``ResidualTol``         | | ``<real>``   | | Residual tolerance                         |
   +-------------------------+----------------+----------------------------------------------+
   | ``ResidualProj``        | | ``yes``      | | Enable residual projection                 |
   |                         | | ``no``       |                                              |
   +-------------------------+----------------+----------------------------------------------+
   | ``WriteToFieldFile``    | | ``yes``      | | Write to field file                        |
   |                         | | ``no``       |                                              |
   +-------------------------+----------------+----------------------------------------------+



.. _tab:temperatureparams:

.. table:: ``TEMPERATURE`` keys in the ``.par`` file

   +--------------------------+--------------+----------------------------------------------+
   |   Key                    | | Value(s)   | | Description                                |
   +==========================+==============+==============================================+
   | ``ResidualTol``          | | ``<real>`` | | Residual tolerance                         |
   +--------------------------+--------------+----------------------------------------------+
   | ``ResidualProj``         | | ``yes``    | | Enable residual projection                 |
   |                          | | ``no``     |                                              |
   +--------------------------+--------------+----------------------------------------------+
   |``ConjugatedHeatTransfer``| | ``yes``    | | Enable conjugate heat transfer             |
   |                          | | ``no``     |                                              |
   +--------------------------+--------------+----------------------------------------------+
   | ``WriteToFieldFile``     | | ``yes``    | | Write to field file                        |
   |                          | | ``no``     |                                              |
   +--------------------------+--------------+----------------------------------------------+
   | ``Advection``            | | ``yes``    | | Enable advection (default is yes)          |
   |                          | | ``no``     |                                              |
   +--------------------------+--------------+----------------------------------------------+
   | ``Conductivity``         | | ``<real>`` | | Thermal conductivity                       |
   |                          |              | | (required only if                          |
   |                          |              |   ``VariableProperties=no``)                 |
   +--------------------------+--------------+----------------------------------------------+
   | ``RhoCp``                | | ``<real>`` | | Rho*cp                                     |
   |                          |              | | (required only if                          |
   |                          |              |   ``VariableProperties=no``)                 |
   +--------------------------+--------------+----------------------------------------------+



.. _tab:scalarparams:

.. table:: ``SCALAR%%`` keys in the ``.par`` file

   +--------------------------+----------------+----------------------------------------------+
   |   Key                    | | Value(s)     | | Description                                |
   +==========================+================+==============================================+
   | ``Solver``               | | ``helm``     | | Specify solver (Helmholtz, CVODE, or none) | 
   |                          | | ``cvode``    |                                              |  
   |                          | | ``none``     |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``Advection``            | | ``yes``      | | Enable advection (default is yes)          |
   |                          | | ``no``       |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``ResidualTol``          | | ``<real>``   | | Residual tolerance if ``Solver=helm``      |
   +--------------------------+----------------+----------------------------------------------+
   |``ConjugatedHeatTransfer``| | ``yes``      | | Enable conjugate heat transfer             |
   |                          | | ``no``       |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``Density``              | | ``<real>``   | | Density (required only if                  |
   |                          |                |   ``VariableProperties=no``)                 |
   +--------------------------+----------------+----------------------------------------------+
   | ``Diffusivity``          | | ``<real>``   | | Diffusivity (required only if              | 
   |                          |                |   ``VariableProperties=no``)                 |
   +--------------------------+----------------+----------------------------------------------+
   | ``WriteToFieldFile``     | | ``yes``      | | Write to field file                        |
   |                          | | ``no``       |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``AbsoluteTol``          | | ``<real>>``  | | Absolute tolerance if ``Solver=cvode``     |
   +--------------------------+----------------+----------------------------------------------+



.. _tab:cvodeparams:

.. table:: ``CVODE`` keys in the ``.par`` file

   +--------------------------+----------------+----------------------------------------------+
   |   Key                    | | Value(s)     | | Description                                |
   +==========================+================+==============================================+
   | ``RelativeTol``          | | ``<real>``   | | Relative tolerance (applies to all scalars)|
   +--------------------------+----------------+----------------------------------------------+
   | ``Stiff``                | | ``yes``      | | If ``Stiff=yes`` use BDF timestepper,      |
   |                          | | ``no``       | | if ``Stiff=no`` use Adams Moulton,         |
   |                          |                | | default is no                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``PreConditioner``       | | ``yes``      | | Enable user-supplied preconditioner        |
   |                          | | ``no``       |                                              |
   +--------------------------+----------------+----------------------------------------------+
   | ``DtMax``                | | ``<real>``   | | Maximum CVODE timestep allowed             |
   +--------------------------+----------------+----------------------------------------------+


----------------------
Case Setup File (.usr)
----------------------

.....................
Contents of .usr File
.....................


The most important interface to Nek5000 is the set of Fortran subroutines that are contained in the
``.usr`` file.  This file allows direct access to all runtime variables.  Here, the user may
specify spatially varying properties (e.g., viscosity), volumetric heating sources, body forces,
and so forth.  One can also specify arbitrary initial and boundary conditions through the routines
``useric`` and ``userbc``.  The routine ``userchk`` allows the user to interrogate the
solution at the end of each timestep for diagnostic purposes.   The ``.usr`` files provided in
the ``Nek5000/short_tests/`` directory illustrate several of the more common analysis tools.  For
instance, there are utilities for computing the time average of :math:`u`, :math:`u^2`, etc. so that one
can analyze mean and rms distributions with the postprocessor.  There are routines for computing
the vorticity or the scalar :math:`\lambda_2` for vortex identification, and so forth.

.....................
Routines in .usr File
.....................



The routine ``uservp`` specifies the variable properties of the governing equations.  This
routine is called once per processor, and once per discrete point therein. 


+---------------------------------+----------------------+--------------+-------------+
| Equation                        | ``utrans``           | ``udiff``    | ``ifield``  |
+=================================+======================+==============+=============+
| Momentum Eq. :eq:`ns_momentum`  | :math:`\rho`         | :math:`\mu`  | 1           |
+---------------------------------+----------------------+--------------+-------------+
| Energy Eq. :eq:`energy`         | :math:`\rho c_p`     | :math:`k`    | 2           |
+---------------------------------+----------------------+--------------+-------------+
| Passive scalar :eq:`pass_scal`  | :math:`(\rho c_p)_i` | :math:`k_i`  | :math:`i-1` |
+---------------------------------+----------------------+--------------+-------------+

.. code-block:: fortran
 
 subroutine uservp (ix,iy,iz,eg)
 include 'SIZE'
 include 'TOTAL'
 include 'NEKUSE'

 integer iel
 iel = gllel(eg)

 udiff =0.
 utrans=0.

 return
 end

The routine ``userdat`` is called right after the geometry is loaded into NEK5000 and prior to
the distribution of the GLL points. This routine is called once per processor but for all the data
on that processor. At this stage the elements can be modified as long as the topology is preserved.
It is also possible to alter the type of boundary condition that is initially attributed in the
``.rea`` file, as illustrated below (the array ``cbc(face,iel,field``) contains the boundary
conditions per face and field of each element). Note the spacing allocated to each BC string is of
three units.

.. code-block:: fortran

  subroutine usrdat
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'
  integer iel,f

  do iel=1,nelt  !  Force flux BCs
  do f=1,2*ndim
     if (cbc(f,iel,1).eq.'W  ') cbc(f,iel,2) = 'f  ' ! flux BC for temperature
  end do
  end do

  return
  end

The routine ``usrdat2`` is called after the GLL points were distributed and allows at this point only for affine transformations of the geometry.

.. code-block:: fortran

  subroutine usrdat2
  include 'SIZE'
  include 'TOTAL'

  return
  end

The routine ``userf`` is called once for each point and provides the force term in Eq. :eq:`ns_momentum`. Not that according to the dimensionalization in Eq. :eq:`ns_momentum` the force term :math:`\mathbf{f}` is in fact multiplied by the density :math:`\rho`.

.. code-block:: fortran

  subroutine userf  (ix,iy,iz,eg)
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'

  ffx = 0.0
  ffy = 0.0
  ffz = 0.0

  return
  end

Similarly to ``userf`` the routine ``userq`` provides the force term in Eq. :eq:`energy` and the subsequent passive scalar equations according to Eq. :eq:`pass_scal`.

.. code-block:: fortran

  subroutine userq  (ix,iy,iz,eg)
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'

  qvol   = 0.

  return
  end

The boundary conditions are assigned in ``userbc`` for both the fluid, temperature and all other scalars. An extensive list of such possible boundary conditions is available in :ref:`sec:boundary`. 

.. code-block:: fortran

  subroutine userbc (ix,iy,iz,iside,ieg)
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'

  ux=0.0
  uy=0.0
  uz=0.0
  temp=0.0
  flux = 1.0

  return
  end

Initial conditions are attributed in ``useric`` similarly to the boundary conditions

.. code-block:: fortran

  subroutine useric (ix,iy,iz,ieg)
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'

  uy=0.0
  ux=0.0
  uz=1.0

  return
  end

The routine ``userchk`` is called once per processor after each timestep (and once after the initialization is finished). This is the section where the solution can be interrogated and subsequent changes can be made.

.. code-block:: fortran

  subroutine userchk
  include 'SIZE'
  include 'TOTAL'
  include 'NEKUSE'

  call outpost(vx,vy,vz,pr,t,'ext')

  return
  end

The routine ``usrdat3`` is not widely used, however it shares the same properties with ``usrdat2``.

.. code-block:: fortran

        subroutine usrdat3
        include 'SIZE'
        include 'TOTAL'
  c
        return
        end

Nek5000 can solve the dimensional or non-dimensional equations by setting the following parameters

+---------------------------+-------------------------------------+
| Dimensional parameters    | Non-dimensional parameters          |
+===========================+=====================================+
| ``p1`` = :math:`\rho`     | ``p1`` = 1                          |
+---------------------------+-------------------------------------+
| ``p2`` = :math:`\nu`      | ``p2`` = :math:`1/Re` :math:`(-Re)` |
+---------------------------+-------------------------------------+
| ``p7`` = :math:`\rho C_p` | ``p7`` = 1                          |
+---------------------------+-------------------------------------+
| ``p8`` = :math:`k`        | ``p8`` = :math:`1/Pe` :math:`(-Pe)` |
+---------------------------+-------------------------------------+

alternatively the variable properties can be set in the ``uservp`` routine.

**What is a SESSION file?**

To run Nek5000, each simulation must have a SESSION.NAME file. This file is read in by the code and gives the path to the relevant files describing the structure and parameters of the simulation. The SESSION.NAME file is a file that contains the name of the simulation and the full path to supporting files. For example, to run the eddy example from the repository, the SESSION.NAME file would look like

.. code-block:: none

  eddy_uv
  /homes/user_name/nek5_svn/examples/eddy/


------------------------
Problem-Size File (SIZE)
------------------------

SIZE file defines the problem size, i.e. spatial points at which the solution is to be evaluated within each element, number of elements per processor etc.
The SIZE file governs the memory allocation for most of the arrays
in Nek5000, with the exception of those required by the C utilities.
The primary parameters of interest in SIZE are:

* **ldim** = 2 or 3.  This must be set to 2 for two-dimensional or axisymmetric simulations  (the latter only partially supported) or to 3 for three-dimensional simulations.
* **lx1** controls the polynomial order of the approximation, :math:`N = {\tt lx1-1}`.
* **lxd** controls the polynomial order of the integration forconvective terms.  Generally, :math:`{\tt lxd=3 * lx1/2}`.  On some platforms, however,it is important for memory access performance that ``lx1`` and ``lxd`` be even.
* **lx2** = ``lx1`` or ``lx1-2``.  This determines the formulation for the Navier-Stokes  solver (i.e., the choice between the :math:`\mathbb{P}_N - \mathbb{P}_N` or :math:`\mathbb{P}_N - \mathbb{P}_{N-2}` methods) and the approximation order for the pressure, ``lx2-1``.
* **lelt** determines the *maximum* number of elements *per processor*

The total size of the problem is ``lx1*ly1*lz1*lelt``.

...................
Memory Requirements
...................

Per-processor memory requirements for  Nek5000 scale
roughly as 400 8-byte words per allocated gridpoint.  The number
of *allocated* gridpoints per processor is
:math:`n_{\max}` = ``lx1*ly1*lz1*lelt``.
(For 3D, ``lz1=ly1=lx1``; for 2D, ``lz1=1``, ``ly1=lx1``.)
If required for a particular simulation, more memory may be made
available by using additional processors.  For example, suppose
one needed to run a simulation with 6000 elements of order :math:`N=9`.
To leading order, the total memory requirements would be
:math:`{\tt \approx E(N+1)^3 points \times 400 (wds/pt) \times 8 bytes/wd =
6000 \times 10^3 \times 400 \times 8 = 19.2}` GB.  Assuming there
is 400 MB of memory per core available to the user (after accounting
for OS requirements), then one could run this simulation with
:math:`{\tt P \geq 19,200 MB / (400 MB/proc) = 48}` processors.
To do so, it would be necessary to set :math:`{\tt lelt} \geq 6000/48 = 125`.

We note two other parameters of interest in the parallel context:

* **lp**, the maximum number of processors that can be used.
* **lelg**, an upper bound on the number of elements in the simulation.

There is a slight memory penalty associated with these variables, so
one generally does not want to have them excessively large.  It is
common, however, to have lp be as large as anticipated for a given
case so that the executable can be run without recompiling on
any admissible number of processors (:math:`P_{mem} \leq P \leq E`,
where :math:`P_{mem}` is the value computed above).

-----------------------------------
Geometry and Parameters File (.rea)
-----------------------------------

The ``.rea`` file consists of several sections. The mesh specifications  with **geometry**, **curvature** and **boundary conditions** are in the second section.

...............................
Parameters and logical switches
...............................

**parameters** 
    These control the runtime parameters such as viscosity,
    conductivity, number of steps, timestep size, order of the timestepping,
    frequency of output, iteration tolerances, flow rate, filter strength,
    etc.   There are also a number of free parameters that the user can
    use as handles to be passed into the user defined routines in the ``.usr`` file.
**passive scalar data** 
    This information can be specified also in the ``uservp`` routine in the ``.usr``
    file. If specified in the ``.rea`` file then the coefficients for the conductivity 
    term are listed in ascending order for passive scalars ranging ``1..9`` 
    followed by the values for the :math:`\rho c_p` coefficients.

    .. code-block:: none

      4  Lines of passive scalar data follows 2 CONDUCT; 2 RHOCP
         1.00000       1.00000       1.00000       1.00000       1.00000
         1.00000       1.00000       1.00000       1.00000
         1.00000       1.00000       1.00000       1.00000       1.00000
         1.00000       1.00000       1.00000       1.00000

**logicals**  
    These determine whether one is computing a steady or unsteady
    solution, whether advection is turned on, etc.


Next we have the logical switches as follow, a detailed explanation to be found in :ref:`sec:switches` 

.. code-block:: none


           13  LOGICAL SWITCHES FOLLOW
  T     IFFLOW
  T     IFHEAT
  T     IFTRAN
  T T F F F F F F F F F IFNAV & IFADVC (convection in P.S. fields)
  F F T T T T T T T T T T IFTMSH (IF mesh for this field is T mesh)
  F     IFAXIS
  F     IFSTRS
  F     IFSPLIT
  F     IFMGRID
  F     IFMODEL
  F     IFKEPS
  F     IFMVBD
  F     IFCHAR

................................
Mesh and boundary condition info
................................

.. highlight:: none

**geometry**
    The geometry is specified in an arcane format specifying
    the :math:`xyz` locations of each of the eight points for each element,
    or the :math:`xy` locations of each of the four points for each element in 2D.
    A line of the following type may be encountered at the beginning 
    of the mesh section of the ``.rea`` file::

      3.33333       3.33333     -0.833333      -1.16667     XFAC,YFAC,XZERO,YZERO

    This part is to be read by Prenek and provides the origin of the system of 
    coordinates ``XZERO;YZERO`` as well as the size of the cartesian units 
    ``XFAC;YFAC``. This one line has no impact on the mesh as being read in Nek5000.

    The header of the mesh data may have the following representation::

       **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8
            226  3         192           NEL,NDIM,NELV

    The header states first how many elements are available in total (226), what
    dimension is the the problem (here three dimensional), and how many elements 
    are in the fluid mesh (192).

      .. _tab:element:

      .. table:: Geometry description in ``.rea`` file

         +-------------------------------------------------------------------------------------+
         | ``ELEMENT 1 [ 1A] GROUP 0``                                                         |
         +=====================================================================================+
         | ``Face {1,2,3,4}``                                                                  |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`x_{1,\ldots,4}=` | 0.000000E+00 | 0.171820E+00 | 0.146403E+00 | 0.000000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`y_{1,\ldots,4}=` | 0.190000E+00 | 0.168202E+00 | 0.343640E+00 | 0.380000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`z_{1,\ldots,4}=` | 0.000000E+00 | 0.000000E+00 | 0.000000E+00 | 0.000000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+
         | ``Face {5,6,7,8}``                                                                  |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`x_{5,\ldots,8}=` | 0.000000E+00 | 0.171820E+00 | 0.146403E+00 | 0.000000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`y_{5,\ldots,8}=` | 0.190000E+00 | 0.168202E+00 | 0.343640E+00 | 0.380000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+
         | :math:`z_{5,\ldots,8}=` | 0.250000E+00 | 0.250000E+00 | 0.250000E+00 | 0.250000E+00 |
         +-------------------------+--------------+--------------+--------------+--------------+

    Following the header, all elements are listed. The fluid elements are listed 
    first, followed by all solid elements if present. In this case there are (34) 
    solid elements.

    The data following the header is formatted as shown in :numref:`tab:element`. This provides all the coordinates of an element for top and bottom faces. The numbering of the vertices is shown in Fig. :numref:`fig:elorder`. The header for each element as in :numref:`tab:element`, i.e. ``[1A] GROUP`` is reminiscent of older Nek5000 format and does not impact the mesh generation at this stage. (We are inquiring whether other groups still use it.)

      .. _fig:elorder:

      .. figure:: figs/3dcube_1.png
          :align: center
          :figclass: align-center
          :alt: rea-geometry

          Geometry description in ``.rea`` file (sketch of one element ordering - Preprocessor 
          corner notation) 


**curvature**
    This section describes the curvature of the elements. It is expressed as deformation of the linear elements.
    Therefore, if no elements are curved (if only linear elements are present) the section remains empty.

    The section header may look like this::

      640 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE

    Curvature information is provided by edge and element. Therefore up to 12 curvature entries can be present for each element.
    Only non-trivial curvature data needs to be provided, i.e., edges that correspond to linear elements, since they have no curvature, will have no entry.
    The formatting for the curvature data is provided in :numref:`tab:midside`.

      .. _tab:midside:

      .. table:: Curvature information specification

         +-----------+---------+--------------+--------------+--------------+--------------+--------------+------------+
         | ``IEDGE`` | ``IEL`` | ``CURVE(1)`` | ``CURVE(2)`` | ``CURVE(3)`` | ``CURVE(4)`` | ``CURVE(5)`` | ``CCURVE`` |
         +===========+=========+==============+==============+==============+==============+==============+============+
         | 9         | 2       | 0.125713     | -0.992067    | 0.00000      | 0.00000      | 0.00000      | m          |
         +-----------+---------+--------------+--------------+--------------+--------------+--------------+------------+
         | 10        | 38      | 0.125713     | -0.992067    | 3.00000      | 0.00000      | 0.00000      | m          |
         +-----------+---------+--------------+--------------+--------------+--------------+--------------+------------+
         | 1         | 40      | 1.00000      | 0.000000     | 0.00000      | 0.00000      | 0.00000      | C          |
         +-----------+---------+--------------+--------------+--------------+--------------+--------------+------------+

    There are several types of possible curvature information represented by ``CCURVE``. This include:

    - 'C' stands for circle and is given by the radius of the circle,  in ``CURVE(1)``, all other compoentns of the ``CURVE`` array are not used but need to be present.
    - 's' stands for sphere and is given by the radius and the center of the sphere, thus filling the first 4 components of the ``CURVE`` array. The fifth component needs to be present but is not utilized.
    - 'm' is given by the coordinates of the midside-node, thus using the first 3 components of the ``CURVE`` array, and leads to a second order reconstruction of the face.  The fourth and fifth components need to be present but are not utilized.

    Both 'C' and 's' types allow for a surface of as high order as the polynomial used in the spectral method, since they have an underlying analytical description, any circle arc can be fully determined by the radius and end points. However for the 'm' curved element descriptor the surface can be reconstructed only up to second order. This can be later updated to match the high-order polynomial after the GLL points have been distributed across the boundaries. This is the only general mean to describe curvature currrently in Nek5000 and corresponds to a HEX20 representation.

      .. _fig:edges:

      .. figure:: figs/3dcube.png
          :align: center
          :figclass: align-center
          :alt: edge-numbering

          Edge numbering in ``.rea`` file, the edge number is in between parenthesis. The other
          numbers represent vertices.

    For better understanding let us focus on what the data in :numref:`tab:midside` signifies. Edge 9 of element 2 has a edge  midpoint at (0.125713, -0.992067, 0.00000)  and so on. For edge numbering the reader is advised to check Fig. :numref:`fig:edges`, which illustrates the relationship between vertex numbering and edge numbering.

    To maninpulate the geometry in Nek5000 at runtime, it is possible to use  ``usrdat2``. In this subroutine the user can deform the geometry to match the intended surface, followed by a call to the subroutine ``fixgeom`` which can realign the point distribution in the interior of the element.

      .. _fig:ex1:

      .. figure:: figs/base1.png
          :align: center
          :figclass: align-center
          :alt: edge-numbering

          Example mesh - without curvature. Square dots represent example vertices.

    We also note, that, unlike the geometry data, each curvature entry (as shown in :numref:`tab:midside`) is formatted and the format is **dependent on the total number of elements**. Three cases exist as shown in the code below:

      .. code-block:: fortranfixed

                       if (nelgt.lt.1000) then
                          write(10,'(i3,i3,5g14.6,1x,a1)') i,eg,
       $                  (vcurve(k,i,kb),k=1,5),cc
                       else if (nelgt.lt.1000000) then
                          write(10,'(i2,i6,5g14.6,1x,a1)') i,eg,
       $                  (vcurve(k,i,kb),k=1,5),cc
                       else
                          write(10,'(i2,i12,5g14.6,1x,a1)') i,eg,
       $                  (vcurve(k,i,kb),k=1,5),cc

    The fortran format is as follows:

    - For a total number of elements below 1,000 the format is ``(i3,i3,5g14.6,1x,a1)``.
    - For a total number of elements 1,000 - 999,999 the format is ``(i2,i6,5g14.6,1x,a1)``.
    - For a total number of elements above 999,999 the format is ``(i2,i12,5g14.6,1x,a1)``.

    .. _fig:ex2:

    .. figure:: figs/modified1.png
        :align: center
        :figclass: align-center
        :alt: edge-numbering

        Example mesh - with curvature. Circular dots represent example midsize points.

    To further illustrate the usage of curvature data, let us examine an example of ``.rea`` file with and wiuthout curvature information and the corresponding mesh representation. :numref:`fig:ex1` represents a 12 element box mesh (2x2x3, with periodic conditions in :math:`z`) without curvature, while :numref:`fig:ex2` presents the same mesh with a sinusoidal deformation in direction :math:`y`. Only two edges per element are curved.

    The input for the mesh without curvature is:

    .. include:: mesh_example.txt
        :literal:

    The input for the mesh with curvature is:

    .. include:: mesh_curv_example.txt
        :literal:

    Note that element and boundary condition information are identical between the two cases.

**boundary conditions**
    Boundary conditions (BCs) are specified for each field in sequence: velocity, temperature and passive scalars. The section header for each field will be as follows (example for the velocity)::

      ***** FLUID   BOUNDARY CONDITIONS *****

    and the data is stored as illustarted in :numref:`tab:bcs`. For each field boundary conditions are listed for each face of each element.

    Boundary conditions are given in order per each element, see :numref:`tab:bcs` column ``IEL``, and faces listed in ascending order 1-6 in column ``IFACE``. Note that the header in :numref:`tab:bcs` does not appear in the actual ``.rea``.

    The ordering for faces each element is shown in :numref:`fig:forder`. A total equivalent to :math:`6N_{field}` boundary conditions are listed for each field, where :math:`N_{field}` is the number of elements for the specific field. :math:`N_{field}` is equal to the total number of fluid elements for the velocity and equal to the total number of elements (including solid elements) for temperature. For the passive scalars it will depend on the specific choice, but typically scalars are solved on the temeprature mesh (solid+fluid).

      .. _fig:forder:

      .. figure:: figs/3dcube_2.png
          :align: center
          :figclass: align-center
          :alt: edge-numbering

          Face ordering for each element.

    Each BC letter condition is formed by three characters. Common BCs include:

    - ``E`` - internal boundary condition. No additional information needs to be provided.
    - ``SYM`` - symmetry boundary condition. No additional information needs to be provided.
    - ``P`` - periodic boundary conditions,  which indicates that an element face is connected to another element to establish a periodic BC. The connecting element and face need be  to specified in ``CONN-IEL`` and ``CONN-IFACE``.
    - ``v`` - imposed velocity boundary conditions (inlet). The value is specified in the user subroutines. No additional information needs to be provided in the ``.rea`` file.
    - ``W`` - wall boundary condition (no-slip) for the velocity. No additional information needs to be provided.
    - ``O`` - outlet boundary condition (velocity). No additional information needs to be provided.
    - ``t`` - imposed temperature  boundary conditions (inlet). The value is specified in the user subroutines. No additional information needs to be provided in the ``.rea`` file.
    - ``f`` - imposed heat flux  boundary conditions (temperature). The value is specified in the user subroutines. No additional information needs to be provided in the ``.rea`` file.
    - ``I`` - adiabatic boundary conditions (temeperature). No additional information needs to be provided.

    Many of the BCs support either a constant specification or a user defined specification which may be an arbitrary function.   For example, a constant Dirichlet BC for velocity is specified by ``V``, while a user defined BC is specified by ``v``.   This upper/lower-case distinction is  used for all cases.   There are about 70 different types of boundary conditions in all, including free-surface, moving boundary, heat flux, convective cooling, etc. The above cases are just the most used types.

      .. _tab:bcs:

      .. table:: Formatting of boundary conditions input.

         +---------+---------+-----------+--------------+----------------+---------+---------+---------+
         | ``CBC`` | ``IEL`` | ``IFACE`` | ``CONN-IEL`` | ``CONN-IFACE`` |         |         |         |
         +=========+=========+===========+==============+================+=========+=========+=========+
         | E       | 1       | 1         | 4.00000      | 3.00000        | 0.00000 | 0.00000 | 0.00000 |
         +---------+---------+-----------+--------------+----------------+---------+---------+---------+
         | ``..``  | ``..``  | ``..``    | ``..``       | ``..``         | ``..``  | ``..``  | ``..``  |
         +---------+---------+-----------+--------------+----------------+---------+---------+---------+
         | W       | 5       | 3         | 0.00000      | 0.00000        | 0.00000 | 0.00000 | 0.00000 |
         +---------+---------+-----------+--------------+----------------+---------+---------+---------+
         | ``..``  | ``..``  | ``..``    | ``..``       | ``..``         | ``..``  | ``..``  | ``..``  |
         +---------+---------+-----------+--------------+----------------+---------+---------+---------+
         | P       | 23      | 5         | 149.000      | 6.00000        | 0.00000 | 0.00000 | 0.00000 |
         +---------+---------+-----------+--------------+----------------+---------+---------+---------+

    As in the case of the curvature entries, the boundary conditions entries are formatted and **the format is dependent on the total number of elements**.
    The code below shows an example of writing statement for boundary conditions:

      .. code-block:: fortranfixed

                        if (nlg.lt.1000) then
                           write(10,'(a1,a3,2i3,5g14.6)')
           $               chtemp,s3,eg,i,(vbc(ii,i,kb),ii=1,5)
                        else if (nlg.lt.100000) then
                           write(10,'(a1,a3,i5,i1,5g14.6)')
           $               chtemp,s3,eg,i,(vbc(ii,i,kb),ii=1,5)
                        else if (nlg.lt.1000000) then
                           write(10,'(a1,a3,i6,5g14.6)')
           $               chtemp,s3,eg,(vbc(ii,i,kb),ii=1,5)
                        else
                           write(10,'(a1,a3,i12,5g18.11)')
           $               chtemp,s3,eg,(vbc(ii,i,kb),ii=1,5)
                        end if

    The fortran format is as follows:

    - For a total number of elements below 1,000 the format is ``(a1,a3,2i3,5g14.6)``.
    - For a total number of elements 1,000 - 99,999 the format is ``(a1,a3,i5,i1,5g14.6)``.
    - For a total number of elements 100,000 - 999,999 the format is ``(a1,a3,i6,5g14.6)``.
    - For a total number of elements above 999,999 the format is ``(a1,a3,i12,5g18.11)``.

    We note that:

    - The first item in the format for each of the four cases is a string containing a space.
    - The second item in the format for each of the four cases is a string specifying the boundary condition type.
    - In cases where the total number of elements is bigger than 99,999, the ``IFACE`` item is omitted. Given that Nek5000 already knows the ordering of the actual faces within each element in column ``IFACE`` is in fact not needed.
    - The number of significant digits increases in the fourth case. This is needed for periodic boundary conditions.

...........
Output info
...........

**restart conditions**
    Here, one can specify a file to use as an initial condition.
    The initial condition need not be of the same polynomial order
    as the current simulation.   One can also specify that, for example,
    the velocity is to come from one file and the temperature from another.
    The initial time is taken from the last specified restart file, but
    this can be overridden.

**history points**
    The following section defines history points in the ``.rea`` file, see example ``vortex/r1854a.rea``, or ``shear4/shear4.rea``::

       0 PACKETS OF DATA FOLLOW
       ***** HISTORY AND INTEGRAL DATA *****
           56 POINTS. H code, I,J,H,IEL
       UVWP    H     31     31   1   6
       UVWP    H     31     31   31  6
       UVWP    H     31     31   31  54
        "      "      "      "    "   "

    The ``"56 POINTS"`` line needs to be followed by 56 lines of the type shown. However, in each of the following lines, which have the ``UVWP`` etc., location is CRUCIAL, it
    must be layed out exactly as indicated above (these lines contain character strings, they use formatted reads), it is therefore advisable to refer to the examples ``vortex, shear4``.  If you want to pick points close to the center of element 1 and are running with ``lx1=10``, say, you might choose ``UVWP H 5 5 5 1``. (the indicated point would really be at the middle of the element only if ``lx1=9``)

    The ``UVWP`` tells the code to write the 3 velocity components and pressure to the ``.sch`` file at
    each timestep (or, more precisely, whenever ``mod(istep,iohis)=0``, where ``iohis=param(52))``.
    Note that if you have more than one history point then they are written sequentially at each
    timestep. Thus 10 steps in the first example with ``param(52)=2`` would write ``(10/2)*56 = 280``
    lines to the ``.sch`` file, with 4 entries per line. The "H" indicates that the entry corresponds to a requested history point. A note of caution: if the ``ijk`` values (5 5 5 in the preceding example line) exceed ``lx1,ly1,lz1`` of your SIZE file, then they are truncated to that value. For example, if ``lx1=10`` for the data at the top (31 31 31) then the code will use ``ijk`` of (10 10 10), plus the given element number, in identifying the history point. It is often useful to set ``ijk`` to large values (i.e., > ``lx1``) because the endpoints of the spectral element mesh are invariant when ``lx1`` is changed.

**output specifications**
    Outputs are discussed in a separate section of the manual, available online.

It is important to note that Nek5000 currently supports two input file
formats, ASCII and binary.   The ``.rea`` file format
described above is ASCII.  For the binary format, all sections
of the ``.rea`` file having storage requirements that scale with
number of elements (i.e., geometry, curvature, and boundary
conditions) are moved to a second, ``.re2``, file and
written in binary.   The remaining sections continue to
reside in the ``.rea`` file.   The distinction between
the ASCII and binary formats is indicated in the ``.rea``
file by having a negative number of elements.
There are converters, ``reatore2`` and ``re2torea``, in the Nek5000
tools directory to change between formats.   The binary file
format is the default and important for ``I/O`` performance when the
number of elements is large ( :math:`>100000`, say).

..........
Parameters
..........

- :math:`\rho`, the density, is taken to be time-independent and
  constant; however, in a multi-fluid system
  different fluids can have different value of constant density.
- :math:`\mu`, the dynamic viscosity can vary arbitrarily in
  time and space; it can also be a function of temperature
  (if the energy equation is included) and strain rate
  invariants (if the stress formulation is selected).
- :math:`\sigma`, the surface-tension coefficient can vary
  arbitrarily in
  time and space; it can also be a function of temperature
  and passive scalars.
- :math:`\overline{\beta}`, the effective thermal expansion
  coefficient, is
  assumed time-independent and constant.
- :math:`{\bf f}(t)`, the body force per unit mass term can
  vary with time, space, temperature and passive scalars.
- :math:`\rho c_{p}`, the volumetric specific heat, can vary
  arbitrarily with time, space and temperature.
- :math:`\rho L`, the volumetric latent heat of fusion at a front,
  is taken to be time-independent and constant; however,
  different constants can be assigned to different fronts.
- :math:`k`, the thermal conductivity, can vary with time,
  space and temperature.
- :math:`q_{vol}`, the volumetric heat generation, can vary with
  time, space and temperature.
- :math:`h_{c}`, the convection heat transfer coefficient, can vary
  with time, space and temperature.
- :math:`h_{rad}`, the Stefan-Boltzmann constant/view-factor product,
  can vary with time, space and temperature.
- :math:`T_{\infty}`, the environmental temperature, can vary
  with time and space.
- :math:`T_{melt}`, the melting temperature at a front, is taken
  with time and space; however, different melting temperature
  can be assigned to different fronts.

In the solution of the governing equations together with
the boundary and initial conditions, Nek5000 treats the
above parameters as pure numerical values; their
physical significance depends on the user's choice of units.
The system of units used is arbitrary (MKS, English, CGS,
etc.). However, the system chosen must be used consistently
throughout. For instance, if the equations and geometry
have been non-dimensionalized, the :math:`\mu / \rho` in the fluid
momentum equation is in fact
the inverse Reynolds number, whereas if the equations are
dimensional, :math:`\mu / \rho` represents the kinematic viscosity with
dimensions of :math:`length^{2}/time`.

-----------
Data Layout
-----------

Nek5000 was designed with two principal performance criteria in mind,
namely, *single-node* performance and *parallel* performance.

A key precept in obtaining good single node performance was to use,
wherever possible, unit-stride memory addressing, which is realized by
using contiguously declared arrays and then accessing the data in
the correct order.   Data locality is thus central to good serial
performance.   To ensure that this performance is not compromised
in parallel, the parallel message-passing data model is used, in which
each processor has its own local (private) address space.  Parallel
data, therefore, is laid out just as in the serial case, save that there
are multiple copies of the arrays---one per processor, each containing
different data.  Unlike the shared memory model, this distributed memory
model makes data locality transparent and thus simplifies the task of
analyzing and optimizing parallel performance.

Some fundamentals of Nek5000's internal data layout are given below.

1. Data is laid out as  :math:`u_{ijk}^e = u(i,j,k,e)`

   .. |br| raw:: html

      <br />

   ``i=1,...,nx1``   (``nx1 = lx1``) |br|
   ``j=1,...,ny1``   (``ny1 = lx1``) |br|
   ``k=1,...,nz1``   (``nz1 = lx1`` or 1, according to ndim=3 or 2)

   ``e=1,...,nelv``, where ``nelv`` :math:`\leq` ``lelv``, and ``lelv`` is the upper
   bound on number of elements, *per processor*.
2. Fortran data is stored in column major order (opposite of C).
3. All data arrays are thus contiguous, even when :math:`{\tt nelv} < {\tt lelv}`.
4. Data accesses are thus primarily unit-stride (see chap.8 of DFM
   for importance of this point), and in particular, all data on
   a given processor can be accessed as, e.g.,

      .. code-block:: fortran

         do i=1,nx1*ny1*nz1*nelv
            u(i,1,1,1) = vx(i,1,1,1)
         end do

   which is equivalent but superior (WHY?) to:

      .. code-block:: fortran

         do e=1,nelv
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            u(i,j,k,e) = vx(i,j,k,e)
         end do
         end do
         end do
         end do

   which is equivalent but vastly superior (WHY?) to:

      .. code-block:: fortran

         do i=1,nx1
         do j=1,ny1
         do k=1,nz1
         do e=1,nelv
            u(i,j,k,e) = vx(i,j,k,e)
         end do
         end do
         end do
         end do
5. All data arrays are stored according to the SPMD programming
   model, in which address spaces that are local to each processor
   are private --- not accessible to other processors except through
   interprocessor data-transfer (i.e., message passing).  Thus

      .. code-block:: fortran

         do i=1,nx1*ny1*nz1*nelv
            u(i,1,1,1) = vx(i,1,1,1)
         end do

   means different things on different processors and ``nelv`` may
   differ from one processor to the next.
6. For the most part, low-level loops such as above are expressed in
   higher level routines only through subroutine calls, e.g.,:

      .. code-block:: fortran

         call copy(u,vx,n)

   where ``n:=nx1*ny1*nz1*nelv``.   Notable exceptions are in places where
   performance is critical, e.g., in the middle of certain iterative
   solvers. 
