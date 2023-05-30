.. _case_files_usr:

=========================
User Routines File (.usr)
=========================

The user file is a special case file that implements the user interface to *Nek5000*. 
What follows is a brief description of the available subroutines and a primer on writing your own user code.
An empty user file template can be found in the source directory at ``Nek5000/core/zero.usr``.
When setting up a new case, it is recommended to start by copying this file to your run directory when you don't already have a similar case to use as a template.

To aid in implementing custom physics, the *Nek5000* user file is partitioned into subroutines which provide access to specific, commonly used terms in the solver. 
These subroutines can be broadly divided into 3 categories:

- Local routines - called every time step for every field and GLL point
- Global routines - called once per time step
- Initialization routines - only called once during initialization

All of the routines in the ``.usr`` file have access to a common set of global variables and arrays. 
These are included in the ``TOTAL`` common block, included at the beginning of each subroutine.
Some of the more commonly used variables are shown in :numref:`tab:Globalvars`.
For a more complete list see :ref:`sec:commonvars`.

.. _tab:Globalvars:

.. csv-table:: Commonly used global variables available in ``TOTAL``
   :header: Variable,Description,Note
   :widths: 20,55,20

   ``pi``,:math:`\pi`,``pi=4.0atan(1.0)``
   ``time``,physical time,
   ``nelv``,number of elements in the velocity mesh on the local MPI rank,
   ``nelt``,number of elements in the temperature mesh on the local MPI rank, see :ref:`conjht`
   ``nelgv``,total number of elements in the velocity mesh, 
   ``nelgt``,total number of elements in the temperature mesh, see :ref:`conjht`
   ``idsess``,session ID for NekNek, see :ref:`neknek`
   ``ifield``,active solution field,see :numref:`tab:ifield`

.. Note:: 

  ``nelt`` is always greater than or equal to ``nelv``. 
  For non-conjugate heat transfer cases, ``nelv`` and ``nelt`` are identical. 
  Similarly for ``nelgv`` and ``nelgt``.

Also included in the ``TOTAL`` block is the field coefficient array.
This array contains reference values for fluid properties that are inherited directly from the ``.par`` file.
For a simulation with constant properties, this array corresponds directly to density, viscosity, etc. 
as the values from this array are copied directly into the ``vtrans`` and ``vdiff`` arrays internally by *Nek5000* (See :ref:`sec:uservp` for more information on ``vtrans`` and ``vdiff``).
For simulations with variable properties, this array can be useful as it retains the values assigned in the ``.par`` file.

.. _tab:cpfld:

.. csv-table:: The field coefficient array
   :header: Variable,Entry in ``.par`` file,Description

   "``cpfld(1,1)``",``velocity:viscosity``,Reference viscosity
   "``cpfld(1,2)``",``velocity:density``,Reference density
   "``cpfld(2,1)``",``temperature:conductivity``,Reference conductivity
   "``cpfld(2,2)``",``temperature:rhocp``,Reference rho-cp
   "``cpfld(3,1)``",``scalar01:diffusivity``,Reference diffusivity for scalar 1
   "``cpfld(3,2)``",``scalar01:density``,Reference density for scalar 1
   "``cpfld(i+2,1)``",``scalar i:diffusivity``,Reference diffusivity for :math:`i`
   "``cpfld(i+2,2)``",``scalar i:density``,Reference density for scalar :math:`i`

.. Note::
  The entries in the field coefficient array for velocity and temperature correspond directly to parameters 1, 2, 7, and 8 from the old ``.rea`` format.
  However, no corresponding parameters exist for the passive scalars.

.. _local_routines:

--------------
Local Routines
--------------

The local subroutines ``uservp``, ``userf``, and ``userq`` are called at every GLL point, for every field at every time step.
The subroutines ``userbc`` and ``useric`` are similar, but ``userbc`` is only called for GLL points on a boundary face with certain boundary conditions and ``useric`` is only called during initialization.
These subroutines take ``ix``, ``iy``, ``iz``, and ``eg`` as arguments, which correspond to the local GLL indexing and global element number.
Note that the local element number can be accessed via the ``gllel`` array.

.. _sec:NEKUSE:

......
NEKUSE
......

The ``ix``, ``iy``, etc. indices can be used to cross-reference the local solution arrays, however, the ``NEKUSE`` common block is provided for easier access.
This block contains solver variables that may be useful for defining custom models, such as variable properties, or a localized heating rate.
Immediately before a local routine is called by *Nek5000*, the subroutine ``nekasgn`` is called.
This routine sets many of the commonly used values in ``NEKUSE``, making them available for use in the local subroutines.
:numref:`tab:NEKUSEpre` describes the variables that are assigned in ``nekasgn``.

.. _tab:NEKUSEpre:

.. csv-table:: Prepopulated ``NEKUSE`` variables
   :header: Variable,Description,Solution Array,Note
   :widths: 15,50,20,15

   ``x``,x-coordinate,"``xm1(ix,iy,iz,ie)``",
   ``y``,y-coordinate,"``ym1(ix,iy,iz,ie)``",
   ``z``,z-coordinate,"``zm1(ix,iy,iz,ie)``",
   ``r``,r-coordinate,N/A,:math:`\sqrt{x^2+y^2}`
   ``theta``,:math:`\theta`-coordinate,N/A,"``theta=atan2(y,x)``"
   ``ux``,x-velocity,"``vx(ix,iy,iz,ie)``",
   ``uy``,y-velocity,"``vy(ix,iy,iz,ie)``",
   ``uz``,z-velocity,"``vz(ix,iy,iz,ie)``",
   ``temp``,temperature,"``t(ix,iy,iz,ie,1)``",
   ``ps(i)``,"passive scalar 'i', :math:`\phi_i`","``t(ix,iy,iz,ie,i+1)``",":math:`i=(` ``ifield`` :math:`-2)`"
   ``pa``,pressure,"``pr(ix,iy,iz,ie)``",:math:`P_N/P_N` only
   ``p0``,thermodynamic pressure,``p0th``,
   ``udiff``,diffusion coeffcient,"``vdiff(ix,iy,iz,ie,ifield)``","See :numref:`tab:uservp`"
   ``utrans``,convective coefficient,"``vtrans(ix,iy,iz,ie,ifield)``","See :numref:`tab:uservp`"

..   ``si2``,strain rate invariant II,"``sii(ix,iy,iz,ie)``",
     ``si3``,strain rate invarient III,"``siii(ix,iy,iz,ie)``",

The active solution field -- velocity, temperature, etc. -- is indicated by the global variable ``ifield``. 
The value that corresponds to each field is described in :numref:`tab:ifield`.
``ifield`` is not explicitly passed as an input to any of the user subroutines, rather it is controlled at a higher level directly by *Nek5000*.
It is set by the solver for all of the local routines, but not any of the initialization or global routines.

.. _tab:ifield:

.. csv-table:: Corresponding solution fields for ``ifield``
   :header: ``ifield`` value,Field name,Solution variable

   1,velocity,:math:`\mathbf u`
   2,temperature,:math:`T`
   3,passive scalar 1,:math:`\phi_1`
   4,passive scalar 2,:math:`\phi_2`
   :math:`\ge 5`,passive scalar :math:`(` ``ifield`` :math:`-2)`,:math:`\phi_{ifld-2}`

.. _sec:uservp:

...................
uservp
...................

This function is only called if the ``variableProperties`` key under ``[PROBLEMTYPE]`` in the ``.par`` file is set to true (see :ref:`here<case_files_par>`).
It can be used to specify customized or solution dependent material properties.
It is called for every GLL-point for every field at every time step.
The diffusion and transport coefficients should be set using the variables described in the table below.
The diffusion coefficient refers the viscosity, thermal conductivity, or diffusivity for passive scalars.
The transport coefficient refers to the coefficient attached to the convective term, typically density for velocity and passive scalars and the product of density and specific heat for temperature.

.. _tab:uservp:

.. table:: Terms set in ``uservp``

   +------------+-----------------------+-------------------------------------------------------------------------------+----------------------------+
   |            | Description           | Variable in governing equations                                               | solution field(s)          |
   +============+=======================+===============================================================================+============================+
   |            |                       | :math:`\mu` in the :ref:`momentum equation <intro_ns>`                        | ``ifield = 1``             |
   |            |                       +-------------------------------------------------------------------------------+----------------------------+ 
   | ``udiff``  | diffusion coefficient | :math:`\lambda` in the :ref:`energy equation <intro_energy>`                  | ``ifield = 2``             |
   |            |                       +-------------------------------------------------------------------------------+----------------------------+
   |            |                       | :math:`\Gamma_i` in the :ref:`passive scalar equations <intro_pass_scal>`     | ``ifield = 3 .. npscal+2`` |
   +------------+-----------------------+-------------------------------------------------------------------------------+----------------------------+
   |            |                       | :math:`\rho` in the :ref:`momentum equation <intro_ns>`                       | ``ifield = 1``             |
   |            |                       +-------------------------------------------------------------------------------+----------------------------+ 
   | ``utrans`` | transport coefficient | :math:`(\rho c_p)` in the :ref:`energy equation <intro_energy>`               | ``ifield = 2``             |
   |            |                       +-------------------------------------------------------------------------------+----------------------------+ 
   |            |                       | :math:`\rho_i` in the :ref:`passive scalar equations <intro_pass_scal>`       | ``ifield = 3 .. npscal+2`` |
   +------------+-----------------------+-------------------------------------------------------------------------------+----------------------------+

.. Warning::

  The corresponding entries in ``vdiff`` and ``vtrans`` are overwritten by whatever is assigned to ``udiff`` and ``utrans``. Setting ``vdiff`` and ``vtrans`` directly is not supported.

:Example:
  The code block below shows how to implement a variable viscosity as a function of temperature, with the density, rho-cp, and thermal conductivity set from the values in the ``.par`` file using the :ref:`field coefficient array <tab:cpfld>`.

.. literalinclude:: examples/uservp.txt
   :language: fortran

.. _userf:

...................
userf
...................

This functions sets the source term (which will be subsequently be multiplied by the density) for the momentum equation.
It allows the user to effectively add an acceleration term.


:Example:
  The code block below shows how to implement a body force proportional to the temperature, similar to what would be done for a Boussinesq model for buoyancy.

.. literalinclude:: examples/userf.txt
   :language: fortran

.. _sec:userq:

...................
userq
...................

This functions sets the source term for the :ref:`energy<intro_energy>` and :ref:`passive scalar<intro_pass_scal>` equations.
An explicit source term can be set using ``qvol``.
In the latest version available from the master branch on github, an implicit source term can be set using ``avol``.

A source term that has the form

.. math::
   q'''=\alpha -\beta T

can be implemented either entirely explicitly, or semi-implicitly.
In general, the implicit term should be used wherever possible as it tends to stabilize the solution.
Both approaches are shown below.
It is not necessary for :math:`\alpha` and :math:`\beta` to be constants.
They can vary with time, position, or any of the solution variables -- including temperature.
However, using a solution variable may impose limits on the stability of the solution.

:Example:
  In the first example, the source term is set entirely explicitly

.. literalinclude:: examples/userq1.txt
   :language: fortran

:Example:
  In the second example, the implicit source term is leveraged.
  Both implementations will yield the same solution for a converged result.

.. literalinclude:: examples/userq2.txt
   :language: fortran

.. _sec:userbc:

...................
userbc
...................

This functions sets boundary conditions. 
It is only called for special boundary condition types -- any lowercase value in the ``cbc`` array -- and only for points on the boundary surface.
It includes an additional argument compared to the other Local Routines.
The ``iside`` variable refers to which side of the element the boundary condition is on. 
This can be used for accessing the appropriate entry in the ``boundaryID`` or ``cbc`` arrays.
For a more complete list of all supported boundary conditions, see :ref:`boundary-conditions`.

The available values that can be set for velocity are listed in :numref:`tab:velbcs` along with their definitions and the corresponding entry in the ``cbc`` array, where :math:`\mathbf{\hat e}` denotes a unit vector.

.. _tab:velbcs:

.. csv-table:: Velocity boundary conditions set in ``userbc``
   :widths: 10,45,30,15
   :header:  ,Description,Definition,``cbc`` value

   ``ux``,x-velocity,":math:`\mathbf u\cdot\mathbf{\hat e_x}`",``v``
   ``uy``,y-velocity,":math:`\mathbf u\cdot\mathbf{\hat e_y}`",``v``
   ``uz``,z-velocity,":math:`\mathbf u\cdot\mathbf{\hat e_z}`",``v``
   ``un``,velocity normal to the boundary face,":math:`\mathbf u\cdot\mathbf {\hat e_n}`",``vl``   
   ``u1``,velocity tangent* to the boundary face,":math:`\mathbf u\cdot\mathbf {\hat e_t}`",``vl``   
   ``u2``,velocity bitangent* to boundary face,":math:`\mathbf u\cdot\mathbf {\hat e_b}`",``vl``   
   ``trx``,"traction in the x-direction",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_x}`","``s``, ``sh``"
   ``try``,"traction in the y-direction",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_y}`","``s``, ``sh``"
   ``trz``,"traction in the z-direction",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_z}`","``s``, ``sh``"
   ``trn``,"traction normal to the boundary face",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_n}`",``sl``
   ``tr1``,"traction tangent* to the boundary face",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_t}`","``sl``, ``shl``"
   ``tr2``,"traction bitangent* to the boundary face",":math:`(\boldsymbol{\underline\tau}\cdot\mathbf{\hat e_n})\cdot\mathbf{\hat e_b}`","``sl``, ``shl``"

.. Warning::

  \*The tangent and bitangent directions are not guaranteed to be consistent between elements in 3D domains.

The available values that can be set for temperature are listed in :numref:`tab:Tbcs` along with their definitions and the corresponding entry in the ``cbc`` array.
These correspond to standard Dirichlet, Neumann, and Robin boundary conditions.

.. _tab:tbcs:

.. table:: Temperature boundary conditions set in ``userbc``

   +----------+-----------------------------------------+---------------------------------------------------------------+---------------+
   |          | Description                             | Definition                                                    | ``cbc`` value |
   +==========+=========================================+===============================================================+===============+
   | ``temp`` | temperature                             | :math:`T`                                                     | ``t``         |
   +----------+-----------------------------------------+---------------------------------------------------------------+---------------+
   | ``flux`` | heat flux, :math:`q''`                  | :math:`\lambda\nabla T\cdot\mathbf{\hat e_n}=q''`             | ``f``         |
   +----------+-----------------------------------------+---------------------------------------------------------------+---------------+
   | ``hc``   | heat transfer coefficient, :math:`h`    | :math:`\lambda\nabla T\cdot\mathbf{\hat e_n}=h(T-T_{\infty})` | ``c``         |
   +----------+-----------------------------------------+                                                               |               |
   | ``tinf`` | ambient temperature, :math:`T_{\infty}` |                                                               |               |
   +----------+-----------------------------------------+---------------------------------------------------------------+---------------+

.. Note::
  Both heat transfer coefficient and ambient temperature must be specified for a Robin boundary condition.

A few examples are shown next to demonstrate how to set simple boundary conditions.

:Example: 
  In the example below, the code sets a parabolic inlet velocity with a constant inlet temperature of zero and a constant wall temperature of one. 
  The temperature field has the same BC of ``t``  on both the inlet and the wall, so the velocity BC is accessed to differentiate between the two. 
  Also note that this routine will not be called for ``ifield=1`` for the ``W`` boundary, but it will be called for ``ifield=2`` for the ``t`` boundary co-located with the ``W`` boundary.

.. literalinclude:: examples/userbc1.txt
   :language: fortran

Boundary conditions are only applied to boundary faces with the corresponding ``cbc`` array value.
In the example above, the ``cbc`` array does not need to be checked to assign the parabolic velocity profile as it will be ignored for any BC that is not ``v``.

.. _userbc_ex2:

:Example:
  In this example, the ``boundaryID`` array is used to differentiate between the inlet and two different walls.
  The inlet (ID = 1) has a velocity profile and constant temperature value.
  The walls (IDs 2 and 3 respectively) are set as a positive heat flux on wall 2 and a negative (cooling) heat flux on wall 3.
  This example corresponds to the example setup shown in :ref:`usrdat <usrdat_ex>`.

.. literalinclude:: examples/userbc2.txt
   :language: fortran

The temperature and passive scalars use the same form of governing equation and are handled practically identically by *Nek5000*.
The boundary conditions available for passive scalars are therefore the same as those available for temperature.
These are listed in :numref:`tab:psbcs` along with their definitions and the corresponding entry in the ``cbc`` array.

.. _tab:psbcs:

.. table:: Passive scalar boundary conditions set in ``userbc``

   +----------+-----------------------------------------+----------------------------------------------------------------------------------+---------------+
   |          | Description                             | Definition                                                                       | ``cbc`` value |
   +==========+=========================================+==================================================================================+===============+
   | ``temp`` | Dirichlet                               | :math:`\phi_i`                                                                   | ``t``         |
   +----------+-----------------------------------------+----------------------------------------------------------------------------------+---------------+
   | ``flux`` | Neumann, :math:`\psi_i`                 | :math:`\Gamma_i\nabla \phi_i\cdot\mathbf{\hat e_n}=\psi_i`                       | ``f``         |
   +----------+-----------------------------------------+----------------------------------------------------------------------------------+---------------+
   | ``hc``   | transfer coefficient, :math:`\eta`      | :math:`\Gamma_i\nabla \phi_i\cdot\mathbf{\hat e_n}=\eta(\phi_i-\phi_{i,\infty})` | ``c``         |
   +----------+-----------------------------------------+                                                                                  |               |
   | ``tinf`` | ambient value, :math:`\phi_{i,\infty}`  |                                                                                  |               |
   +----------+-----------------------------------------+----------------------------------------------------------------------------------+---------------+

.. Note::
  Temperature and passive scalar boundary conditions are contextual and will set the boundary condition for the active solution field. See :numref:`tab:ifield`.

.. Warning::
  The ``ps(i)`` variable array provided by ``NEKUSE`` is **NOT** used to set the passive scalar boundary conditions.

:Example:
  In this example, the ``ifield`` variable is used to distinguish between a constant inlet temperature of 0 and a transient concentration for passive scalar 1 which is distributed normally in time.

.. literalinclude:: examples/userbc3.txt
   :language: fortran

It is not uncommon to use multiple methods of discerning between different boundary conditions.
The user may need to implement specific conditions for different fields on different boundaries at different times.
This can quickly lead to complex, embedded if-then statements.
Use care to ensure you're setting the correct boundary conditions!

.. _sec:useric:

...................
useric
...................

This functions sets the initial conditions and behaves similarly to ``userbc``.
It is called only during initialization after ``usrdat3`` for every solution field on every GLL point in the domain.

.. Warning::
  ``useric`` is **NOT** called for any fields loaded from a restart file.

.. .. literalinclude:: examples/useric.txt
   :language: fortran

.. _global_routines:

---------------
Global Routines
---------------

.. _sec:userchk:

...................
userchk
...................

This is a general purpose routine that is executed both during initialization and after every time step.
It can be used for a wide variety of purposes, such as monitoring the solution, post-processing results, or implementing entirely new physical solvers.

It is helpful for users to familiarize themselves with the included utility subroutines in *Nek5000* as they can make performing complex calculations on the entire dataset much easier.
For a list of some of the commonly used subroutines, see :ref:`here <append_subroutines>`.

:Example: 
  Most of the code shown below will only be executed if *Nek5000* is run in "post-processing" mode, i.e. ``numSteps = 0`` is specified in the ``.par`` file. 
  The solution is rescaled by reference length, velocity, and pressure and written to an output file. 
  The name of the output file will be prepended with the characters ``dim`` and it will contain dimensional quantities that can be loaded into Paraview or VisIt for further post processing.

.. literalinclude:: examples/userchk.txt
   :language: fortran

.. _sec:userqtl:

...................
userqtl
...................

This function is only called if the ``equation`` key under ``[PROBLEMTYPE]`` in the ``.par`` file is set to ``lowMachNS`` (see :ref:`here<case_files_par>`).
This function can be used to specify a customized thermal divergence for the low Mach solver.
The thermal divergence refers to a non-zero right hand side of the divergence constraint (see :ref:`intro_low_mach`).

*Nek5000* includes a thermal divergence model for a single component ideal gas as denoted by the call to ``userqtl_scig``.
To use this model, simply leave ``userqtl`` as it is in the template file, ``zero.usr``.
The implementation of any other model is left to the user.

.. _initialization_routines:

-----------------------
Initialization Routines
-----------------------

.. _sec:usrdat:

...................
usrdat
...................

This function can be used to modify the element vertices and is called before the spectral element mesh (GLL points) has been laid out.
It can be used to fill the ``cbc`` array based on ``BoundaryID`` for 3rd party meshes.

.. _usrdat_ex:

:Example: 
  In the code below, the ``cbc`` array is filled for a 3rd party mesh. 
  The ``boundaryID`` array is filled with the Boundary ID values set in Gmsh or the sideset numbers specified in an exodus mesh file.
  This example corresponds with the setup shown in the :ref:`second example for userbc <userbc_ex2>`.

.. literalinclude:: examples/usrdat.txt
   :language: fortran

.. _sec:usrdat2:

...................
usrdat2
...................

This function can be used to modify the spectral element mesh.
The geometry information (mass matrix, surface normals, etc.) will be rebuilt after this routine is called.
Any changes to the ``cbc`` array must be made before or during this call.

:Example:
  In the code below, the mesh is scaled by a factor of :math:`1/D_h`, where :math:`D_h` is the hydraulic diameter.
  It is typically convenient to generate a mesh with dimensions, this allows the user to non-dimensionalize the mesh at runtime.

.. literalinclude:: examples/usrdat2.txt
   :language: fortran

.. _sec:usrdat3:

...................
usrdat3
...................

This function can be used to initialize any additional case/user specific data.
The GLL mesh should not be further modified in this routine as the geometric factors are not automatically recomputed.

