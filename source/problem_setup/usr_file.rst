.. _case_files_usr:

=========================
User Routines File (.usr)
=========================

The user file is a special case file that implements the the user interface to *Nek5000*. 
What follows is a brief description of the available subroutines and a primer on writing your own user code.
An empty user file template can be found in the source directory ``Nek5000/core/zero.usr``.
When setting up a brand new case, it is recommended to start by copying this file to your run directory.

To aid in implementing custom physics, the *Nek5000* user file is partitioned into subroutines which provide access to specific, commonly used terms in the solver. 
These subroutines can be broadly divided into 3 categories 

- Local routines - called every time step for every field and GLL point
- Global routines - called once per time step
- Intialization routines - only called once during intialization

.. _local_routines:

--------------
Local Routines
--------------

In many of the subroutines available in the ``.usr`` file include the ``NEKUSE`` common block. 
This block contains solver variables that may be useful for defining custom models, such as variable properties, or a localized heating rate.

The following variables are assigned in the subroutine ``nekasgn``, which is called before the subroutines ``uservp``, ``userf``, ``userq``, ``userbc``, and ``useric``.
These subroutines take ``ix``, ``iy``, ``iz``, and ``eg`` as arguments, which correspond to the local GLL indexing and global element number.
The global element number is translated into a local element number and ``nekasgn`` fills the variable with the corresponding entry from the solution array.

.. _tab:Globalvars:

.. csv-table:: Commonly used global variables available in ``TOTAL``
   :header: Variable,Description,Note
   :widths: 20,55,20

   ``pi``,:math:`\pi`,``pi=4.0atan(1.0)``
   ``time``,physical time,
   ``dt``,time step size,

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
   ``ps(i)``,passive scalar \"i\","``t(ix,iy,iz,ie,i+1)``",
   ``pa``,pressure,"``pr(ix,iy,iz,ie)``",not recommended for use with :math:`P_N/P_{N-2}`
   ``p0``,thermodynamic pressure,``p0th``,
   ``udiff``,diffusion coeffcient,"``vdiff(ix,iy,iz,ie,ifield)``","viscosity, conductivity, or diffusivity"
   ``utrans``,convective coefficient,"``vtrans(ix,iy,iz,ie,ifield)``","density, rho-cp, etc."

..   ``si2``,strain rate invariant II,"``sii(ix,iy,iz,ie)``",
     ``si3``,strain rate invarient III,"``siii(ix,iy,iz,ie)``",

.. _tab:NEKUSEvar:

.. table:: ``NEKUSE`` common block variables

   +-----------------------------+-----------------------------------------------------------------+
   |   Variable                  | | Description                                                   |
   +=============================+=================================================================+
   | ``ux`` , ``uy`` , ``uz``    | | ``x`` , ``y`` , ``z`` velocity components                     |
   +-----------------------------+-----------------------------------------------------------------+
   | ``un`` , ``u1`` , ``u2``    | | ``x`` , ``y`` , ``z`` velocity component of face unit normal  |
   +-----------------------------+-----------------------------------------------------------------+
   | ``trx`` , ``try`` , ``trz`` | | ``x`` , ``y`` , ``z`` traction components                     |
   +-----------------------------+-----------------------------------------------------------------+
   | ``trn`` , ``tr1`` , ``tr2`` | | ``x`` , ``y`` , ``z`` traction component of face unit normal  |
   +-----------------------------+-----------------------------------------------------------------+
   | ``pa`` , ``p0``             | | Outlet pressure, system pressure                              |
   +-----------------------------+-----------------------------------------------------------------+
   | ``ffx`` , ``ffy`` , ``ffz`` | | ``x`` , ``y`` , ``z`` acceleration                            |
   +-----------------------------+-----------------------------------------------------------------+
   | ``temp``                    | | Temperature                                                   |
   +-----------------------------+-----------------------------------------------------------------+
   | ``flux``                    | | Heat flux                                                     |
   +-----------------------------+-----------------------------------------------------------------+
   | ``hc`` , ``hrad``           | | Heat transfer coefficient (convective, radiative)             |
   +-----------------------------+-----------------------------------------------------------------+
   | ``tinf``                    | | Temperature at infinity                                       |
   +-----------------------------+-----------------------------------------------------------------+
   | ``qvol`` , ``avol``         | | Source terms for temperature and passive scalars              |
   +-----------------------------+-----------------------------------------------------------------+
   | ``sigma``                   | | Surface-tension coefficient                                   |
   +-----------------------------+-----------------------------------------------------------------+
   | ``ps``                      | | Passive scalars                                               |
   +-----------------------------+-----------------------------------------------------------------+

.. _case_files_uservp:

...................
uservp
...................

This function can be used to specify customized or solution dependent material properties.
It is called for every GLL-point for every field at every time step.
The diffusion and transport coefficients should be set using the variables described in the table below.
The diffusion coefficient refers to the variable :math:`\mu` in the :ref:`intro_ns`, the variable :math:`\lambda` in the :ref:`intro_energy`, and the variable :math:`\Gamma_i` in the :ref:`intro_pass_scal` transport equation.
The transport coefficient refers to the coefficient attached to the convective term, i.e., variable :math:`\rho` in the :ref:`intro_ns`, the variables :math:`(\rho c_p)` in the :ref:`intro_energy`, and the variables :math:`(\rho c_p)_i` in the :ref:`intro_pass_scal` transport equation.

.. csv-table::
   :header: Variable,Description,Note
   :widths: 20,55,20

   ``udiff``,diffusion coefficient,"viscosity, conductivity, or diffusivity"
   ``utrans``,transport coefficient,"density, rho-cp, etc. "

:Example:

.. code-block:: fortran

      if (ifield.eq.1) then
         udiff  = a * exp(-b*temp) ! dynamic viscosity
         utrans = 1.0              ! density
      else if (ifield.eq.2) then
         udiff  = 1.0              ! conductivity
         utrans = 1.0              ! rho*cp
      endif

...................
userf
...................

This functions sets the source term (which will be subsequently be multiplied by he density) for the momentum equation.
It allows the user to effectively add an acceleration term.


:Example:
  The code block below shows how to implement gravity in the z-direction

.. code-block:: fortran

      real g
      parameter(g = 9.81)

      ffx = 0.0
      ffy = 0.0
      ffz = -g ! gravitational acceleration

...................
userq
...................

This functions sets the source term for the energy (temperature) and passive scalar equations.
An explicit source term can be set using ``qvol``.
In the latest version availble from the master branch on github, an implicit source term can be set using ``avol``.

...................
userbc
...................

This functions sets boundary conditions. 
Note, this function is only called for special boundary condition types and only for points on the boundary surface.
It includes an additional argument compared to the other Local Routines.
The ``iside`` variables refers to which side of the element the boundary condition is on. 
This can be used for accessing the appropriate entery in the ``boundaryID`` or ``cbc`` arrays.

:Example: 
  In the example below, the code sets a parabolic inlet velocity with a constant inlet temperature of 0.0 and a constant wall temperature of 1.0. 
  The temperature field has the same BC of ``t``  on both the inlet and the wall, so the velocity BC is accessed to differentiate between the two. 
  Also note that this routine will not be called for ``ifield=1`` for the ``W`` boundary, but it will be called for ``ifield=2`` for the ``t`` boundary colocated with the ``W`` boundary.

.. code-block:: fortran

  integer ie
  character*3 cb3

  ie=gllel(eg) !get local element number 
  cb3=cbc(iside,ie,1) !access the velocity boundary condition

  uz = 3./2. (1.0-(2.0*y-1.0)**2

  if(cb3.eq.'v  ')
    temp = 0.0 !set inlet temperature to 0.0
  elseif(cb3.eq.'W  ')
    temp = 1.0 !set wall temperature to 1.0
  endif

:Example:
  In this example, the ``boundaryID`` array is used to set a positive heat flux on wall 1 and a negative (cooling) heat flux on wall 2.

.. code-block:: fortran

  integer ie
  
  ie=gllel(eg)  !get local element number

  if(boundaryID(iside,ie).eq.1)
    flux = 1.0
  elseif(boundaryID(iside,ie).eq.2)
    flux = -1.0
  endif

...................
useric
...................

This functions sets the initial conditions.

.. _global_routines:

---------------
Global Routines
---------------

...................
userchk
...................

This is a general purpose routine that gets executed both during intialization and after every time
step.

...................
userqtl
...................

This function can be used  to specify a cutomzized thermal diveregence for the low Mach solver.
step.

.. _initialization_routines:

-----------------------
Initialization Routines
-----------------------

...................
usrdat
...................

This function can be used to modify the element vertices and is called before the spectral element mesh (GLL points) has been laid out.

...................
usrdat2
...................

This function can be used to modify the spectral element mesh.
The geometry information (mass matrix, surface normals, etc.) will be rebuilt after this routine is called.

...................
usrdat3
...................

This function can be used to initialize case/user specific data.

