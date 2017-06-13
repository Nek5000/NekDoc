.. _user_files:

==========
User Files
==========

Each simulation is defined by three files: the .rea file, the .usr file, and the SIZE file.  In
addition, there is a derived .map file that is generated from the .rea file by running *genmap*,
which will determine how the elements will be split across processors in the case of a parallel
run.  SIZE controls (at compile time) the polynomial degree used in the simulation, as well as the
space dimension :math:`d=2` or :math:`3`.

The SESSION.NAME file provides the name and path of the .rea file and the path to it.  It does not
however need to correspond to a .usr file of an identical name. This allows for different test
cases (.usr files) that use the same geometry and boundary conditions (.rea files).

This chapter provides an introduction to the basic files required to set up a Nek5000 simulation.

.. _user_files_session:

------------
SESSION File
------------

To run NEK5000, each simulation must have a SESSION.NAME file. This file is read in by the code and
gives the path to the relevant files describing the structure and parameters of the simulation. The
SESSION.NAME file is a file that contains the name of the simulation and the full path to
supporting files. For example, to run the eddy example from the repository, the SESSION.NAME file
would look like:

  eddy_uv
  /home/user_name/Nek5000/short_tests/eddy/ 

.. _user_files_usr:

----------------------
Case Setup File (.usr)
----------------------

.....................
Contents of .usr File
.....................


The most important interface to Nek5000 is the set of Fortran subroutines that are contained in the
*.usr* file.  This file allows direct access to all runtime variables.  Here, the user may
specify spatially varying properties (e.g., viscosity), volumetric heating sources, body forces,
and so forth.  One can also specify arbitrary initial and boundary conditions through the routines
``useric`` and ``userbc``.  The routine ``userchk`` allows the user to interrogate the
solution at the end of each timestep for diagnostic purposes.   The *.usr* files provided in
the *Nek5000/short_tests/* directory illustrate several of the more common analysis tools.  For
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
  enddo
  enddo

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

The boundary conditions are assigned in ``userbc`` for both the fluid, temperature and all other scalars. An extensive list of such possible boundary conditions is available in Section.~\ref{sec:boundary}. 

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

alternatively the variable properties can be set in the USERVP routine.

**What is a SESSION file?**

To run NEK5000, each simulation must have a SESSION.NAME file. This file is read in by the code and gives the path to the relevant files describing the structure and parameters of the simulation. The SESSION.NAME file is a file that contains the name of the simulation and the full path to supporting files. For example, to run the eddy example from the repository, the SESSION.NAME file would look like::

  eddy_uv\\
  /homes/user\_ name/nek5\_ svn/examples/eddy/


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
* **lelt** determines the *maximum* number of elements *per processor}*

The total size of the problem is ``lx1*ly1*lz1*lelt``.

...................
Memory Requirements
...................

.. highlight:: bash

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
    use as handles to be passed into the user defined routines in the .usr file.
**passive scalar data** 
    This information can be specified also in the ``.uservp`` routine in the .usr 
    file. If specified in the .rea file then the coefficients for the conductivity 
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


Next we have the logical switches as follow, a detailed explanation to be found in Sec:\ref{sec:switches} 

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

**geometry**
    The geometry is specified in an arcane format specifying
    the :math:`xyz` locations of each of the eight points for each element,
    or the :math:`xy` locations of each of the four points for each element in 2D.
    A line of the following type may be encountered at the beginning 
    of the mesh section of the area file::

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

    Following the header, all elements are listed. The fluid elements are listed 
    first, followed by all solid elements if present. In this case there are (34) 
    solid elements.

    The data following the header is formatted as shown in Table :numref:`tab:element`. This provides all the coordinates of an element for top and bottom faces. The numbering of the vertices is shown in Fig. :numref:`fig:elorder`. The header for each element as in Table. :numref:`tab:element`, i.e. ``[1A] GROUP`` is reminiscent of older Nek5000 format and does not impact the mesh generation at this stage. (We are inquiring whether other groups still use it.)

    .. _tab:element:

    .. table:: Geometry description in .rea file

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

(NECESSARY TEXT FOR SOME REASON?)

    .. _fig:elorder:

    .. figure:: figs/3dcube_1.png
        :align: center
        :figclass: align-center
        :alt: rea-geometry

        Geometry description in .rea file (sketch of one element ordering - Preprocessor 
        corner notation) 

-----------
Data Layout
-----------

