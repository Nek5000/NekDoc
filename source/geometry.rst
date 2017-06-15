========
Geometry
========

-----------------------------
Generating a Mesh with Genbox
-----------------------------

Note that in case of any changes in the SIZE file, a recompilation is necessary.

..........................
Uniformly Distributed Mesh
..........................

Suppose you wish to simulate flow through an axisymmetric pipe,
of radius :math:`R=0.5` and length :math:`L=4`.  You estimate that you will
need 3 elements in radial :math:`(y)` direction, and 5 in the :math:`x` direction,
as depicted in :numref:`fig:mesh_axi1`.
This would be specified by the following input file (called *pipe.box*)
to genbox:

.. code-block:: none

   axisymmetric.rea
   2                      spatial dimension
   1                      number of fields
   #
   #    comments:   This is the box immediately behind the
   #                refined cylinder in Ugo's cyl+b.l. run.
   #
   #
   #========================================================
   #
   Box 1                         Pipe
   -5 -3                         Nelx  Nely
   0.0   4.0   1.0               x0  x1   ratio
   0.0   0.5   1.0               y0  y1   ratio
   v  ,O  ,A  ,W  ,   ,          BC's:  (cbx0, cbx1, cby0, cby1, cbz0, cbz1)

.. _fig:mesh_axi1:

.. figure:: figs/mesh_axi1.png
    :align: center
    :figclass: align-center
    :alt: axis-pipe-mesh

    Axisymmetric pipe mesh.

- The first line of this file supplies the name of an existing 2D .rea file that has the appropriate run parameters (viscosity, timestep size, etc.). These parameters can be modified later, but it is important that axisymmetric.rea be a 2D file, and not a 3D file.
- The second line indicates the number of fields for this simulation, in this case, just 1, corresponding to the velocity field (i.e., no heat transfer).
- The next set of lines just shows how one can place comments into a genbox input file.
- The line that starts with "Box" indicates that a new box is starting, and that the following lines describe a typical box input.  Other possible key characters (the first character of Box, "B") are "C" and "M", more on those later.
- The first line after "Box" specifies the number of elements in the
  :math:`x` and :math:`y` directions.   The fact that these values are negative indicates
  that you want genbox to automatically generate the element distribution
  along each axis, rather than providing it by hand.  (More on this below.)
- The next line specifies the distribution of the 5 elements in the :math:`x` direction.
  The mesh starts at :math:`x=0` and ends at :math:`x=4.0`.  The *ratio* indicates the
  relative size of each element, progressing from left to right.
- The next line specifies the distribution of the 3 elements in the :math:`y` direction,
  starting at :math:`y=0` and going to :math:`y=0.5`.  Again,
  *ratio* =1.0 indicates that the elements will be of uniform height.
- The last line specifies boundary conditions on each of the 4 sides of the
  box:

  - Lower-case *v* indicates that the left :math:`(x)` boundary is to be a velocity
    boundary condition, with a user-specified distribution determined by
    routine *userbc* in the .usr file.  (Upper-case :math:`V` would indicate that
    the velocity is constant, with values specified in the .rea file.)
  - *O* indicates that the right :math:`(x)` boundary is an outflow boundary -- the
    flow leaves the domain at the left and the default exit pressure is :math:`p=0`.
  - *A* indicates that the lower :math:`(y)` boundary is the axis---this condition
    is mandatory for the axisymmetric case, given the fact that the lower domain
    boundary is at :math:`y=0`, which corresponds to :math:`r=0`.
  - *W* indicates that the upper :math:`(y)` boundary is a wall.  This would be
    equivalent to a *v* or *V* boundary condition, with :math:`{\bf u}=0`.

...........
Graded Mesh
...........

.. _fig:mesh_axi2:

.. figure:: figs/mesh_axi2.png
    :align: center
    :figclass: align-center
    :alt: axis-pipe-mesh-graded

    Axisymmetric pipe mesh, graded

Suppose you wish to have the mesh be graded,
that you have increased resolution near the wall.
In this case you change *ratio* in the :math:`y`-specification
of the element distribution.  For example, changing the 3 lines
in the above genbox input file from

.. code-block:: none

   -5 -3                         Nelx  Nely
   0.0   4.0   1.0               x0  x1   ratio
   0.0   0.5   1.0               y0  y1   ratio

to

.. code-block:: none

   -5 -4                         Nelx  Nely
   0.0   4.0   1.0               x0  x1   ratio
   0.0   0.5   0.7               y0  y1   ratio

yields the mesh shown in :numref:`fig:mesh_axi2`.

...........................
User-Specified Distribution
...........................

.. _fig:mesh_axi3:

.. figure:: figs/mesh_axi3.png
    :align: center
    :figclass: align-center
    :alt: axis-pipe-mesh-user

    Axisymmetric pipe mesh, user specified.

You can also specify your own, precise, distribution of element
locations.   For example, another graded mesh similar to the
one of the preceding example could be built by changing the
genbox input file to contain:

.. code-block:: none

   -5  4                                               Nelx  Nely
   0.0   4.0   1.0                                     x0  x1   ratio
   0.000    0.250    0.375    0.450    0.500           y0  y1 ... y4

Here, the positive number of elements for the :math:`y` direction indicates
that genbox is expecting ``Nely+1`` values of :math:`y` positions on the
:math:`y`-element distribution line.   This is the genbox default, which
explains why it corresponds to ``Nely`` :math:`>` 0.  The corresponding mesh
is shown in :numref:`fig:mesh_axi3`.

............................
Mesh Modification in Nek5000
............................

For complex shapes, it is often convenient to modify the mesh
direction in the simulation code, Nek5000.  This can be done
through the usrdat2 routine provided in the .usr file.
The routine usrdat2 is called by nek5000 immediately after
the geometry, as specified by the .rea file, is established.
Thus, one can use the existing geometry to map to a new geometry
of interest.

For example, suppose you want the above pipe geometry to have
a sinusoidal wall.  Let :math:`{\bf x} := (x,y)` denote the old geometry,
and :math:`{\bf x}' := (x',y')` denote the new geometry.  For a domain
with :math:`y\in [0,0.5]`, the following function will map the straight
pipe geometry to a wavy wall with amplitude :math:`A`, wavelength :math:`\lambda`:

.. math::

    y'(x,y) = y  + y A \sin( 2 \pi x / \lambda ).

Note that, as :math:`y \longrightarrow 0`, the perturbation,
:math:`yA \sin( 2 \pi x / \lambda )`, goes to zero.  So, near the axis,
the mesh recovers its original form.

In nek5000, you would specify this through usrdat2 as follows

.. code-block:: fortran

   subroutine usrdat2
   include 'SIZE'
   include 'TOTAL'

   real lambda

   ntot = nx1*ny1*nz1*nelt

   lambda = 3.
   A      = 0.1

   do i=1,ntot
      argx         = 2*pi*xm1(i,1,1,1)/lambda
      ym1(i,1,1,1) = ym1(i,1,1,1) + ym1(i,1,1,1)*A*sin(argx)
   end do

   param(59) = 1.  ! Force nek5 to recognize element deformation.

   return
   end

Note that, since nek5000 is modifying the mesh, postx will not
recognize the current mesh unless you tell it to, because postx
looks to the .rea file for the mesh geometry.  The only way for
nek5000 to communicate the new mesh to postx is via the .fld
file, so you must request that the geometry be dumped to the
.fld file.   This is done by modifying the OUTPUT SPECIFICATIONS,
which are found near the bottom of the .rea file.  Specifically,
change

.. code-block:: none

   ***** OUTPUT FIELD SPECIFICATION *****
   6 SPECIFICATIONS FOLLOW
   F      COORDINATES
   T      VELOCITY
   T      PRESSURE
   T      TEMPERATURE
   F      TEMPERATURE GRADIENT
   0      PASSIVE SCALARS

to

.. code-block:: none

   ***** OUTPUT FIELD SPECIFICATION *****
   6 SPECIFICATIONS FOLLOW
   T      COORDINATES                       <------  CHANGE HERE
   T      VELOCITY
   T      PRESSURE
   T      TEMPERATURE
   F      TEMPERATURE GRADIENT
   0      PASSIVE SCALARS

The result of above changes is shown in :numref:`fig:wavypipe`.

.. _fig:wavypipe:

.. figure:: figs/wavypipe.png
    :align: center
    :figclass: align-center
    :alt: axis-pipe-mesh-wavy

    Axisymmetric pipe mesh.

.......................................
Cylindrical/Cartesian-transition Annuli
.......................................

.. _fig:cylbox_2d:

.. figure:: figs/cylbox_2d.png
    :align: center
    :figclass: align-center
    :alt: annuli-mesh-1

    Cylinder mesh

.. _fig:cylbox_2da:

.. figure:: figs/cylbox_2da.png
    :align: center
    :figclass: align-center
    :alt: annuli-mesh-2

    Cylinder mesh

An updated version of genb6, known as genb7, is currently under development
and designed to simply/automate the construction of cylindrical annuli,
including *basic* transition-to-Cartesian elements.   More sophisticated
transition treatments may be generated using the GLOBAL REFINE options in
prenek or through an upgrade of genb7, as demand warrants.
Example 2D and 3D input files are provided in the nek5000/doc files
*box7.2d* and *box7.3d*.
:numref:`fig:cylbox_2d` shows a 2D example generated using
the *box7.2d* input file, which reads:

.. code-block:: none

   x2d.rea
   2                      spatial dimension
   1                      number of fields
   #
   #    comments
   #
   #
   #========================================================
   #
   Y                   cYlinder
   3 -24 1             nelr,nel_theta,nelz
   .5 .3               x0,y0 - center of cylinder
   ccbb                descriptors: c-cyl, o-oct, b-box (1 character + space)
   .5 .55 .7 .8        r0 r1 ... r_nelr
   0  1  1             theta0/2pi theta1/2pi  ratio
   v  ,W  ,E  ,E  ,    bc's (3 characters + comma)
    
An example of a mesh is shown in :numref:`fig:cylbox_2d`.   The mesh has been quad-refined
once with oct-refine option of prenek. The 3D counterpart to this
mesh could joined to a hemisphere/Cartesian transition built with
the spherical mesh option in prenek.

-----------------------
Extrusion and Mirroring
-----------------------

---------------
Moving Geometry
---------------

-------------------------------
Boundary and Initial Conditions
-------------------------------

--------------------------------------
Parallel Mesh Partitioning with Genmap
--------------------------------------
