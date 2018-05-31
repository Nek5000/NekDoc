.. _mesh_gen: 

-----------------------------
Generating a Mesh with Genbox
-----------------------------

..........................
Uniformly Distributed Mesh
..........................

Suppose you wish to simulate flow through an axisymmetric pipe,
of radius :math:`R=0.5` and length :math:`L=4`.  You estimate that you will
need 3 elements in radial :math:`(y)` direction, and 5 in the :math:`x` direction,
as depicted in :numref:`fig:mesh_axi1`.
This would be specified by the following input file (called ``pipe.box``)
to ``genbox``:

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

- The first line of this file supplies the name of an existing 2D ``.rea`` file that has the appropriate run parameters (viscosity, timestep size, etc.). These parameters can be modified later, but it is important that ``axisymmetric.rea`` be a 2D file, and not a 3D file.
- The second line indicates the number of fields for this simulation, in this case, just 1, corresponding to the velocity field (i.e., no heat transfer).
- The next set of lines just shows how one can place comments into a ``genbox`` input file.
- The line that starts with "Box" indicates that a new box is starting, and that the following lines describe a typical box input.  Other possible key characters (the first character of Box, "B") are "C" and "M", more on those later.
- The first line after "Box" specifies the number of elements in the
  :math:`x` and :math:`y` directions.   The fact that these values are negative indicates
  that you want ``genbox`` to automatically generate the element distribution
  along each axis, rather than providing it by hand.  (More on this below.)
- The next line specifies the distribution of the 5 elements in the :math:`x` direction.
  The mesh starts at :math:`x=0` and ends at :math:`x=4.0`.  The ``ratio`` indicates the
  relative size of each element, progressing from left to right.
- The next line specifies the distribution of the 3 elements in the :math:`y` direction,
  starting at :math:`y=0` and going to :math:`y=0.5`.  Again,
  ``ratio`` =1.0 indicates that the elements will be of uniform height.
- The last line specifies boundary conditions on each of the 4 sides of the
  box:

  - Lower-case *v* indicates that the left :math:`(x)` boundary is to be a velocity
    boundary condition, with a user-specified distribution determined by
    routine ``userbc`` in the ``.usr`` file.  (Upper-case :math:`V` would indicate that
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
In this case you change ``ratio`` in the :math:`y`-specification
of the element distribution.  For example, changing the 3 lines
in the above ``genbox`` input file from

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
``genbox`` input file to contain:

.. code-block:: none

   -5  4                                               Nelx  Nely
   0.0   4.0   1.0                                     x0  x1   ratio
   0.000    0.250    0.375    0.450    0.500           y0  y1 ... y4

Here, the positive number of elements for the :math:`y` direction indicates
that ``genbox`` is expecting ``Nely+1`` values of :math:`y` positions on the
:math:`y`-element distribution line.   This is the ``genbox`` default, which
explains why it corresponds to ``Nely`` :math:`>` 0.  The corresponding mesh
is shown in :numref:`fig:mesh_axi3`.

-----------------------
Mesh Modification
-----------------------

For complex shapes, it is often convenient to modify the mesh
direction in the simulation code, Nek5000.  This can be done
through the ``usrdat2`` routine provided in the ``.usr`` file.
The routine ``usrdat2`` is called by Nek5000 immediately after
the geometry, as specified by the ``.rea`` file, is established.
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

In Nek5000, you would specify this through ``usrdat2`` as follows

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

Note that, since Nek5000 is modifying the mesh, ``postx`` will not
recognize the current mesh unless you tell it to, because ``postx``
looks to the ``.rea`` file for the mesh geometry.  The only way for
Nek5000 to communicate the new mesh to ``postx`` is via the ``.fld``
file, so you must request that the geometry be dumped to the
``.fld`` file.  
The result of above changes is shown in :numref:`fig:wavypipe`.

.. _fig:wavypipe:

.. figure:: figs/wavypipe.png
    :align: center
    :figclass: align-center
    :alt: axis-pipe-mesh-wavy

    Axisymmetric pipe mesh.

.. _sec:genbox:

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

More sophisticated
transition treatments may be generated using the GLOBAL REFINE options in
``prenek`` or through an upgrade of ``genb7``, as demand warrants.
Example 2D and 3D input files are provided in the ``nek5000/doc`` files
``box7.2d`` and ``box7.3d``.
:numref:`fig:cylbox_2d` shows a 2D example generated using
the ``box7.2d`` input file, which reads:

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
once with oct-refine option of ``prenek``. The 3D counterpart to this
mesh could joined to a hemisphere/Cartesian transition built with
the spherical mesh option in ``prenek``.

-----------------------
Mesh Extrusion and Mirroring
-----------------------

In ``nek5000/tools``, there is a code ``n2to3.f`` that can be compiled with your
local fortran compiler (preferably not g77).
By running this code, you can extend two dimensional domains to
three dimensional ones with a user-specified number of levels in the
:math:`z`-direction.  Such a mesh can then be modified using the mesh modification
approach. Assuming you have a valid two-dimensional mesh, ``n2to3`` is straightforward
to run.  Below is a typical session, upon typing ``n2to3`` the user is prompted at the command line

.. code-block:: none

    Input old (source) file name:
   h2e
    Input new (output) file name:
   h3e
    input number of levels: (1, 2, 3,... etc.?):
   16
    input z min:
   0
    input z max:
   16
    input gain (0=custom,1=uniform,other=geometric spacing):
   1
    This is for CEM: yes or no:
   n
    Enter Z (5) boundary condition (P,v,O):
   v
    Enter Z (6) boundary condition (v,O):
   0
    this is cbz: v  O   <---

         320 elements written to h3e.rea
   FORTRAN STOP

In this context CEM stands for computational electromagnetics, a spin-off of the original Nek5000 code.

The domain in which the fluid flow/heat transfer
problem is solved consists of two distinct subdomains. The
first subdomain is that part of the region occupied by
fluid, denoted :math:`\Omega_f`, while the second subdomain is that part
of the region occupied by a solid, denoted :math:`\Omega_s`. These two
subdomains are depicted in :numref:`fig-walls`. The entire domain is denoted as :math:`D=\Omega_f \cup \Omega_s`.
The fluid problem is solved in the domain :math:`\Omega_f`, while the
temperature in the energy equation is solved in the
entire domain; the passive scalars can be solved in either
the fluid or the entire domain.

We denote the entire boundary of :math:`\Omega_f` as :math:`\partial \Omega_f`, that part
of the boundary of :math:`\Omega_f` which is not shared by :math:`\Omega_s` as
:math:`\overline{\partial \Omega_f}`, and
that part of the boundary of :math:`\Omega_f` which is shared by :math:`\Omega_s`.
In addition, :math:`\partial \Omega_{s}, \overline{\partial \Omega_s}` are analogously defined.
These distinct portions of the
domain boundary are illustrated in :numref:`fig-walls`.
The restrictions on the domain for Nek5000 are itemized below.

- The domain :math:`\Omega=\Omega_f \cup \Omega_s` must correspond either to a
  planar (Cartesian) two-dimensional geometry, or to the
  cross-section of an axisymmetric region specified by
  revolution of the cross-section about a specified axis, or
  by a (Cartesian) three-dimensional geometry.
- For two-dimensional and axisymmetric geometries, the
  boundaries of both subdomains, :math:`\partial \Omega_f` and
  :math:`\partial \Omega_s`, must be
  representable as (or at least approximated by) the union of
  straight line segments, splines, or circular arcs.
- Nek5000 can interpret a two-dimensional image as either
  a planar Cartesian geometry, or
  the cross-section of an axisymmetric body. In the case of
  the latter, it is assumed that the :math:`y`-direction is the radial
  direction, that is, the axis of revolution is at :math:`y=0`.
  Although an axisymmetric geometry is, in fact,
  three-dimensional, Nek5000 can assume that the field variables
  are also axisymmetric ( that is, do not depend on azimuth,
  but only :math:`y`, that is, radius, :math:`x`, and :math:`t` ), thus reducing the
  relevant equations to "two-dimensional" form.

Fully general three-dimensional meshes generated by other softwares
packages can be input to ``prenek`` as import meshes.
