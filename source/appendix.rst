==========
Appendices
==========

-----------------------------------
rea File (legacy)
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

-------------------------------
Internal Input Parameters/Switches
-------------------------------

....................
Parameters
....................

.. raw:: html

    <style> .red {color:red} </style>

.. role:: red

| **P001  DENSITY** density for the case of constant properties (for variable density see parameter ``P030``).
|
| **P002  VISCOS**  dynamic viscosity (if :math:`<0 \rightarrow Re`, otherwise :math:`1/Re`).
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
  | ``vtrans``    | (lx1,ly1,lz1,lelt,ldimt1)          | real    | convective coefficient                   |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vdiff``     | (lx1,ly1,lz1,lelt,ldimt1)          | real    | diffusion coefficient                    |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vxlag``     | (lx1,ly1,lz1,lelv,2)               | real    | x-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vylag``     | (lx1,ly1,lz1,lelv,2)               | real    | y-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``vzlag``     | (lx1,ly1,lz1,lelv,2)               | real    | z-velocity at previous time steps        |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``prlag``     | (lx2,ly2,lz2,lelv,lorder2)         | real    | pressure at previous time steps          |
  +---------------+------------------------------------+---------+------------------------------------------+
  | ``tlag``      | (lx1,ly1,lz1,lelv,lorder-1,ldimt1) | real    | t and ps at previous time steps          |
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
  | ``cbc``       | (6,lelt,ldimt1)           | character*3 | boundary condition                        |
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


