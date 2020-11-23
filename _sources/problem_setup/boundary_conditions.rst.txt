-------------------------------
Boundary Conditions
-------------------------------

.. The boundary conditions for Nek5000 are stored as part of the mesh, i.e. either part of the ``.rea`` or ``.re2`` file.
.. Any mesh generated with either ``prenek`` or ``genbox`` will include the assigned boundary conditions.
.. These are available at runtime in the ``cbc(iface,iel,ifld)`` array, indexed by face number, local element number, and field number.
.. For meshes converted from exodus format via the ``exo2nek`` script, the sideset numbers will be converted.
.. These are available at runtime in the ``bc(5,iface,iel,1)`` array, indexed by face number and local element number.
.. All sidesets will need to be translated into appropriate boundary conditions.
.. It is recommended to do this in ``usrdat``.
.. The available boundary conditions for velocity are listed in :numref:`tab:BCf`, and for temperature and passive scalars in :numref:`tab:BCt`.
.. 

The boundary conditions can be imposed in various ways:

- when the mesh is generated, e.g. with ``genbox``, as is explained in :ref:`sec:genbox`
- when an ``.rea`` file is read in ``prenek``
- translated from side-set numbers in ``usrdat`` when using ``exotonek`` or similar

.. TODO: add exotonek tutorial

The general convention for boundary conditions is

- uppercase letters correspond to primitive boundary conditions, as given in :numref:`tab:BCf`, :numref:`tab:BCt`
- lowercase letters correspond to user defined boundary conditions, see :numref:`tab:uBCf` , :numref:`tab:userBCt`
- lowercase letters ending with ``l``, i.e. ``'vl '``, are specified in face-local coordinates, i.e. normal, tangent and bitangent directions.

Uppercase boundary conditions which require assigned values in the ``.rea`` file are considered legacy and are not recommended for use.

..............
Fluid Velocity
..............

Two types of boundary conditions are applicable to the fluid velocity : essential (Dirichlet) boundary condition in which the velocity is specified and natural (Neumann) boundary condition in which the traction is specified.
For segments that constitute the boundary :math:`\partial \Omega_f`, see :numref:`fig-walls`, one of these two types of boundary conditions must be assigned to each component of the fluid velocity.
The fluid boundary condition can be *all Dirichlet* if all velocity components of :math:`{\bf u}` are specified; or it can be *all Neumann* if all traction components :math:`{\bf t} = [-p {\bf I} + \mu (\nabla {\bf u} + (\nabla {\bf u})^{T})] \cdot {\bf n}`, where :math:`{\bf I}` is the identity tensor, :math:`{\bf n}` is the unit normal and :math:`\mu` is the dynamic viscosity, are specified; or it can be *mixed Dirichlet/Neumann* if Dirichlet and Neumann conditions are selected for different velocity components.
If the nonstress formulation is selected, then traction is not defined on the boundary.
In this case, any Neumann boundary condition imposed must be homogeneous, i.e. equal to zero.
.. In addition, mixed Dirichlet/Neumann boundaries must be aligned with one of the Cartesian axes.
For flow geometry which consists of a periodic repetition of a particular geometric unit, the periodic boundary conditions can be imposed, as illustrated in :numref:`fig-walls` .
The available primitive boundary conditions for the fluid are given in :numref:`tab:BCf` , with the user-specified boundary conditions in :numref:`tab:uBCf` .

.. _tab:BCf:

.. table:: Primitive boundary conditions for velocity

   +------------+----------------------------+--------------+---------------------------------------------------+
   | Identifier | Description                | Type         | Note                                              |
   +============+============================+==============+===================================================+
   | P          | Periodic                   | --           | Standard periodic boundary condition              |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | p          | Periodic                   | --           | For periodicity within a single element           |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | O          | Outflow                    | Neumann      | Open boundary condition, zero pressure            |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | ON         | Outflow, Normal            | Mixed        | Zero velocity in non-normal directions            |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | W          | Wall                       | Dirichlet    | No slip, :math:`{ \bf{u} = 0}`                    | 
   +------------+----------------------------+--------------+---------------------------------------------------+
   | SYM        | Symmetry                   | Mixed        |                                                   | 
   +------------+----------------------------+--------------+---------------------------------------------------+
   | A          | Axisymmetric boundary      | --           |                                                   |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | E          | Interior boundary          | --           |                                                   |
   +------------+----------------------------+--------------+---------------------------------------------------+
   
.. _tab:uBCf:

.. table:: User defined boundary conditions for velocity

   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | Identifier | Description                 | Type         | Note                                                                 |
   +============+=============================+==============+======================================================================+
   | v          | Velocity                    | Dirichlet    | Standard velocity boundary condition                                 |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | vl         | Velocity, local             | Dirichlet    | Face-local coordinates (normal, tangnent, bitangent)                 |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | o          | Outflow                     | Neumann      | Open boundary condition, specified pressure                          |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | on         | Outflow, Normal             | Mixed        | Zero velocity in non-normal directions                               |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | s          | Traction                    | Neumann      | Specified traction in all directions                                 |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | sl         | Traction, local             | Neumann      | Face-local coordinates (normal, tangent, bitangent)                  | 
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | sh         | Traction, horizontal        | Mixed        | Specified traction with zero normal velocity                         |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | shl        | Traction, horizontal, local | Mixed        | Zero normal velocity, traction in tangent and bitangent              |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | int        | Interpolated (NEKNEK)       | Dirichlet    | Interpolated from the adjacent overset mesh, see: :ref:`neknek`      |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | mm         | Moving mesh                 | --           |                                                                      | 
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | ms         | Moving surface              | --           |                                                                      |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | msi        | Moving internal surface     | --           |                                                                      |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+   
   | mv         | Moving boundary             | Dirichlet    |                                                                      |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | mvn        | Moving boundary, normal     | Dirichlet    | Zero velocity in non-normal directions                               |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+ 
.. | ms         | Moving surface              | --           |                                                         |
.. +------------+-----------------------------+--------------+---------------------------------------------------------+
.. | msi        | Moving internal surface     | --           |                                                         |
.. +------------+-----------------------------+--------------+---------------------------------------------------------+    

.. _tab:LBCf:

.. .. table:: Legacy boundary conditions for velocity
..
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | Identifier | Description             | Type            | Note                                    |
..   +============+=========================+=================+=========================================+
..   | V          | Velocity                | Dirichlet       |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | VL         | Velocity, local         | Dirichlet       |                                         | 
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | S          | Traction                | Neumann         |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | SL         | Traction, local         | Neumann         |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | MM         | Moving mesh             | Dirichlet       |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | MS         | Moving surface          | Dirichlet       |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
..   | MSI        | Moving interior surface | --              |                                         |
..   +------------+-------------------------+-----------------+-----------------------------------------+
.. | MF         |                         | --              |                                         |
.. +------------+-------------------------+-----------------+-----------------------------------------+
.. | WS         |                         | --              |                                         |
.. +------------+-------------------------+-----------------+-----------------------------------------+
.. | WSL        |                         | --              |                                         |
.. +------------+-------------------------+-----------------+-----------------------------------------+


The open(outflow) boundary condition ("O") arises as a natural boundary condition from the variational formulation of Navier Stokes. 
We identify two situations

- In the non-stress formulation, open boundary condition ('Do nothing')

  .. math::

     [-p{\bf I} + \nu(\nabla {\bf u})]\cdot {\bf n}=0

- In the stress formulation, free traction boundary condition

  .. math::

     [-p{\bf I} + \nu(\nabla {\bf u}+\nabla {\bf u}^T)]\cdot {\bf n}=0

- the symmetric boundary condition ("SYM") is given as

  .. math::

     {\bf u} \cdot {\bf n} &= 0\ ,\\
     (\nabla {\bf u} \cdot {\bf t})\cdot {\bf n} &= 0

  where :math:`{\bf n}` is the normal vector and :math:`{\bf t}` the tangent vector. 
  If the normal and tangent vector are not aligned with the mesh the stress formulation has to be used.
- the periodic boundary condition ("P") needs to be prescribed in the ``.rea`` or ``.re2`` file since it already assigns the last point to first via :math:`{\bf u}({\bf x})={\bf u}({\bf x} + L)`, where :math:`L` is the periodic length.
- the wall boundary condition ("W") corresponds to :math:`{\bf u}=0`.

For a fully-developed flow in such a configuration, one can effect great computational efficiencies by considering the problem in a single geometric unit (here taken to be of length :math:`L`), and requiring periodicity of the field variables.
Nek5000 requires that the pairs of sides (or faces, in the case of a three-dimensional mesh) identified as periodic be identical (i.e., that the geometry be periodic).

For an axisymmetric flow geometry, the axis boundary condition is provided for boundary segments that lie entirely on the axis of symmetry.
This is essentially a symmetry (mixed Dirichlet/Neumann) boundary conditionin which the normal velocity and the tangential traction are set to zero.

For free-surface boundary segments, the inhomogeneous traction boundary conditions involve both the surface tension coefficient :math:`\sigma` and the mean curvature of the free surface.

...............................
Temperature and Passive Scalars
...............................

The three types of boundary conditions applicable to the
temperature are: essential (Dirichlet) boundary
condition in which the temperature is specified;
natural (Neumann) boundary condition in which the heat flux
is specified; and mixed (Robin) boundary condition
in which the heat flux is dependent on the temperature
on the boundary.
For segments that constitute the boundary
:math:`\partial \Omega_f' \cup \partial \Omega_s'` (refer to Fig. 2.1),
one of the above three types of boundary conditions must be
assigned to the temperature.

The two types of Robin boundary condition for temperature
are: convection boundary conditions for which the heat
flux into the domain depends on the heat transfer coefficient
:math:`h_{c}` and the difference between the environmental temperature
:math:`T_{\infty}` and the surface temperature; and radiation
boundary conditions for which the heat flux into the domain
depends on the Stefan-Boltzmann constant/view-factor
product :math:`h_{rad}` and the difference between the fourth power
of the environmental temperature :math:`T_{\infty}` and the fourth
power of the surface temperature.

The boundary conditions for the passive scalar fields are analogous to those used for the temperature field.
Thus, the temperature boundary condition menu will reappear for each passive scalar field so that the user can specify an independent set of boundary conditions for each passive scalar field.

.. _tab:BCt:

.. table:: Primitive boundary conditions (Temperature and Passive scalars)

   +------------+----------------------------+--------------+---------------------------------------------------+
   | Identifier | Description                | Type         | Note                                              |
   +============+============================+==============+===================================================+
   | P          | Periodic                   | --           | Standard periodic boundary condition              |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | p          | Periodic                   | --           | For periodicity within a single element           |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | I          | Insolated                  | Neumann      | zero gradient                                     |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | O          | Outflow                    | Neumann      | Identical to "I"                                  |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | SYM        | Symmetry                   | Neumann      | Identical to "I"                                  |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | A          | Axisymmetric boundary      | --           |                                                   |
   +------------+----------------------------+--------------+---------------------------------------------------+
   | E          | Interior boundary          | --           |                                                   |
   +------------+----------------------------+--------------+---------------------------------------------------+

.. _tab:userBCt:

.. table:: User defined boundary conditions for temperature and passive scalars

   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | Identifier | Description                 | Type         | Note                                                                 |
   +============+=============================+==============+======================================================================+
   | t          | Temperature                 | Dirichlet    | Standard Dirichlet boundary condition                                |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | f          | Flux                        | Neumann      | Standard Neumann boundary condition                                  |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | c          | Newton cooling              | Robin        | Specified heat transfer coefficient                                  |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
   | int        | Interpolated (NEKNEK)       | Dirichlet    | Interpolated from the adjacent overset mesh, see: :ref:`neknek`      |
   +------------+-----------------------------+--------------+----------------------------------------------------------------------+
  
- open boundary condition ("O")

  .. math::

     k(\nabla T)\cdot {\bf n} =0

- insulated boundary condition ("I")

  .. math::

     k(\nabla T)\cdot {\bf n} =0

  where :math:`{\bf n}` is the normal vector and :math:`{\bf t}` the tangent vector. If the normal and tangent vector are not aligned with the mesh the stress formulation has to be used.
- the periodic boundary condition ("P") needs to be prescribed in the ``.rea`` file since it already assigns the last point to first via :math:`{\bf u}({\bf x})={\bf u}({\bf x} + L)`, where :math:`L` is the periodic length.
- Newton cooling boundary condition ("c")

  .. math::

     k(\nabla T)\cdot {\bf n}=h(T-T_{\infty})

- flux boundary condition ("f")

  .. math::

     k(\nabla T)\cdot {\bf n} =f


............................
Internal Boundary Conditions
............................

In the spatial discretization, the entire computational
domain is subdivided into macro-elements, the boundary
segments shared by any two of these macro-elements
in :math:`\Omega_f` and :math:`\Omega_s` are denoted as internal boundaries.
For fluid flow analysis with a single-fluid system or heat
transfer analysis without change-of-phase, internal
boundary conditions are irrelevant as the corresponding
field variables on these segments are part of the
solution. However, for a multi-fluid system and for
heat transfer analysis with change-of-phase, special
conditions are required at particular internal
boundaries, as described in the following.

For a fluid system composes of multiple immiscible fluids,
the boundary (and hence the identity) of each fluid must
be tracked, and a jump in the normal traction exists
at the fluid-fluid interface if the surface tension
coefficient is nonzero.
For this purpose, the interface between any two fluids
of different identity must be defined as a special type of
internal boundary, namely, a fluid layer;
and the associated surface tension coefficient also
needs to be specified.

In a heat transfer analysis with change-of-phase, Nek5000 assumes
that both phases exist at the start of the solution, and that
all solid-liquid interfaces are specified as special internal
boundaries, namely, the melting fronts.
If the fluid flow problem is considered, i.e., the energy
equation is solved in conjunction with the momentum and
continuity equations, then only
the common boundary between the fluid and the solid
(i.e., all or portion of :math:`\partial \overline{\Omega}_f'` in :numref:`fig-walls`)
can be defined as the melting front.
In this case, segments on :math:`\partial \overline{\Omega}_f'` that
belong to the dynamic melting/freezing interface need to be
specified by the user.
Nek5000 always assumes that the density of the two phases
are the same (i.e., no Stefan flow); therefore at the melting
front, the boundary condition for the fluid velocity is the
same as that for a stationary wall, that is, all velocity
components are zero.
If no fluid flow is considered, i.e., only the energy equation
is solved, then any internal boundary can be defined as
a melting front.
The temperature boundary condition at the melting front
corresponds to a Dirichlet
condition; that is, the entire segment maintains a constant temperature
equal to the user-specified melting temperature :math:`T_{melt}`
throughout the solution.
In addition, the volumetric latent heat of fusion :math:`\rho L`
for the two phases,
which is also assumed to be constant, should be specified.

