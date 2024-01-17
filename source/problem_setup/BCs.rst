.. _boundary-conditions:

-------------------------------
Boundary Conditions
-------------------------------

For all boundaries, it is necessary to provide *Nek5000* with the boundary condition *type* (e.g. that the velocity is to be specified), and for some types, the boundary condition *value* (e.g. what the velocity is specified as).
Boundary condition values are assigned in the ``.usr`` file in the ``userbc`` subroutine (see :ref:`sec:userbc` for details).
This section focuses on what boundary condition types are available in *Nek5000* and how to assign them.

We first discuss how *Nek5000* stores the boundary condition information in the code itself.
A boundary condition type is stored on every face of every element for every solution field.
These are stored in memory in the ``cbc`` array, which is of type ``character*3``.
It is declared in a common block in ``Nek5000/core/INPUT`` as ``cbc(6,lelt,0:ldimt1)``, which is then included in the user subroutines in the ``.usr`` file as well as any relevant source code files.
Each boundary condition type has a corresponding *character identifier code* which tells *Nek5000* how to treat the boundary on each element.
The ``cbc`` array is set at runtime either from information carried directly in the mesh file (``.rea``/``.re2``), from information provided in the ``.par`` file, or (for advanced users) directly in ``usrdat`` or ``usrdat2``.

The general conventions for boundary condition identifiers are:

- uppercase letters correspond to boundary conditions which do not require any additional user input.
- lowercase letters correspond to user defined boundary conditions which require values to be set in the :ref:`user file <sec:userbc>`.
- lowercase letters ending with ``l``, e.g. ``vl``, are specified in face-local coordinates, i.e. normal, tangent and bitangent directions.

The available boundary condition types, along with the identifier codes, are described in the following sections for the fluid velocity and pressure, and the temperature and passive scalars.

.. _sec:velbcs:

...........................
Fluid Velocity and Pressure
...........................

Two kinds of boundary conditions are applicable to the fluid velocity: essential (Dirichlet) boundary conditions in which the velocity is specified, and natural (Neumann) boundary conditions in which the traction is specified.
For segments that constitute the boundary :math:`\partial \Omega_f`, see :numref:`fig-walls`, one of these two types of boundary conditions must be assigned to each component of the fluid velocity.

The fluid boundary condition can be *all Dirichlet* if all velocity components of :math:`{\bf u}` are specified, or it can be *all Neumann* if all three traction components (:math:`\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}`) are specified on the boundary face. 
Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
It can also be *mixed Dirichlet/Neumann* if Dirichlet and Neumann conditions are selected for different velocity components.
If the :ref:`no-stress formulation <sec:nostress>` is selected, then traction is not defined on the boundary.
In this case, any Neumann boundary condition imposed must be homogeneous, i.e. equal to zero, and mixed Dirichlet/Neumann boundaries must be aligned with one of the Cartesian axes.
These conditions are not required for the :ref:`full-stress formulation <sec:fullstress>`.

.. For flow geometries which consist of a periodic repetition of a particular geometric unit, periodic boundary conditions can be imposed, as illustrated in :numref:`fig-walls` .

Neumann boundary conditions for velocity are assigned differently depending on if the no-stress or full-stress formulation is used.
In the no-stress formulation, the non-diagonal terms are neglected and we define the stress tensor as:

 .. math:: 

  \boldsymbol{\underline \tau} \equiv \mu \nabla \bf u

and in the full-stress formulation as:

 .. math::

   \boldsymbol{\underline \tau} \equiv \mu\left[\nabla {\bf u} + \left(\nabla {\bf u}\right)^T\right]
  
Inlet (Dirichlet), ``v``
````````````````````````

Standard Dirichlet boundary condition for velocity.

 .. math::

     \nabla p \cdot {\bf \hat e_n} &= \mu\left(\nabla\times\boldsymbol{\omega}\right) \cdot {\bf \hat e_n}\\
     {\bf u} \cdot {\bf \hat e_x} &= u_x\\
     {\bf u} \cdot {\bf \hat e_y} &= u_y\\
     {\bf u} \cdot {\bf \hat e_z} &= u_z
    
Where, :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`\boldsymbol{\omega}` is the vorticity, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`u_x`, :math:`u_y`, and :math:`u_z` are set in the :ref:`user file <sec:userbc>`.

Inlet (Dirichlet) - local ``vl``
````````````````````````````````

Standard Dirichlet boundary condition for velocity in local coordinates.

 .. math::

     \nabla p \cdot {\bf \hat e_n} &= \mu\left(\nabla\times\boldsymbol{\omega}\right) \cdot {\bf \hat e_n}\\
     {\bf u} \cdot {\bf \hat e_n} &= u_n\\
     {\bf u} \cdot {\bf \hat e_t} &= u_1\\
     {\bf u} \cdot {\bf \hat e_b} &= u_2
    
Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, :math:`\boldsymbol{\omega}` is the vorticity, and :math:`u_n`, :math:`u_1`, and :math:`u_2` are set in the :ref:`user file <sec:userbc>`.


Outlet, ``O``
`````````````

The open (outflow) boundary condition arises as a natural boundary condition from the variational formulation of Navier Stokes. 

  .. math::

     p &= 0\\
     \boldsymbol{\underline \tau} \cdot {\bf \hat e_n} &= 0

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face.
The ``userbc`` subroutine is not called for this boundary condition type.

Pressure Outlet, ``o``
``````````````````````

Similar to a standard outlet, but with a specified pressure.

  .. math::

     p &= p_a\\
     \boldsymbol{\underline \tau} \cdot {\bf \hat e_n} &= 0

Where :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`p_a` is set in the :ref:`user file <sec:userbc>`.

Outlet - Normal, ``ON``
```````````````````````

Open boundary with zero velocity in the tangent and bitangent directions.

  .. math::
     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} &= 0\\
     {\bf u} \cdot {\bf \hat e_t} &= 0\\
     {\bf u} \cdot {\bf \hat e_b} &= 0

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.
The ``userbc`` subroutine is not called for this boundary condition type.

Pressure Outlet - normal, ``on``
````````````````````````````````

Similar to an outlet - normal boundary, but with a specified pressure.

  .. math::

     p &= p_a\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right) \cdot {\bf \hat e_n} &= 0\\
     {\bf u} \cdot {\bf \hat e_t} &= 0\\
     {\bf u} \cdot {\bf \hat e_b} &= 0

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`p_a` is set in the :ref:`user file <sec:userbc>`.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.

Symmetry, ``SYM``
`````````````````

Symmetric face or a slip wall.

  .. math::

     \nabla p \cdot {\bf \hat e_n} &= 0\\
     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= 0

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face.
If the surface normal vector is not aligned with a principal Cartesian axis, the :ref:`full-stress formulation <sec:fullstress>` must be used.
The ``userbc`` subroutine is not called for this boundary condition type.

Traction, ``s``
```````````````

Full Neumann boundary conditions for velocity.

  .. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_x} &= tr_x\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_y} &= tr_y\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_z} &= tr_z

Where, :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face, :math:`{\bf \hat e_x}`, :math:`{\bf \hat e_y}`, and :math:`{\bf \hat e_z}` are unit vectors aligned with the Cartesian axes and :math:`tr_x`, :math:`tr_y`, and :math:`tr_z` are set in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for the this boundary type.

Traction - local, ``sl``
````````````````````````

Similar to traction, but in local coordinates.

  .. math::

     p &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_n} &= tr_n\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, and :math:`tr_n`, :math:`tr_1`, and :math:`tr_2` are set in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for the this boundary type.

Traction - horizontal, local, ``shl``
`````````````````````````````````````

Similar to symmetry, but with specified non-zero traction in the tangent and bitangent directions.

  .. math::

     \nabla p \cdot {\bf \hat e_n} &= \mu\left(\nabla\times\boldsymbol{\omega}\right) \cdot {\bf \hat e_n}\\
     {\bf u} \cdot {\bf \hat e_n} &= 0\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_t} &= tr_1\\
     \left(\boldsymbol{\underline \tau} \cdot {\bf \hat e_n}\right)\cdot {\bf \hat e_b} &= tr_2

Where, :math:`{\bf \hat e_n}`, :math:`{\bf \hat e_t}`, and :math:`{\bf \hat e_b}` are the normal, tangent, and bitangent unit vectors on the boundary face, :math:`\boldsymbol{\omega}` is the vorticity, and :math:`tr_1` and :math:`tr_2` are set in the :ref:`user file <sec:userbc>`.
The :ref:`full-stress formulation <sec:fullstress>` must be used for the this boundary type.

Wall, ``W``
```````````

Dirichlet boundary condition corresponding to a no-slip wall.

  .. math::

     \nabla p \cdot {\bf \hat e_n} &= \mu\left(\nabla\times\boldsymbol{\omega}\right) \cdot {\bf \hat e_n}\\
     {\bf u} &= 0

Where, :math:`{\bf \hat e_n}` is the unit vector normal to the boundary face and :math:`\boldsymbol{\omega}` is the vorticity.
The ``userbc`` subroutine is not called for this boundary condition type.
  
Other BCs
`````````

.. _tab:BCf:

.. csv-table:: Primitive boundary conditions for velocity
   :header: Identifier,Description,Type,Note
   :widths: 5,30,10,55

   ``P`` , "Periodic", --, "Standard periodic boundary condition"
   ``p`` , "Periodic", --, "For periodicity within a single element"
   ``A`` , "Axisymmetric boundary", --, "Can only be used on face 1, treated as ``SYM``"
   ``E`` , "Interior boundary", --, "--"
   ``'   '`` , "Empty", --, "Treated as an interior boundary"
   ``sh`` , "Traction, horizontal",        Mixed,     "Specified traction with zero normal velocity"
   ``int``, "Interpolated (NEKNEK)",       Dirichlet, "Interpolated from the adjacent overset mesh, see: :ref:`neknek`"
   ``mm`` , "Moving mesh",                 --,        "--"
   ``ms`` , "Moving surface",              --,        "--"
   ``msi``, "Moving internal surface",     --,        "--"
   ``mv`` , "Moving boundary",             Dirichlet, "--"
   ``mvn``, "Moving boundary, normal",     Dirichlet, "Zero velocity in non-normal directions"

.. Note::

   To use periodic boundary conditions, ``P``, in third-party meshes the face meshes must be conformal and must have a corresponding pair of boundary ID values which need to be provided during conversion, i.e. to ``exo2nek``, ``gmsh2nek``, or ``cgns2nek``. 
   Additionally, the mesh must be at least 3 elements thick in the direction normal to the periodic boundaries.
 
The periodic boundary condition (``P``) needs to be prescribed in the ``.rea`` or ``.re2`` file since it already assigns the last point to first via :math:`{\bf u}({\bf x})={\bf u}({\bf x} + L)`, where :math:`L` is the periodic length. 
For a fully-developed flow in such a configuration, one can effect great computational efficiencies by considering the problem in a single geometric unit (here taken to be of length :math:`L`), and requiring periodicity of the field variables. 
*Nek5000* requires that the pairs of sides (or faces, in the case of a three-dimensional mesh) identified as periodic be identical (i.e., that the geometry be periodic).

For an axisymmetric flow geometry, the axis boundary condition (``A``) is provided for boundary segments that lie entirely on the axis of symmetry. This is essentially a symmetry (mixed Dirichlet/Neumann) boundary condition in which the normal velocity and the tangential traction are set to zero.
This requires a 2D mesh where the x-axis is the axis of rotation.

For free-surface boundary segments, the inhomogeneous traction boundary conditions involve both the surface tension coefficient :math:`\sigma` and the mean curvature of the free surface.

.. _sec:tempbcs:

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

.. csv-table:: Primitive boundary conditions (Temperature and Passive scalars)
   :widths: 5,10,10,75
   :header: Identifier,Description,Type,Note

   ``P``, Periodic, --, "Standard periodic boundary condition"
   ``p``, Periodic, --, "For periodicity within a single element"
   ``I``, Insulated, Neumann, "zero gradient"
   ``O``, Outflow, Neumann, "Identical to ``I``"
   ``SYM``, Symmetry, Neumann, "Identical to ``I``"
   ``A``, Axisymmetric boundary, --, "treated as ``I``"
   ``E``, Interior boundary, --, "--"

.. _tab:userBCt:

.. csv-table:: User defined boundary conditions for temperature and passive scalars
   :widths: 5,10,10,75
   :header: Identifier,Description,Type,Note

   ``t``, "Temperature", "Dirichlet", "Standard Dirichlet boundary condition"
   ``f``, "Flux", "Neumann", "Standard Neumann boundary condition"
   ``c``, "Newton cooling", "Robin", "Specified heat transfer coefficient"
   ``int``, "Interpolated (NEKNEK)", "Dirichlet", "Interpolated from the adjacent overset mesh, see: :ref:`neknek`"
  
- open boundary condition ("O")

  .. math::

     k(\nabla T)\cdot {\bf \hat e_n} =0

- insulated boundary condition ("I")

  .. math::

     k(\nabla T)\cdot {\bf \hat e_n} =0

where :math:`{\bf \hat e_n}` is the unit normal vector, :math:`{\bf \hat e_t}` the unit tangent vector and :math:`{\bf \hat e_b}` is the unit bitangent vector.
If the normal, tangent, and bitangent vectors are not aligned with the mesh the stress formulation has to be used.
- the periodic boundary condition ("P") needs to be prescribed in the ``.rea`` file since it already assigns the last point to first via :math:`{\bf u}({\bf x})={\bf u}({\bf x} + L)`, where :math:`L` is the periodic length.
- Newton cooling boundary condition ("c")

  .. math::

     k(\nabla T)\cdot {\bf \hat e_n}=h(T-T_{\infty})

- flux boundary condition ("f")

  .. math::

     k(\nabla T)\cdot {\bf \hat e_n} =f


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

.. _sec:settingbcs:

..........................................................
Setting Boundary Conditions Types
..........................................................

Assigning boundary condition types in *Nek5000* is handled differently depending on if you are using a third-party meshing tool such as *Gmsh*, *ICEM*, *Cubit*, etc. and importing the mesh with ``exo2nek``, ``gmsh2nek``, or ``cgns2nek``, or if you are using a Nek-native tool such as *preNek* or ``genbox`` (see :ref:`tools_genbox`).
In either case, the boundary condition types are set by assigning the corresponding character identifier code in the character boundary condition array, ``cbc``.
The character boundary condition array itself is described :ref:`here <sec:probvars>` and the supported character codes were described in the sections above for :ref:`momentum <sec:velbcs>` and :ref:`temperature and passive scalars <sec:tempbcs>`.
The differences between Nek-native tools and third-party meshing tools are only in how this array gets set.
For Nek-native tools, this array is read directly from the ``.rea`` or ``.re2`` file, which is set based on input provided to the tool itself.
For third-party meshing tools, the boundary *ID* is set in the tool -- e.g. as a *sideset ID* in *ICEM* -- and this information is propagated to the ``.re2`` (mesh) file.
The ``cbc`` array is later filled at runtime based on the boundary IDs.

The recommended method of setting the boundary condition type from the boundary ID is through the ``.par`` file.
This is done through the ``boundaryTypeMap`` key, which is available for the ``VELOCITY``, ``TEMPERATURE``, and ``SCALARXX`` directives.
By default, *Nek5000* assumes the boundary IDs are sequential and start from 1.
If this is not the case, the optional ``boundaryIDMap`` key is available for the ``MESH`` directive.
See :ref:`here <case_files_par>` for more information on the ``.par`` file.
A few simple examples of setting the BC types via the ``.par`` file for a mesh with boundary IDs assigned in a third-party mesher are below.

.. warning::

   Setting the boundary condition types in the ``.par`` file is **NOT** supported in V19 or earlier versions. 

In the simplest example, the mesh has 4 boundaries each with a sequentially numbered boundary ID.

.. csv-table:: Desired Boundary Types
   :header: Boundary ID, Velocity, Temperature

   1,``v``,``t``
   2,``O``,``I``
   3,``W``,``f``
   4,``SYM``,``I``

To set the boundary condition types, the ``boundaryTypeMap`` key is used in the ``.par`` file.
The ``boundaryTypeMap`` key is a comma-separated list of the boundary condition types to be assigned to the domain and is avaialble for the velocity, temperature and passive scalar fields.
The character identifiers can always be used for assignment.
Additionally, some of the common boundary types can be assigned using plain-English equivalents in the ``.par`` file only.
For a list of these see :ref:`here <sec:engidentifiers>`.
By default, *Nek5000* assumes the boundary IDs in your mesh start with 1 and are numbered sequentially.
Due to the sequential ordering of the boundary IDs in this example, these boundary types can be set using only the ``boundaryTypeMap`` keys in the ``VELOCITY`` and ``TEMPERATURE`` directives:

.. code-block:: ini

   [VELOCITY]
   boundaryTypeMap = v, O, W, SYM

   [TEMPERATURE]
   boundaryTypeMap = t, I, f, I  

If your boundary IDs are not sequential or do not start with 1, they can be explicitly declared using the ``boundaryIDMap`` key in the ``MESH`` directive.
The ``boundaryIDMap`` key is a comma-separated list of integers corresponding to the boundary IDs in your mesh.
When using the ``boundaryIDMap`` key, *Nek5000* makes no assumptions regarding the boundary ID values.

.. code-block:: ini

   [MESH]
   boundaryIDMap = 3, 4, 1, 2

   [VELOCITY]
   boundaryTypeMap = W, SYM, v, O  

   [TEMPERATURE]
   boundaryTypeMap = f, I, t, I

