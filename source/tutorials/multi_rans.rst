.. _multi_rans:

-----------------------------------
Multi-Component RANS 
-----------------------------------

This is an advanced tutorial comprising RANS simulation in a multi-component geometry. 
The components include an inlet nozzle coupled with a wire-pin bundle. 
Given the scale of the geometry and mesh used in this  tutorial, it is assumed that the user has access to a computing cluster for running the case. 
The NekNek module is used for interfacing between the two components, allowing for runtime transfer of field data and coupled simultaneous solution of the flow and temperature equations. 
The tutorial employs several advanced  Nek5000 tools, procedures and user routines including

 * ``exo2nek`` for importing third party mesh file to Nek5000
 * ``reatore2`` for converting ASCII Nek5000 mesh to binary format
 * ``gencon`` for generating mesh connectivity file
 * ``wiremesher`` for generating wire-wrapped pin bundle mesh
 * ``nekamg_setup`` for generating algebraic multi-grid solver (AMG) files
 * :math:`k-\tau` RANS simulation setup 
 * ``NekNek`` setup

..........................
Before You Begin
..........................

It is highly recommended that new users familiarize themselves with the basic Nek5000 simulation
setup files and procedures outlined in the :ref:`fdlf` and :ref:`perhill` tutorials. Further, some 
sections of the tutorial will assume familiarity with :ref:`tutorial_rans` and :ref:`neknek` tutorials.

.. Warning::
  The RANS models rely on setting an implicit source term for robustness. This is not supported in V19 or earlier versions of *Nek5000*. Any implementation of RANS models should use the latest master branch from github.

..........................
Case Overview
..........................

All required files for this tutorial can be downloaded using the following link:

* :download:`inlet_bundle.tar <multi_rans/inlet_bundle.tar.gz>`

Extract the above tar file and navigate to the ``inlet_bundle`` (parent) directory

.. code-block:: console

 $ tar -xzvf inlet_bundle.tar.gz
 $ cd inlet_bundle
	
It has the following directories and files, discussed in detail in the following sections:

* ``inletMesh`` --> Directory containing inlet mesh related files.
* ``bundleMesh`` --> Directory containing pin-wire bundle mesh related files.
* ``SIZE`` --> Parameter file for defining problem size.
* ``inlet_bundle.usr`` --> Fortran file for user specified subroutines.
* ``inlet.par`` --> Nek5000 parameter file.
* ``bundle.par`` --> Nek5000 parameter file.
* ``SPLINE`` --> Common block Fortran file for spline interpolation.
* ``InletProf.dat`` --> Inlet condition data file.
* ``neksaw`` --> Sample script for launching NekNek job (on Sawtooth cluster).

..............................
Geometry Description
..............................

The geometry under consideration consists of two components including an inlet nozzle connected to a wire-wrapped pin bundle.
These are shown below with the inlet nozzle domain in red and the pin bundle domain in light gray. 

.. _fig:sfr_geom:

.. figure:: multi_rans/geom.png
   :align: center
   :figclass: align-center

   Component geometries including an inlet nozzle (red) connected to a wire-wrapped pin bundle (gray)
   
The geometry is non-dimensionalized with respect to the pin diameter, :math:`D`. 
The axial extent of the inlet nozzle component is :math:`z/D=[-4.75,0.25]` and the wire-pin bundle axial dimensions
are :math:`z/D=[0,12.5]`.  The two components, therefore, have an axial overlap of :math:`\Delta z/D= 0.25`,
while the lateral dimensions are conformal. 

.. Note::
  A reasonable overlap in the computational domains is imperative for a robust ``NekNek`` simulation.
  As a general rule, increasing the extent of overlap will render more stability to the coupled solver.

The Pitch-to-diameter ratio of the wire-pin bundle is 1.13 and the wire diameter is :math:`D_w=0.12875D`. 
To prevent sharp corners, a fillet of diameter :math:`0.05 D` is introduced between the pin and the wire and the wire is slightly submerged in the pin by a distance of :math:`0.025D_w`.
The side length of the hexagonal bundle is :math:`1.865D` and the diameter of the inlet nozzle boundary is :math:`1.875D`.
The length of the bundle component is equal to the wire pitch, :math:`12.5 D`.

.....................................................
Thermal-hydraulic Parameters and Boundary Conditions
.....................................................

All flow and thermal properties are listed in :numref:`tab:nb_BC`. 
A non-dimensional bulk velocity of  :math:`U_{in}=1` is specified at the inlet and inlet non-dimensional diameter is :math:`D_{in}=1.875 D`. 
The hydraulic diameter of the pin-wire bundle is :math:`D_b\approx 0.398 D`. 
A no-slip boundary condition is  specified on all walls and the corresponding values of :math:`k` and :math:`\tau` on the walls is zero. 
A constant non-dimensional temperature :math:`T^*=1` is specified at the inlet and a non-dimensional heat flux of :math:`q'' =1` is applied on the pin walls from half axial length of the bundle. 
An insulated wall condition is specified on all remaining boundaries. 
A Prandtl number of 0.005 was chosen corresponding to liquid Sodium.

.. _tab:nb_BC:

.. csv-table:: Flow and Thermal Properties 
   :header: Parameter, Variable, Value

   Inlet Reynolds number, :math:`Re_{in}`, 60000                           
   Bundle Reynolds number, :math:`Re_{b}`, 10000 (approx)                  
   Dimensionless Kinematic visocity, :math:`\nu^*`, 3.125e-5                        
   Prandtl Number, :math:`Pr`, 0.005                           
   Peclet Number, :math:`Pe`, 160                             
   Turbulent Prandtl number, :math:`Pr_t`, 1.5                             

..............................
Mesh Generation
.............................. 

########################
Inlet Nozzle (exo2nek)
########################

A third-party meshing tool (e.g., ANSYS ICEM) is required for generating the mesh for the inlet component.
In this case, the mesh was saved as an ``EXODUS II (.exo)`` mesh file. 
*Nek5000* offers the ``exo2nek`` mesh conversion tool for converting an ``.exo`` mesh file to ``.re2`` format, as well as ``gmsh2nek`` for converting ``.msh`` files from Gmsh.

By default, all Nektools support only 150,000 elements. 
To run this case, ``exo2nek`` and ``gencon`` will need to be compiled with support for more than 1,000,000 elements. 
To do this, go to the ``Nek5000/tools`` directory and edit the ``maketools`` file.

.. code-block:: console

  $ cd ~/Nek5000/tools
  $ vi maketools

Uncomment the line specifying the maximum number of elements and change it to a large value, such as that shown on the highlighted line below.

.. literalinclude:: multi_rans/maketools
   :emphasize-lines: 4

The ``exo2nek`` and ``gencon`` tools can then be recompiled with

.. code-block:: console

  $ ./maketools exo2nek gencon

.. Note::

  If you still get errors about the tools not supporting enough elements, you may need to delete the binary in ``Nek5000/bin`` and corresponding ``*.o`` files in the tool's subdirectory manually before recompiling.
  Executing ``$ ./makenek clean`` will also work, but this will delete all of the tools.

The above command will compile the ``exo2nek`` and ``genmap`` tools. 
The latter is required for :ref:`generating the mesh connectivity <multi_gencon>` file at a later step.  


Currently, ``exo2nek`` supports the following mesh elements,

 * ``TET4``
 * ``WEDGE6``
 * ``HEX8``
 * ``HEX20``

The user must ensure that the third-party mesh comprises only the above listed element types. Navigate to the folder 
containing inlet mesh and run ``exo2nek``.

.. code-block:: console
	
  $ cd inletMesh
  $ exo2nek

The output from ``exo2nek`` is shown below, where the expected user input is highlighted

.. literalinclude:: multi_rans/exo2nek.output
   :language: none
   :emphasize-lines: 2,4,28,59,61
	
Following the above steps will generate the file ``inlet.re2`` in the current directory. Note that the 
``sideSet ID`` for all mesh boundaries must be specified in the ``.exo`` file using the third-party 
meshing software of user's choice. For the inlet nozzle component, the following IDs are assigned:

 * 2 --> Nozzle inlet
 * 3 --> Interfacing surface (bundle surface)
 * 4 --> Walls

Return and move the mesh file to the parent directory:

.. code-block:: console

   $ cd ../
   $ mv inletMesh/inlet.re2 .

#############################
Wire-pin Bundle (wiremesher)
#############################

To generate the pin-wire bundle mesh, navigate to the ``wireMesh`` folder:

.. code-block:: console

   $ cd wireMesh
	
It contains two sub-directories, viz., ``wire2nek`` and ``matlab``. Input parameters for the meshing 
script are specified in the header of the ``matlab/wire_mesher.m`` file, as shown:

.. literalinclude:: multi_rans/inlet_bundle/bundleMesh/matlab/wire_mesher.m
	:language: matlab
	:lines: 1-20
	
All input variables are suitably annotated in the above code snippet. Although all input dimensions shown are in ``mm``, the 
script eventually produces the wiremesh in non-dimensional units, normalized with pin diameter ``D``. It will usually
take some heuristic experimentation to specify optimum parameters based on user requirements (such as resolution,
number of elements) and to ensure that the mesh does not have Jacobian related errors. Critical parameters, that 
control the mesh resolution and contribute towards a successful mesh, include: 

 * ``Df`` --> Higher fillet diameter will be less likely to cause any errors and produce a smoother mesh. Should be adjusted to a reasonable value.
 * ``T`` --> Trims off a portion of the wire to avoid pinching between the wire and neighboring pin. Typically set to ``0.05*Dw``.
 * ``S`` --> Submerges the wire slightly into the pin to avoid sharp corners. Typically set to ``0.025*Dw``.
 * ``Adjust`` --> Ensures trimming only occurs if wire passes neighboring pins. Typically set to 1.
 * ``iFtF`` --> (Deprecated) Adds a layer next to outer wall. Set to 0.
 * ``FtF_rescale`` --> (Deprecated) Inactive if ``iFtF=0``. 
 * ``G`` --> Controls gap between peripheral wire and outer wall. Adjust to needed value.
 * ``ne`` --> controls the number of pins in the radial direction. ``ne=2`` will produce a bundle with 7 pins, ``ne=3`` will produce 19 pins, and so on. 
 * ``Col`` --> Controls the resolution in the azimuthal direction.
 * ``Row`` --> Controls the resolution in radial direction. 
 * ``Rowdist`` --> Controls layer width distribution percentage in radial direction, from interior to pin wall. Must add up to 100 and entries must be equal to ``Row``.
 * ``Lay`` --> Controls resolution in axial direction. Specifies number of elements in axial length equal to 60 degree rotation of wire.  

The mesher is initiated by simply running the ``doall.sh`` bash script. Ensure that both Matlab and Python with numpy
(tested with Python 3.8) are active before launching the script and that Fortran compilers are available. 

.. code-block:: console

   $ ./doall.sh
	
The script can take a while to complete. Upon completion it generates ``wire_out.rea`` mesh file, which is the 
ASCII mesh file for Nek5000. Convert this into the binary format by running ``reatore2`` tool. Follow the prompts:

.. code-block:: console

   Input .rea name:
   wire_out
	
.. code-block:: console

   Input .rea/.re2 output name:
   bundle
	
We finally obtain the ``bundle.re2`` file which contains the pin-wire bundle mesh for Nek5000 run. Boundary IDs
are assigned by the ``wiremesher`` as:

 * 1 --> Fuel pin walls
 * 2 --> Bundle hexagonal (outer) walls
 * 3 & 4 --> Axial end surfaces

Return and move the mesh file to the parent directory:

.. code-block:: console

   $ cd ../
   $ mv bundleMesh/bundle.re2 .


.. _multi_gencon:
 
####################################
Generating Connectivity file (.co2)
####################################

After generating the mesh files for both components, it is necessary to generate the corresponding connectivity files
using ``gencon`` tool. Note that using ``gencon`` tool instead of ``genmap``, which generates  map (``.ma2``) file, is the 
recommended procedure for large meshes. See :ref:`build_pplist` for details. 

Users must include ``PARRSB`` in the ``PPLIST`` in ``makenek`` file (location: Nek5000/bin/makenek) which instructs Nek5000
to partition the mesh during run-time and requires ``.co2`` file instead of ``.ma2`` for running the case. 

Run ``gencon`` from the parent folder for each mesh file. Users will be prompted to specify the mesh file name and tolerance.
Use 0.01 for inlet and 0.2 (default) for bundle mesh:

.. code-block:: console

	Input .rea / .re2 name:
	inlet
	reading inlet.re2                                                                   
	Input mesh tolerance (default 0.2):
	0.01
	
.. code-block:: console

	Input .rea / .re2 name:
	bundle
	reading bundle.re2                                                                   
	Input mesh tolerance (default 0.2):
	0.2
	
The above will generate ``inlet.co2`` and ``bundle.co2`` connectivity files, respectively.

.........................................
Parameter File (.par)
.........................................

``NEKNEK`` requires separate ``.par`` file for each of the components. The files are included in the parent folder and shown below:

.. literalinclude:: multi_rans/inlet_bundle/inlet.par
	:language: ini
	
.. literalinclude:: multi_rans/inlet_bundle/bundle.par
	:language: ini
	
Both parameter files are identical, except for one important difference. To restart the case from
any given time, separate restart file names should be specified to the ``startFrom`` parameter. It is critical that the
properties and time step size are identical for both ``.par`` files. Values are assigned in dimensionless form;
density is set to unity, viscosity and diffusivity are set to :math:`1/\nu^*` and conductivity to Peclet number 
(-ve sign indicates that solver is run in dimensionless form).

Note that given the large size of meshes, the ``preconditioner`` must be set to ``semg_amg_hypre``. This invokes the algebraic
multigrid (AMG) solver for pressure instead of the default ``XXT`` solver. The AMG preconditioner requires third party ``HYPRE``
libraries which must be included in the preprocessor list option in ``makenek`` file (location: Nek5000/bin/makenek). 
Also including the ``PARRSB`` option, as mentioned earlier, the ``PPLIST`` therefore is:

.. code-block:: console
	
	PPLIST = "HYPRE PARRSB"

Further details on all parameters of ``.par`` file can be found :ref:`here <case_files_par>`.
 
.........................................
User Routines (.usr file)
......................................... 

Basics of the required setup routines for a NekNek simulation can be found in the :ref:`neknek` turorial, while for a RANS simulation
in the :ref:`tutorial_rans` tutorial. Although this section decribes all user routines required for a NekNek RANS simulation in detail, 
a comprehensive understanding of routines from these simpler cases is recommended before proceeding.

Following headers are required at the beginning of ``.usr`` file for loading RANS related subroutines:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 17-18
	
``NekNek`` related parameters are specified in ``usrdat`` routine:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 220-238
		
``ngeom`` specifies the number of overlapping Schwarz-like iterations, while ``ninter`` controls the time 
extrapolation order of boundary conditions at the overlapping interface. ``ninter=1`` is unconditionally 
stable, while a higher temporal order will typically require more iterations for stability (``ngeom>2``). 
For computational savings, we maintain first order temporal extrapolation for this tutorial. 
``nfld_neknek`` specifies the number of total field arrays that are transferred between the two meshes 
and must be equal to 7 for 3D RANS cases (3 velocity, 1 pressure and 3 scalar field arrays - temperature,
:math:`k` and :math:`\tau`).

:Note:
	Ensure that proper common block headers are included in subroutines. ``NEKNEK`` header is required 
	for routines where ``idsess`` needs to be accessed, as shown below.
	
Boundary Condition specification and RANS initialization is performed in ``usrdat2``:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 240-320
	
``NekNek`` solver launches two Nek5000 sessions simultaneously and field data transfer is performed between the
two sessions on each time iteration. Each session is assigned a unique id, stored in the variable ``idsess``.
Here, ``idsess=0`` is assigned to the inlet component solve and ``idsess=1`` to bundle component. Boundary 
conditions are assigned using this variable for each component, as shown above. 

Recall the boundary IDs assigned to each component during the mesh generation process, described in the preceding section.
Character codes for different boundary conditions are stored in the ``cbc`` array. Their detailed description can be found in
:ref:`boundary-conditions`. For each component, the nested loops go through all elements and their faces, populating
``cbc`` array for all fields based on mesh assigned boundary IDs. Note that ``int`` boundary condition
must be assigned to the overlapping surfaces of the inlet and bundle components. ``int`` condition is replaced 
internally with Dirichlet boundary conditions subsequently by Nek5000. Flux boundary condition, ``f``, is assigned to
pin walls for temperature field while insulated, ``I``, is assigned to all other walls.

With regards to RANS initialization; ``m_id=4`` selects the :math:`k-\tau` RANS model and ``w_id=2`` selects the
wall distance computing algorithm. :math:`k` and :math:`\tau` fields are stored in the 3rd and 4th index, respectively, 
specified with ``ifld_k`` and ``ifld_omega``. Set ``ifcoeffs`` to ``.true.`` only if user specified RANS coefficients
are required. For details on the RANS related parameters, refer :ref:`tutorial_rans` tutorial.  

:Note:
	``rans_init`` must be called after populating ``cbc`` array

For RANS simulation, diffusion coefficients are assigned in the ``uservp`` routine. The routine used here remains 
nearly identical to the :ref:`tutorial_rans` tutorial:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 20-51
	
Only turbulent Prandtl number is changed to ``Pr_t=1.5`` for this tutorial. This value is more appropriate for 
molten sodium salts as compared to the default value of 0.85 (for air), which is assigned through ``rans_turbPrandtl()`` 
function call.

Source terms for the temperature and scalar equations are assigned through ``userq``. The routine here is identical
to the basic :ref:`tutorial_rans` case:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 76-103

Note that either component does not have any volumetric source heat source and hence ``qvol=0`` for temperature
field (``ifield.eq.2``).

Initial conditions are specified in ``useric``. Similar values are assigned to both components and, thus, the routine
implementation is straightforward. Temperature is initalized to 1 for both components.

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 180-203

Boundary conditions are assigned in ``userbc``. For the inlet component, inlet conditions are assigned using data
generated from RANS simulation in a pipe with identical diameter as the inlet surface. The data is stored in the 
``InletProf.dat`` file which contains axial velocity, :math:`k` and :math:`\tau` information as a function of
radial wall distance. Two plugin subroutines are required, which perform spline interpolation of the data to the 
inlet mesh, viz., ``getInletProf`` and ``init_prof``. These are provided for the user in the ``inlet_bundle.usr``
file and can be used without modification. The usage is shown below:

.. literalinclude:: multi_rans/inlet_bundle/inlet_bundle.usr
	:language: fortran
	:lines: 105-178
	
Note that the diameter of the inlet surface is ``din=1.875``. Spline interpolation routine, ``init_prof``, requires
the wall distance array, ``wd``, which is populated in ``usrdat2`` (in ``rans_init`` call). The distance should be
limited to inlet radius to avoid spline from extrapolating. Inlet component is identified with ``idsess.eq.0`` and
inlet surface with its boundary ID, ``id_face.eq.2``.

Temperature flux must also be assigned in ``userbc`` on the pin surface walls. As mentioned earlier, non-dimensional
unit heat flux is assigned from half axial length of the bundle (``zmid``). A smoothed axial flux profile is imposed using
``tanh`` step function as shown. Flux on all remaining walls is zero.

Field data transfer between ``NekNek`` sessions is performed using ``valint`` array and the grid point locations of 
of overlapping boundaries are identified using the ``imask`` array. Thus, Dirichlet boundary condition values stored in 
``valint``, which contains spectral interpolation field values from overlapping mesh, are imposed at the boundary of 
current mesh, where ``imask = 1``.
 
..............................
SIZE file
..............................

The ``SIZE`` file used for this tutorial is included in the provided tar file. The user needs to ensure that the
auxiliary fields specified in the SIZE file is at minimum ``ldimt=3`` for RANS. Further, ``nsessmax`` must be
set to 2 for ``NEKNEK`` simulation. Other details on the contents of the ``SIZE`` file can be found 
:ref:`here<case_files_SIZE>`.

..............................
Compilation and Running
..............................

Compile from the parent directory with the usual command ``makenek``. 

A sample script for running the case on a cluster computing environment is included in the tar file (``neksaw``). 
The command in the script that launches the ``NEKNEK`` job is

.. code-block:: console

	neknek inlet bundle $((ntpn*nodes/2)) $((ntpn*nodes/2))
	
Here, ``nodes`` variable is the user input on number of nodes assigned for the job. ``ntpn`` is the number of processors/threads 
per node. First two parameters are the names of the component meshes and the following two parameters specify the number of total
threads used for each session, respectively. We use equal number of threads for this turorial, but the user may modify the
distribution of threads as needed. The script can be adopted suitably for any cluster being used. On Sawtooth cluster, 
the script is launched as follows:

.. code-block:: console 
	
	./neksaw inlet_bundle 40 4 30
	
The above runs ``NEKNEK`` job on 40 nodes (20 dedicated to each session) for 4 hours and 30 minutes. Remember to specify
the project name before launching, assigned to ``prj`` variable in the ``neksaw`` script file.

..............................
Helpful Tips
..............................

The following tips may be helpful to make the simulations more tractable:

 * Commence with a small time step size and high viscosity value (low Re) to stabilize the pressure solver
   during initial transients.
 * Accelerate the simulation by running standalone case for inlet component, allowing flow to evolve 
   before using ``NEKNEK`` solver for coupled simulation. Replace the ``int`` boundary condition with ``O`` (outlet)
   for the inlet component. The standalone case setup, if opted for, is left to the user as an exercise.
 * Use OIFS solver to run the simulation at larger time steps (CFL>1). This requires the following entries in the ``.par``
   file (uncomment to use):
 
.. literalinclude:: multi_rans/inlet_bundle/inlet.par
	:language: ini
	:lines: 11-12
	
It is necessary to specify target CFL for the OIFS solver. It calculates the number of extrapolation iterations
based on ``targetCFL`` value.

..............................
Results
..............................

For reference, normalized velocity magnitude and turbulent kinetic energy (TKE) contour plots are shown below along lateral
cross-sections of the inlet and wire-pin bundle components. Note that the overlap plane is located at :math:`z=0` axial location,
where the polar orientation of wire is :math:`7.2^{\circ}`. Contour plots along axial sections at successive downstream locations
are also shown, corresponding to wire polar orientiations of :math:`\{7.2,90,180,360\}^{\circ}`.

Dimensionless temperature axial distribution is shown in :numref:`fig:multi_rans_temp` along :math:`x`-section. Temperature contours
along axial cross-section near the bundle outlet is also included in the figure.
   
.. _fig:multi_rans_v:

.. figure:: multi_rans/results/umag_section.png
   :align: center
   :figclass: align-center

   Normalized velocity magnitude along :math:`x` and :math:`y` cross-sections. Corresponding zoomed views near overlap region (:math:`z=0`) shown.

.. _fig:multi_rans_k:

.. figure:: multi_rans/results/k_section.png
   :align: center
   :figclass: align-center

   Normalized TKE along :math:`x` and :math:`y` cross-sections. Corresponding zoomed views near overlap region (:math:`z=0`) shown.

.. _fig:multi_rans_vpolar:

.. figure:: multi_rans/results/umag_polar.png
   :align: center
   :figclass: align-center

   Normalized velocity magnitude in bundle component at successive downstream locations (corresponding polar orientation annotated).
  
.. _fig:multi_rans_kpolar:

.. figure:: multi_rans/results/k_polar.png
   :align: center
   :figclass: align-center

   Normalized TKE in bundle component at successive downstream locations (corresponding polar orientation annotated).

.. _fig:multi_rans_temp:

.. figure:: multi_rans/results/temp.png
   :align: center
   :figclass: align-center

   Dimensionless Temperature contour plot along :math:`x` cross-section. Axial cross-section near the outlet also shown.
  
