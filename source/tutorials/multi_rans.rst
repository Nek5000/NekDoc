.. _multi_rans:

-----------------------------------
Multi-Component RANS 
-----------------------------------

This is an advanced tutorial comprising RANS simulation in a multi-component geometry. The components
include an inlet nozzle coupled with a wire-pin bundle. Given the scale of the geometry and mesh used in this 
tutorial, it is assumed that the user has access to a computing cluster for running the case.
NekNek module is used for interfacing between the two components, allowing for runtime transfer of field data
and coupled simultaneous solve of flow and temperature equations. The tutorial employs several advanced 
Nek5000 tools, procedures and user routines including,

 * ``exo2nek`` for importing third party mesh file to Nek5000
 * ``reatore2`` for converting ASCII Nek5000 mesh to binary format
 * ``gencon`` for generating mesh connectivity file
 * ``wiremesher`` for generating wire-wrapped pin bundle
 * ``nekamg_setup`` for generating algebraic multi-grid solver (AMG) files
 * :math:`k-\tau` RANS simulation setup 
 *  ``NekNek`` setup

..........................
Before You Begin
..........................

It is highly recommended that new users familiarize themselves with the basic Nek5000 simulation
setup files and procedures outlined in the :ref:`fdlf` and :ref:`perhill` tutorials. Further, some 
sections of the tutorial will assume familiarity with :ref:`rans` and :ref:`neknek` tutorials.

..........................
Case Overview
..........................

All required files for this tutorial can be downloaded using the following link:

 * :download:`inlet_bundle.tar <multi_rans/inlet_bundle.tar>`

Extract the above tar file and navigate to the ``inlet_bundle`` (parent) directory

.. code-block:: console

	tar -xvf inlet_bundle.tar
	cd inlet_bundle
	
It has the following structure, discussed in detail in the following sections:

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

The geometry under consideration consists of two components including an inlet nozzle connected to a
wire-wrapped pin bundle, as shown: 

.. _fig:sfr_geom:

.. figure:: multi_rans/geom.png
   :align: center
   :figclass: align-center

   Component geometries including an inlet nozzle (red) connected to a wire-wrapped pin bundle
   
The geometry is non-dimensionalized with respect to the pin diameter, :math:`D`. 
Axial extent of the inlet nozzle component is :math:`z/D=[-4.75,0.25]` and the wire-pin bundle axial dimensions
are :math:`z/D=[0,12.5]`.  The two components, therefore, have an axial overlap of :math:`\Delta z/D= 0.25`,
while the lateral dimensions are conformal. Users should note that a reasonable overlap in the computational domains is imperative for a 
robust ``NekNek`` simulation. As a general rule, increasing the extent of overlap will render more stability
to the coupled solver. Pitch-to-diameter ratio of the wire-pin bundle is 1.13 and the wire diameter is 
:math:`D_w=0.12875D`. To prevent sharp corners, a fillet of diameter :math:`0.05 D` is introduced between the 
pin and the wire and the wire is slightly submerged in the pin by a distance of :math:`0.025D_w`.
The side length of the hexagonal bundle is :math:`1.865D` and the diameter of the inlet nozzle boundary 
is :math:`1.875D`. Length of the bundle component is equal to the wire pitch.

.....................................................
Thermal-hydraulic Parameters and Boundary Conditions
.....................................................

All flow and thermal properties are listed in :numref:`tab:nb_BC`. A non-dimensional bulk velocity of 
:math:`U_{in}=1` is specified at the inlet and inlet non-dimensional diameter is :math:`D_{in}=1.875`.
Hydraulic diameter of the pin-wire bundle is :math:`D_b\approx 0.398`. No-slip boundary condition is 
specified on all walls and the corresponding values of :math:`k` and :math:`\tau` on the walls is zero.

A constant non-dimensional temperature :math:`T^*=1` is specified at the inlet and a non-dimensional heat flux
of :math:`q'' =1` is applied on the pin walls from half axial length of the bundle. Insulated wall condition
is specified on all remaining boundaries.

.. _tab:nb_BC:

.. table:: Flow and Thermal Properties 

   +----------------------------+--------------------------------+
   |:math:`Re_{in}`             |60000                           |
   +----------------------------+--------------------------------+
   |:math:`Re_{b}`              |10000 (approx)                  |
   +----------------------------+--------------------------------+
   |:math:`\nu^*`               |3.125e-5                        |
   +----------------------------+--------------------------------+
   |:math:`Pr`                  |0.005                           |
   +----------------------------+--------------------------------+
   |:math:`Pe`                  |160                             |
   +----------------------------+--------------------------------+
   |:math:`Pr_t`                |1.5                             |
   +----------------------------+--------------------------------+

..............................
Mesh Generation
.............................. 

########################
Inlet Nozzle (exo2nek)
########################

A third-party meshing tool (e.g., ANSYS ICEM) is required for generating the mesh for the inlet component,
which must be saved as an ``EXODUS II (.exo)`` mesh file. Nek5000 offers the ``exo2nek`` mesh conversion tool
for converting an ``.exo`` mesh file to ``.re2`` format. Ensure that the ``exo2nek`` tool is compiled, available
in the  ``Nek5000`` directory:

.. code-block:: console

	cd ~/Nek5000/tools
	./maketools all

The above commands will compile all available Nek5000 tools including ``exo2nek``. 

:Note:
	Ensure that ``MAXNEL`` parameter in ``maketools.inc`` file is set to a high value. Default value is 150000. 
	Set to a number greater than the number of elements in the mesh (1000000) before running the above commands.

Currently, ``exo2nek`` supports the following mesh elements,

 * ``TET4``
 * ``WEDGE6``
 * ``HEX8``
 * ``HEX20``

User must ensure that the third-party mesh comprises only the above listed element types. Navigate to the folder 
containing inlet mesh and run ``exo2nek``. Provide inputs to the prompts as shown:

.. code-block:: console
	
	cd inletMesh
	exo2nek
	
	please input number of fluid exo files: 
	1

.. code-block:: console

	please input exo file: 
	inlet
	
.. code-block:: console

	inlet.exo                        is an EXODUSII file; version 0.00
	I/O word size 8

	database parameters:

	title         =  Created by ICEMCFD - EXODUS II Interface                            

	num_dim       =        3
	num_nodes     =    36527
	num_elem      =   167965
	num_elem_blk  =        1
	num_side_sets =        3

	element block id   =        1
	element type       =    TETRA
	num_elem_in_block  =   167965
	num_nodes_per_elem =        4

	TETRA4 is valid element in a 3D mesh.
	assume linear hybrid mesh (tetra-hex-wedge)
	one TETRA4 divide into 4 Nek hex elements
	please input number of solid exo files for CHT problem (input 0 for no solid mesh): 
	0
	
.. code-block:: console

	done pre-read exo files
	now converting to nek mesh


	Store SideSet information from EXO file
	Sideset  2 ...
	Sideset  3 ...
	Sideset  4 ...

	Converting elements ...
	 flag1
	 flag2
	 Converting elements in block            1
	 nvert,                     4
	 Converted elements in nek:               671860
	Done :: Converting elements
	 Domain max xyz:   1.8650400000000000        1.6151700000000000        0.0000000000000000
	 Domain min xyz:  -1.8650400000000000       -1.6151700000000002       -5.0000000000000000
	 total element now is                671860
	 fluid exo file            1  has elements       671860
	 calling: gather_bc_info()
	 done: gather_bc_info()
	 ******************************************************
	 Boundary info summary
	 sideSet ID
           2
           3
           4
	 ******************************************************
	Enter number of periodic boundary surface pairs: 
	0
	

.. code-block:: console

	please give re2 file name: 
	inlet
	
Following the above steps will generate the file ``inlet.re2`` in the current directory. Note that the 
``sideSet ID`` for all mesh boundaries must be specified in the ``.exo`` file using the third-party 
meshing software of user's choice. For the inlet nozzle component, the following IDs are assigned:

 * 2 --> Nozzle inlet
 * 3 --> Interfacing surface (bundle surface)
 * 4 --> Walls

Return and move the mesh file to the parent directory:

.. code-block:: console

	cd ../
	mv inletMesh/inlet.re2 .

#############################
Wire-pin Bundle (wiremesher)
#############################

To generate the pin-wire bundle mesh, navigate to the ``wireMesh`` folder:

.. code-block:: console

	cd wireMesh
	
It contains two sub-directories, viz., ``wire2nek`` and ``matlab``. Input parameters for the meshing 
script are specified in the header of the ``matlab/wire_mesher.m`` file, as shown:

.. code-block:: console

	% Mesh parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	D    = 8.00;           % pin diameter  (mm)
	P    = 9.04;           % pin pitch
	Dw   = 1.03;           % wire diameter (mm)
	Df   = 0.4;            % fillet diameter (mm) (making this too small can cause bad elements)
	H    = 100.0;            % wire pitch (mm)
	T    = 0.05*Dw         % trimmed off of wire (mm) (cuts off the tip of the wire)
	S    = 0.025*Dw;       % Wire submerged (mm) (how sunken in the wire is into the pin)
	Adjust = 1;            % if 1, Adjust flattness of wire when away from pins (if trimmed is on, only trim when wire passes pin)
	iFtF   = 0;            % if 1, add layer next to outer can
	G  = 0.0525;           % gap between wire and wall (mm)
	FtF_rescale = 1.0086;  % rescaling of outer FtF, the difference is given by an additional mesh layer - MUST BE BIGGER THAN 1
	ne=2;                  % Number of edge pins. e.g., ne=3 for 19 pin assembly. MINIMUM is 2
	Col=12;                % Number of columns per block (5 or 6 blocks surround each pin)
	Row=3;                 % Number of rows per block (2 blocks fit between neighboring pins)
	Rowdist=[65 30 5];     % distribution of elements in the row (for creating boundary layer)
	Lay=20;                % Number of layers for 60 degree turn (this should be a reasonable number - above 10, tested typically for ~20)
	rperiodic=1.0;         % Less than zero if periodic, Greater than zero if inlet/outlet
	ipolar=0.25*(D/H)*360. % Starting polar orientation of wire (in degrees)

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
 * ``Lay`` --> Controls resolution in axial direction. Specifies number of axial elements in length equal to 60 degree rotation of wire.  

The mesher is initiated by simply running the ``doall_binary`` bash file. Ensure that both Matlab and Python 
(tested with Python 3.8) are active before launching the script and that Fortran compilers are available. 

.. code-block:: console

	./doall_binary
	
The script can take a while to complete. Upon completion it generates ``wire_out.rea`` mesh file, which is the 
ASCII mesh file for Nek5000. Convert this into the binary format by running ``reatore2`` tool. Follow the prompts:

.. code-block:: console

	Input .rea name:
	wire_out
	
.. code-block:: console

	Input .rea/.re2 output name:
	bundle
	
We finally obtain the ``bundle.re2`` file which contains the pin-wire bundle mesh for Nek5000 run. Boundary IDs
are assigned by the wire mesher as:

 * 1 --> Fuel pin walls
 * 2 --> Bundle hexagonal (outer) walls
 * 3 & 4 --> Axial end surfaces

Return and move the mesh file to the parent directory:

.. code-block:: console

	cd ../
	mv bundleMesh/bundle.re2 .
 
####################################
Generating Connectivity file (.co2)
####################################

After generating the mesh files for both components, it is necessary to generate the corresponding connectivity files
using ``gencon`` tool. Note that using ``gencon`` tool instead of ``genmap``, which generates  map (``.ma2``) file, is the 
recommended procedure for large meshes. See :ref:`build_pplist` for details. 

Users must specify ``PPLIST=PARRSB`` in ``makenek`` file (location: Nek5000/bin/makenek) which instructs Nek5000
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

.. code-block:: console

	#
	# nek parameter file - inlet.par
	#
	[GENERAL]
	#startFrom = inlet.fld
	stopAt = numSteps
	numSteps = 20000
	dt = 1.0e-5
	writeInterval = 5000
	timeStepper = BDF2
	#targetCFL=4.0
	#extrapolation = OIFS

	[PROBLEMTYPE]
	variableProperties = yes
	stressFormulation = yes

	[PRESSURE]
	preconditioner = semg_amg
	residualTol = 1.0e-5
	residualProj = yes

	[VELOCITY]
	density = 1.0
	viscosity = -32000.0
	residualTol = 1.0e-6

	[TEMPERATURE]
	solver = none
	rhoCp = 1.0
	conductivity = -160.0
	residualTol = 1.0e-6

	[SCALAR01] 
	density = 1.0
	diffusivity = -32000.0
	residualTol = 1.0e-6

	[SCALAR02] 
	density = 1.0
	diffusivity = -32000.0
	residualTol = 1.0e-6

.. code-block:: console

	#
	# nek parameter file - bundle.par
	#
	[GENERAL]
	#startFrom = bundle.fld
	stopAt = numSteps
	numSteps = 20000
	dt = 1.0e-5
	writeInterval = 5000
	timeStepper = BDF2
	#targetCFL = 4.0
	#extrapolation = OIFS

	[PROBLEMTYPE]
	variableProperties = yes
	stressFormulation = yes

	[PRESSURE]
	preconditioner = semg_amg
	residualTol = 1.0e-5
	residualProj = yes

	[VELOCITY]
	density = 1.0
	viscosity = -32000.0
	residualTol = 1.0e-6

	[TEMPERATURE]
	solver = none
	rhoCp = 1.0
	conductivity = -160.0
	residualTol = 1.0e-6

	[SCALAR01] 
	density = 1.0
	diffusivity = -32000.0
	residualTol = 1.0e-6

	[SCALAR02] 
	density = 1.0
	diffusivity = -32000.0
	residualTol = 1.0e-6
	
Both parameter files are identical, except for one important difference. If the user wants to restart the case from
any given time, separate restart file names should be specified to the ``startFrom`` parameter. It is critical that the
properties and time step size are identical for both ``.par`` files.

Note that given the large size of meshes, the ``preconditioner`` must be set to ``semg_amg``. This invokes the algebraic
multigrid solver for pressure instead of the default ``XXT`` solver.

Further details on all parameters of ``.par`` file can be found :ref:`here <case_files_par>`.
 
.........................................
User Routines (.usr file)
......................................... 

Basics of the required setup routines for a NekNek simulation can be found in the :ref:`neknek` turorial, while for a RANS simulation
in the :ref:`rans` tutorial. Although this section decribes all user routines required for a NekNek RANS simulation in detail, 
a comprehensive understanding of routines from these simpler cases is recommended before proceeding.

Following headers are required at the beginning of ``.usr`` file for loading RANS related subroutines:

.. code-block:: console

	include "experimental/rans_komg.f"
	include "experimental/rans_wallfunctions.f"

``NekNek`` related parameters are specified in ``usrdat`` routine:

.. code-block:: console

	subroutine usrdat() 
	include 'SIZE'
	include 'TOTAL'

	!   ngeom - parameter controlling the number of iterations,
	!   set to ngeom=2 by default (no iterations) 
	!   One could change the number of iterations as
	ngeom = 2

	!   ninter - parameter controlling the order of interface extrapolation for neknek,
	!   set to ninter=1 by default
	!   Caution: if ninter greater than 1 is chosen, ngeom greater than 2 
	!   should be used for stability
	ninter = 1
	
	nfld_neknek=7 !velocity+pressure+t+sc1+sc2 

	return
	end
	
``ngeom`` specifies the number of overlapping Schwarz-like iterations, while ``ninter`` controls the time 
extrapolation order of boundary conditions of the overlapping interface. ``ninter=1`` is unconditionally 
stable, while a higher temporal order will typically require more iterations for stability (``ngeom>2``). 
For computational savings, we maintain first order temporal extrapolation for this tutorial. 
``nfld_neknek`` specifies the number of total field arrays that are transferred between the two meshes 
and must be equal to 7 for 3D RANS cases (3 velocity, 1 pressure and 3 scalar field arrays - temperature,
:math:`k` and :math:`\tau`).

:Note:
	Ensure that proper common block headers are included in subroutines. ``NEKNEK`` header is required 
	for routines where ``idsess`` needs to be accessed, as shown below.
	
Boundary Condition specification and RANS initialization is performed in ``usrdat2``:

.. code-block:: console

	subroutine usrdat2()
	implicit none
	include 'SIZE'
	include 'TOTAL'
	include 'NEKNEK'

	real wd
	common /walldist/ wd(lx1,ly1,lz1,lelv)

	common /rans_usr/ ifld_k, ifld_omega, m_id
	integer ifld_k,ifld_omega, m_id

	integer w_id,imid,i
	real coeffs(30) !array for passing your own coeffs
	logical ifcoeffs
    
	integer ifc,iel,id_face
	
	if(idsess.eq.0)then              !BCs for inlet mesh
	  do iel=1,nelt
	  do ifc=1,2*ndim
	     id_face = bc(5,ifc,iel,1)
	     if (id_face.eq.2) then      !inlet
		 cbc(ifc,iel,1) = 'v  '
		 cbc(ifc,iel,2) = 't  '
		 cbc(ifc,iel,3) = 't  '
		 cbc(ifc,iel,4) = 't  '
	     elseif (id_face.eq.3) then  !interface
		 cbc(ifc,iel,1) = 'int'
		 cbc(ifc,iel,2) = 'int'
		 cbc(ifc,iel,3) = 'int'
		 cbc(ifc,iel,4) = 'int'
	     elseif (id_face.eq.4) then  !walls
		 cbc(ifc,iel,1) = 'W  '
		 cbc(ifc,iel,2) = 'I  '
		 cbc(ifc,iel,3) = 't  '
		 cbc(ifc,iel,4) = 't  '
	     endif
	  enddo
	  enddo
	else                             !BCs for pin-wire bundle mesh
	  do iel=1,nelt
	  do ifc=1,2*ndim
	     id_face = bc(5,ifc,iel,1)
	     if (id_face.eq.3) then      !interface 
		 cbc(ifc,iel,1) = 'int'
		 cbc(ifc,iel,2) = 'int'
		 cbc(ifc,iel,3) = 'int'
		 cbc(ifc,iel,4) = 'int'
	     elseif (id_face.eq.4) then  !outlet 
		 cbc(ifc,iel,1) = 'O  '
		 cbc(ifc,iel,2) = 'O  '
		 cbc(ifc,iel,3) = 'O  '
		 cbc(ifc,iel,4) = 'O  '
	     elseif (id_face.eq.1) then  !pin walls
		 cbc(ifc,iel,1) = 'W  '
		 cbc(ifc,iel,2) = 'f  '
		 cbc(ifc,iel,3) = 't  '
		 cbc(ifc,iel,4) = 't  '
	     elseif (id_face.eq.2) then  !outer walls
		 cbc(ifc,iel,1) = 'W  '
		 cbc(ifc,iel,2) = 'I  '
		 cbc(ifc,iel,3) = 't  '
		 cbc(ifc,iel,4) = 't  '
	     endif
	  enddo
	  enddo
	endif
	
	! RANS initialization
	ifld_k     = 3 
	ifld_omega = 4
	ifcoeffs=.false.

	m_id = 4 ! non-regularized standard k-tau
	w_id = 2 ! distf (coordinate difference, provides smoother function)

	call rans_init(ifld_k,ifld_omega,ifcoeffs,coeffs,w_id,wd,m_id)

	return
	end

``NekNek`` solver launches two Nek5000 sessions simultaneously and field data transfer is performed between the
two sessions after each time iteration. Each session is assigned a unique id, stored in the variable ``idsess``.
Here, ``idsess=0`` is assigned to the inlet component solve and ``idsess=1`` to bundle component. Boundary 
conditions are assigned using this variable for each component, as shown above. 

Recall the boundary IDs assigned to each component, described in the preceding section. Character codes
for different boundary conditions are stored in the ``cbc`` array. Their detailed description can be found in
:ref:`boundary-conditions`. For each component, the nested loops go through all elements and their faces and
populates ``cbc`` array for all fields based on mesh assigned boundary IDs. Note that ``int`` boundary condition
must be assigned to the overlapping surfaces on the inlet and bundle components. ``int`` condition is replaced 
internally with Dirichlet boundary conditions subsequently by Nek5000. Insulated 

With regards to RANS initialization; ``m_id=4`` selects the :math:`k-\tau` RANS model and ``w_id=2`` selects the
wall distance computing algorithm. :math:`k` and :math:`\tau` fields are stored in the 3rd and 4th index, respectively, 
specified with ``ifld_k`` and ``ifld_omega``. Set ``ifcoeffs`` to ``.true.`` only if user specified RANS coefficients
are required. For details on the RANS related parameters, refer :ref:`rans` tutorial.  

:Note:
	``rans_init`` must be called after populating ``cbc`` array

For RANS simulation, diffusion coefficients are assigned in the ``uservp`` routine. The routine used here remains 
nearly identical to the :ref:`rans` tutorial:

.. code-block:: console

	subroutine uservp (ix,iy,iz,eg)
	implicit none
	include 'SIZE'
	include 'TOTAL'
	include 'NEKUSE'

	integer ix,iy,iz,e,eg
	
	common /rans_usr/ ifld_k, ifld_omega, m_id
	integer ifld_k,ifld_omega, m_id
	
	real rans_mut,rans_mutsk,rans_mutso,rans_turbPrandtl
	real mu_t,Pr_t

	e = gllel(eg)

	Pr_t=1.5 !rans_turbPrandtl()
	mu_t=rans_mut(ix,iy,iz,e)

	utrans = cpfld(ifield,2)
	if(ifield.eq.1) then
		udiff = cpfld(ifield,1)+mu_t
	elseif(ifield.eq.2) then
		udiff = cpfld(ifield,1)+mu_t*cpfld(ifield,2)/(Pr_t*cpfld(1,2))
	elseif(ifield.eq.ifld_k) then  
		udiff = cpfld(1,1)+rans_mutsk(ix,iy,iz,e)
	elseif(ifield.eq.ifld_omega) then  
		udiff = cpfld(1,1)+rans_mutso(ix,iy,iz,e)
	endif

	return
	end

Only turbulent Prandtl number is changed to ``Pr_t=1.5`` for this tutorial. This value is more appropriate for 
molten sodium salts as compared to the default value of 0.85 (for air), which is assigned through ``rans_turbPrandtl()`` 
function call.

Source terms for the temperature and scalar equations are assigned through ``userq``. The routine here is identical
to the basic :ref:`rans` case:

.. code-block:: console

	subroutine userq  (ix,iy,iz,ieg)
	implicit none
	include 'SIZE'
	include 'TOTAL'
	include 'NEKUSE'

	common /rans_usr/ ifld_k, ifld_omega, m_id
	integer ifld_k,ifld_omega, m_id

	real rans_kSrc,rans_omgSrc
	real rans_kDiag,rans_omgDiag

	integer ie,ix,iy,iz,ieg
	ie = gllel(ieg)

	if(ifield.eq.2) then
		qvol = 0.0 
		avol = 0.0
	elseif(ifield.eq.ifld_k) then
		qvol = rans_kSrc  (ix,iy,iz,ie)
		avol = rans_kDiag (ix,iy,iz,ie)
	elseif(ifield.eq.ifld_omega) then
		qvol = rans_omgSrc (ix,iy,iz,ie)
		avol = rans_omgDiag(ix,iy,iz,ie)
	endif
	
	return
	end
	
Note that either component does not have any volumetric source heat source and hence ``qvol=0`` for ``ifield .eq. 2``.

Initial conditions are specified in ``useric``. Similar values are assigned to both components and thus, the routine
implementation is straightforward. Temperature is initalized to 1 for both components.

.. code-block:: console

	subroutine useric (ix,iy,iz,eg)
	implicit none
	include 'SIZE'
	include 'TOTAL'
	include 'NEKUSE'

	integer ix,iy,iz,e,eg
	
	common /rans_usr/ ifld_k, ifld_omega, m_id
	integer ifld_k,ifld_omega, m_id
	
	e = gllel(eg)

	ux   = 0.0
	uy   = 0.0
	uz   = 1.0
	temp = 1.0

	if(ifield.eq.2) temp = 1.0
	if(ifield.eq.ifld_k) temp = 0.01
	if(ifield.eq.ifld_omega) temp = 0.2

	return
	end

Boundary conditions are assigned in ``userbc``. For the inlet component, inlet conditions are assigned using data
generated from RANS simulation in a pipe with idential diameter as the inlet surface. The data is stored in the 
``InletProf.dat`` file which contains axial velocity, :math:`k` and :math:`\tau` information as a function of
radial wall distance. Two plugin subroutines are required, which perform spline interpolation of the data to the 
inlet surface, viz., ``getInletProf`` and ``init_prof``. These are provided for the user in the ``inlet_bundle.usr``
file and can be used without modification. The usage is shown below:
 
.. code-block:: console

	subroutine userbc (ix,iy,iz,iside,eg)
	implicit none
	include 'SIZE'
	include 'TOTAL'
	include 'NEKUSE'
	include 'NEKNEK'
	
	integer ix,iy,iz,iside,eg,e
	
	real wd
	common /walldist/ wd(lx1,ly1,lz1,lelv)
	
	common /rans_usr/ ifld_k, ifld_omega, m_id
	integer ifld_k,ifld_omega, m_id
	
	real uin,kin,tauin,wdist
	real din,zmid
	
	integer id_face
	
	integer icalld
	save icalld
	data icalld /0/
	
	e = gllel(eg)
	
	if(icalld.eq.0)then
	  call getInletProf                      !Initialize spline routine
	  icalld = 1
	endif

	ux = 0.0
	uy = 0.0
	
	id_face = bc(5,iside,e,1)
	
	! Inlet condition
	if(idsess.eq.0 .and. id_face.eq.2)then
	  din = 1.875                            !Inlet diameter
	  wdist = min(wd(ix,iy,iz,e),din/2.)     !Wall distance
	  call init_prof(wdist,uin,kin,tauin)    !Spline interpolation from InletProf.dat
	  uz = uin
	  if(ifield.eq.2)temp = 1.0
	  if(ifield.eq.ifld_k)temp = kin
	  if(ifield.eq.ifld_omega)temp = tauin
	endif
	
	! Heat flux on walls
	if(ifield.eq.2)then
	  if(idsess.eq.0)then
	    flux = 0.0
	  else
	    zmid = 12.5/2.0
	    flux=0.0
	      if(id_face.eq.1)then               !Pin walls
	        flux = 0.5*(1.0+tanh(2.0*PI*(zm1(ix,iy,iz,e)-zmid)))
	      endif
	  endif
	endif
	
	return
	end
	
Note that the diameter of the inlet surface is ``din=1.875``. Spline interpolation routine, ``init_prof``, requires
the wall distance array, ``wd``, which is populated in ``usrdat2`` (in ``rans_init`` call). The distance should be
limited to inlet radius to avoid spline from extrapolating. Inlet component is identified with ``idsess.eq.0`` and
inlet surface with its boundary ID, ``id_face.eq.2``.

Temperature flux must also be assigned in ``userbc`` on the pin surface walls. As mentioned earlier, non-dimensional
unit heat flux is assigned from half axial length of the bundle (``zmid``). A smoothed flux profile is imposed using
``tanh`` step function as shown. Flux on all remaining walls is zero. 
 
..............................
SIZE file
..............................

The ``SIZE`` file used for this tutorial is inlcuded in the provided tar file. The user needs to ensure that the
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
	
Here, ``nodes`` variable is the user input on number of nodes. ``ntpn`` is the number of processors per node. First two
parameters are the names of the component meshes and the following two parameters specify the number of total threads
used for each session, respectively. We use equal number of threads for this turorial, but the user may modify the
distribution of threads as needed. The script can be adopted suitably for the cluster being used. On Sawtooth cluster, 
the script is launched as follows:

.. code-block:: console 
	
	./neksaw inlet_bundle 40 4 30
	
The above runs ``NEKNEK`` job on 40 nodes (20 dedicated to each session) for 4 hours and 30 minutes. Remember to specify
the project name before launching, assigned to ``prj`` variable in the script.

On the first run, Nek5000 will generate files for setting up the AMG (algebraic multi-grid) solver with the suffix
``amgdmp_p.dat``, ``amgdmp_i.dat`` and ``amgdmp_j.dat`` for each component. Run the ``nekamg_setup`` tool to 
create the setup files. The prompts will appear as shown:

.. code-block:: console

	Enter name prefix of input file(s):
	inlet
	
.. code-block:: console

	Choose a coarsening method. Available options are:
	- 3: Ruege-Stuben,
	- 6: Falgout (default),
	- 8: PMIS,
	- 10: HMIS,
	Choice:
	10

.. code-block:: console

	Choose an interpolation method. Available options are:
	- 0: classical modified interpolation,
	- 6: extended + i interpolation (default),
	Choice:
	0

User may choose any of the coarsening solver and interpolation methods available which succeffully converges.
For this tutorial we choose HMIS and classical interpolation for both components. 

.. code-block:: console
	
	Enter smoother tolerance [0.5]:
	0.9
	
Repeat the steps for the bundle component. Upon completion three files will be written containing AMG matrices.

On the next run, Nek5000 will run normally and the user may proceed with the simulation.

..............................
Helpful Tips
..............................

The following tips may be helpful to make the simulations tractable:

 * Commence with a small time step size and high viscosity value (low Re) to stabilize the pressure solver
   during initial transients.
 * Accelerate the simulation by running standalone case for inlet component, allowing flow to evolve 
   before using ``NEKNEK`` solver for coupled simulation. Replace the ``int`` boundary condition with ``O`` (outlet)
   for the inlet component. The standalone case setup, if opted for, is left to the user as an exercise.
 * Use OIFS solver to run the simulation at larger time steps (CFL>1). This requires the following entries in the ``.par``
   file:
 
.. code-block:: console
	
	targetCFL=4.0
	extrapolation = OIFS
	
It is necessary to specify target CFL for the OIFS solver to prescribe the number of extrapolation iterations.


..............................
Results
..............................
   

