.. _multi_rans:

-----------------------------------
Multi-Component :math:`k-\tau` RANS 
-----------------------------------

This is an advanced tutorial comprising RANS simulation in a multi-component geometry. The components
include an inlet nozzle coupled with a wire-pin bundle. Given the scale of the geometry and mesh used in this 
tutorial, it is assumed that the user has access to cluster computing system for running the case.
NekNek module is used for interfacing between the two components, allowing for runtime transfer of field data
and coupled simultaneous solve of flow and temperature equations. The tutorial will describe several advanced 
Nek5000 tools, procedures and user routines including,

 * ``exo2nek`` for importing mesh file to Nek5000
 * ``wiremesh`` for generating wire-wrapped pin bundle
 * :math:`k-\tau` RANS simulation setup 
 *  ``NekNek`` setup

..........................
Before You Begin
..........................

It is highly recommended that new users familiarize themselves with the basic Nek5000 simulation
setup files and procedures outlined in the :ref:`fdlf` and :ref:`perhill` tutorials. Further, some 
sections of the tutorial will assume familiarity with :ref:`rans` and :ref:`neknek` tutorials.

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
   
The geometry is non-dimensionalized with respect to the pin diameter, :math:`D`. Accordingly, the 
axial extent of the inlet nozzle component is :math:`z/D=[-5,0]` and the wire-pin bundle axial dimensions
are :math:`z/D=[-0.25,12.250]`. The pitch of the wire wrap is :math:`12.5D`. A portion of the
bundle, :math:`\Delta z/D=0.5`, is included in the inlet component to allow for overlap. The two 
components, therefore, have an axial overlap of :math:`\Delta z/D= 0.25`, while the lateral dimensions are
conformal. Users should note that a reasonable overlap in the computational domains is imperative for a 
robust ``NekNek`` simulation. As a general rule, increasing the extent of overlap will render more stability
to the coupled solver. Pitch-to-diameter ratio of the wire-pin bundle is 1.13 and the wire diameter is 
:math:`0.12875D`. To prevent sharp corners, a fillet of diameter :math:`0.05 D` is introduced between the 
pin and the wire. The side length of the hexagonal bundle is :math:`1.865D` and the diameter of the inlet 
nozzle boundary is :math:`1.875D`. 

..............................
Mesh Generation
.............................. 

A third-party meshing tool (e.g., ANSYS ICEM) is required for generating the mesh for the inlet component,
which must be saved as an ``EXODUS II (.exo)`` mesh file. Nek5000 offers the ``exo2nek`` mesh conversion tool
for converting an ``.exo`` mesh file to ``.re2`` format. Ensure that the ``exo2nek`` tool is compiled, available
in the  ``Nek5000`` directory:

.. code-block:: console

	cd ~/Nek5000/tools
	./maketools exo2nek

Currently, ``exo2nek`` supports the following mesh elements,

 * ``TET4``
 * ``WEDGE6``
 * ``HEX8``
 * ``HEX20``

User must ensure that the third-party mesh comprises only the above listed element types. The ``EXODUS II`` file 
for the inlet component can be downloaded from the following link:

 * :download:`inlet.exo <multi_rans/inlet/inlet.exo>`
 
Run ``exo2nek`` tool from the directory containing the above file and provide inputs to the prompts as shown:

.. code-block:: console

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
	
Following the above steps will generate the file ``fluid.re2`` in the current directory. 