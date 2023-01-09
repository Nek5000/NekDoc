
.. _tools_exo2nek:

-------
exo2nek
-------

*Nek5000* provides the Exodus mesh conversion tool, ``exo2nek``.
This tool converts Exodus meshes (``.exo`` files) into Nek meshes (``.re2`` files).
The source code is located in ``nek5000/tools``, and it can be compiled from that directory with ``./maketools exo2nek``.
In addition to the compilers necessary to use *Nek5000*, ``exo2nek`` requires ``cmake``.
See :ref:`Before You Begin <qstart_before>` for more information on dependencies.
Once it is compilied, the executable will be available in ``nek5000/bin``.

When running this tool, the input steps are as follows:

* It will first ask you for the number of fluid ``.exo`` files, this allows you merge multiple ``.exo`` files into one mesh. 
* Then you input ``.exo`` file name (without extension) accordingly.
* It will then ask for the number of solid ``.exo`` files for a conjugate heat transfer mesh.
  If you do not have solid mesh, then just input 0.
* The tool then starts converting the mesh and prints out some basic information, such as the min and max dimensions, number of elements, and the sideSet IDs.
* Next, it will ask for the number of periodic boundary pairs (see :ref:`exo_pbound` below). 
  If you have none, simply enter 0.
* Finally it will ask for the name of the output file. The ``.re2`` extension is implied and should not be included.

An example is shown below with the expected user input highlighted.

.. literalinclude:: ../tutorials/multi_rans/exo2nek.output
   :language: none
   :emphasize-lines: 2,4,28,59,61

When using a converted Exodus mesh, the user must set up boundary condition in the ``.usr`` file in either the ``usrdat`` or ``usrdat2`` subroutines, e.g.

.. code-block:: fortran

      do iel=1,nelv
      do ifc=1,2*ndim
        id_face = BoundaryID(ifc,iel)
        if (id_face.eq.1) then        ! surface 1 for inlet 
           cbc(ifc,iel,1) = 'v  '
        else if (id_face.eq.2) then    ! surface 2 for outlet
           cbc(ifc,iel,1) = 'O  '
        else if (id_face.eq.3) then    ! surface 3 for wall
           cbc(ifc,iel,1) = 'W  '
        endif
      enddo
      enddo

or 

.. code-block:: fortran

  call setbc(1,1,'v  ') ! set bcID 1 to inlet for field 1 (velocity)
  call setbc(2,1,'O  ') ! set bcID 2 to outlet for field 1 (velocity)
  call setbc(3,1,'W  ') ! set bcID 3 to wall for field 1 (velocity)

.. _exo_pbound:

...................
Periodic Boundaries
...................

Periodic boundaries in *Nek5000* are implemented on a mesh connectivity level.
They must have a conformal element face distribution with a consistent offset vector.
Only translational periodicity is supported in ``exo2nek``.

:Example:
  The example below describes the expected inputs for a mesh with a single pair of periodic boundaries.

First, the user must provide the number of periodic surface pairs. 
In this case, we have 1.

.. code-block:: console

  Enter number of periodic boundary surface pairs
  1

Next the user specifies which surface IDs are periodic.
Here, we set the inlet (sideset 1) to be periodic with the outlet (sideset 2).

.. code-block:: console 

  input surface 1 and  surface 2  sideSet ID
  1 2

The sideset 1 element faces will be mapped to sideset 2 element faces accordingly.
However, this requires that you have conformal meshes on sidesets 1 and 2.
The ``P`` boundary tag will be assigned to the ``cbc`` array, while the sideset ID number is still avaialble in the ``BoundaryID`` array. 

.. Note::

  Mulitple pairs of periodic boundaries are supported

...................
Automatic Tet-2-Hex
...................

As *Nek5000* supports only hexahedral elements, ``exo2nek`` includes a feature that automatically converts tetrahedral and prism meshes to pure hexahedral meshes.
All tetrahedral elements are converted to 4 hexahedral elements and all wedge elements are converted to 3 hexahedral elements. 
These conversions are supported for both 1\ :superscript:`st` and 2\ :superscript:`nd` order elements.

* ``TET4`` + ``WEDGE6`` --> ``HEX8``
* ``TET10`` + ``WEDGE15`` --> ``HEX20``

..............................
Conjugate Heat Transfer Meshes
..............................

Coming soon

.. About Conjugate Heat Transfer (CHT) mesh, you need to create a solid mesh that is conformal to the fluid mesh.
.. 
.. The following example shows how to assign boundary tag in usr file usrdat2() subroutine for CHT mesh
.. 
.. for velocity bc
.. 
..       do iel=1,nelv
..       do ifc=1,2*ndim
..         id_face = bc(5,ifc,iel,1)
..         if (id_face.eq.1) then        ! surface 1 for inlet 
..            cbc(ifc,iel,1) = 'v  '
..         elseif (id_face.eq.2) then    ! surface 2 for outlet
..            cbc(ifc,iel,1) = 'O  '
..         elseif (id_face.eq.3) then    ! surface 3 for wall
..            cbc(ifc,iel,1) = 'W  '
..         endif
..       enddo
..       enddo
.. 
.. for thermal bc
.. 
..       do iel=1,nelt
..       eg = gllel(iel)                 ! get global element number
..       do ifc=1,2*ndim
..         id_face = bc(5,ifc,iel,2)
..         if (eg.le.nelgv) then           ! for fluid domain
..           if (id_face.eq.1) then        ! surface 1 for inlet 
..              cbc(ifc,iel,2) = 't  '
..           elseif (id_face.eq.2) then    ! surface 2 for outlet
..              cbc(ifc,iel,2) = 'O  '
..           elseif (id_face.eq.3) then    ! surface 3 for wall, which connects to solid domain
..              cbc(ifc,iel,2) = 'E  '
..           endif
..         else                            ! for solid domain
..           if (id_face.eq.1) then        ! surface 1 for wall, which connects to fluid domain
..              cbc(ifc,iel,2) = 'E  '
..           elseif (id_face.eq.2) then    ! surface 2 for all external surfaces
..              cbc(ifc,iel,2) = 'I  '
..           endif
..         
..         endif
..       enddo
..       enddo

