
----------------------------
Mesh conversion
----------------------------

In ``nek5000/tools``, there is a tool 'gmsh2nek' that convert Gmsh mesh (version 2, .msh file) into Nek mesh (.re2 file).

By running this tool, it will ask you input for dimensions (2 or 3), and then .msh file name (without extension).

=========================================================================================
set up boundary condition in usr file usrdat2() subroutine

      do iel=1,nelv
      do ifc=1,2*ndim
        id_face = bc(5,ifc,iel,1)
        if (id_face.eq.1) then        ! surface 1 for inlet 
           cbc(ifc,iel,1) = 'v  '
        elseif (id_face.eq.2) then    ! surface 2 for outlet
           cbc(ifc,iel,1) = 'O  '
        elseif (id_face.eq.3) then    ! surface 3 for wall
           cbc(ifc,iel,1) = 'W  '
        endif
      enddo
      enddo

or 

call setbc(1,1,'v  ') ! set bcID 1 to inlet for field 1 (velocity)
call setbc(2,1,'O  ') ! set bcID 2 to outlet for field 1 (velocity)
call setbc(3,1,'W  ') ! set bcID 3 to wall for field 1 (velocity)

=========================================================================================
About periodicity, currently only support translational periodicity.

Example input:

assuming there is only one pair of periodicity. you can have more than one pair. 
"
Enter number of periodic boundary surface pairs
1
"

surface numbers that need periodicity. surface 2 will map to surface 3
"
input surface 1 and  surface 2  sideSet ID
2 3
"

The suface 2 element faces will be mapped to surface 3 element faces accordingly.
However, this requires you have conformal surface 2 and 3 mesh.
 
'P  ' boundary tag will be assigned to cbc array.
bc(5,ifc,ie,1) still stores the surface number. 

==========================================================================================
About Conjugate Heat Transfer (CHT) mesh, you need to create a solid mesh that is conformal to the fluid mesh.

The following example shows how to assign boundary tag in usr file usrdat2() subroutine for CHT mesh

for velocity bc

      do iel=1,nelv
      do ifc=1,2*ndim
        id_face = bc(5,ifc,iel,1)
        if (id_face.eq.1) then        ! surface 1 for inlet 
           cbc(ifc,iel,1) = 'v  '
        elseif (id_face.eq.2) then    ! surface 2 for outlet
           cbc(ifc,iel,1) = 'O  '
        elseif (id_face.eq.3) then    ! surface 3 for wall
           cbc(ifc,iel,1) = 'W  '
        endif
      enddo
      enddo

for thermal bc

      do iel=1,nelt
      eg = gllel(iel)                 ! get global element number
      do ifc=1,2*ndim
        id_face = bc(5,ifc,iel,2)
        if (eg.le.nelgv) then           ! for fluid domain
          if (id_face.eq.1) then        ! surface 1 for inlet 
             cbc(ifc,iel,2) = 't  '
          elseif (id_face.eq.2) then    ! surface 2 for outlet
             cbc(ifc,iel,2) = 'O  '
          elseif (id_face.eq.3) then    ! surface 3 for wall, which connects to solid domain
             cbc(ifc,iel,2) = 'E  '
          endif
        else                            ! for solid domain
          if (id_face.eq.1) then        ! surface 1 for wall, which connects to fluid domain
             cbc(ifc,iel,2) = 'E  '
          elseif (id_face.eq.2) then    ! surface 2 for all external surfaces
             cbc(ifc,iel,2) = 'I  '
          endif
        
        endif
      enddo
      enddo








