
----------------------------
Mesh conversion
----------------------------

In ``nek5000/tools``, there is a tool 'exo2nek' that convert Exodus mesh (.exo file) into Nek mesh (.re2 file).

By running this tool, it will first ask you for number of fluid exo files, this allows you merge multiple exo files into one mesh. 

Then you input .exo file name (without extension) accordingly.

Then it will ask you for number of solid exo files.
if you do not have solid mesh, then just input 0.

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

surface numbers that need periodicity. sideset 2 will map to sideset 3
"
input surface 1 and  surface 2  sideSet ID
2 3
"

The sideset 2 element faces will be mapped to sideset 3 element faces accordingly.
However, this requires you have conformal sideset 2 and 3 mesh.
 
'P  ' boundary tag will be assigned to cbc array.
bc(5,ifc,ie,1) still stores the sideset number. 

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








