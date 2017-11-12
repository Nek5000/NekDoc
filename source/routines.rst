====================
Routines of Interest
====================

The most common routines needed in Nek5000 are

------------------
Naming Conventions
------------------

- ``subroutine f(a,b,c)``

  - ``a``- returned variable
  - ``b,c`` -input data

- ``op[]`` represent operations on operators
- ``c[]``  operations on constants
- ``gl[]`` global operations
- ``col2[] ---``
- ``col3[] --``


-----------
Subroutines
-----------

``subroutine rescale_x(x,x0,x1)``
    Rescales the array ``x`` to be in the range ``(x0,x1)``. This is usually called from ``usrdat2`` in the ``.usr`` file

``subroutine normvc(h1,semi,l2,linf,x1,x2,x3)``
    Computes the error norms of a vector field variable ``(x1,x2,x3)`` defined on mesh 1, the velocity mesh. The error norms are normalized with respect to the volume, with the exception on the infinity norm, ``linf``.

``subroutine comp_vort3(vort,work1,work2,u,v,w)``
    Computes the vorticity (``vort``) of the velocity field, ``(u,v,w)``

``subroutine lambda2(l2)``
    Generates the Lambda-2 vortex criterion proposed by Jeong and Hussain (1995)

``subroutine planar_average_z(ua,u,w1,w2)``
    Computes the r-s planar average of the quantity ``u``.

``subroutine torque_calc(scale,x0,ifdout,iftout)``
    Computes torque about the point ``x0``. Here scale is a user supplied multiplier so that the results may be scaled to any convenient non-dimensionalization. Both the drag and the torque can be printed to the screen by switching the appropriate ``ifdout(drag)`` or ``iftout(torque)`` logical.

``subroutine set_obj``
    Defines objects for surface integrals by changing the value of ``hcode`` for future calculations. Typically called once within ``userchk`` (for ``istep = 0``) and used for calculating torque. (see above)

``subroutine avg1(avg,f, alpha,beta,n,name,ifverbose)``

``subroutine avg2(avg,f, alpha,beta,n,name,ifverbose)``

``subroutine avg3(avg,f,g, alpha,beta,n,name,ifverbose)``
    These three subroutines calculate the (weighted) average of ``f``. Depending on the value of the logical, ``ifverbose``, the results will be printed to standard output along with name. In ``avg2``, the ``f`` component is squared. In ``avg3``, vector ``g`` also contributes to the average calculation.

``subroutine outpost(x,vy,vz,pr,tz,' ')``
    Dumps the current data of ``x``, ``vy``, ``vz``, ``pr``, ``tz`` to an ``.fld`` or ``.f0????`` file for post processing.

``subroutine platform_timer(ivrb)``
    Runs the battery of timing tests for matrix-matrix products,contention-free processor-to-processor ping-pong tests, and ``mpi_all_reduce`` times. Allows one to check the performance of the communication routines used on specific platforms.

``subroutine quickmv``
    Moves the mesh to allow user affine motion.

``subroutine runtimeavg(ay,y,j,istep1,ipostep,s5)``
    Computes, stores, and (for ``ipostep!0``) prints runtime averages of ``j``-quantity ``y`` (along w/ ``y`` itself unless ``ipostep<0``) with ``j`` + '``rtavg_``' + (unique) ``s5`` every ``ipostep`` for ``istep>=istep1``. ``s5`` is a string to append to ``rtavg_`` for storage file naming.

``subroutine lagrng(uo,y,yvec,uvec,work,n,m)``
    Compute Lagrangian interpolant for ``uo``

``subroutine opcopy(a1,a2,a3,b1,b2,b3)``
    Copies ``b1`` to ``a1``, ``b2`` to ``a2``, and ``b3`` to ``a3``, when ``ndim = 3``,

``subroutine cadd(a,const,n)``
    Adds ``const`` to vector ``a`` of size ``n``.

``subroutine col2(a,b,n)``
    For ``n`` entries, calculates ``a=a*b``.

``subroutine col3(a,b,c,n)``
    For ``n`` entries, calculates ``a=b*c``.

---------
Functions
---------

``function glmax(a,n)``

``function glamax(a,n)``

``function iglmax(a,n)``
    Calculates the (absolute) max of a vector that is size ``n``. Prefix ``i`` implies integer type.

``function i8glmax(a,n)``
    Calculates the max of an integer*8 vector that is size ``n``.

``function glmin(a,n)``

``function glamin(a,n)``

``function iglmin(a,n)``
    Calculates the (absolute) min of a vector that is size ``n``. Prefix ``i`` implies integer type.


``function glsc2(a,b,n)``

``function glsc3(a,b,mult,n)``

``function glsc23(z,y,z,n)``
    Performs the inner product in double precision. ``glsc3`` uses a multiplier, ``mult`` and ``glsc23`` performs ``x*x*y*z``.


``function glsum(x,n)``

``function iglsum(x,n)``

``function i8glsum(x,n)``
    Computes the global sum of ``x``, where the prefix, ``i`` specifies type integer, and ``i8`` specifies type integer*8.

---------------------------------------------------------
An Example of Specifying Surface Normals in the .usr File
---------------------------------------------------------

.. code-block:: fortran

   subroutine userbc (ix,iy,iz,iside,eg)
   include 'SIZE'
   include 'TOTAL'
   include 'NEKUSE'

   integer e,eg,f
   real snx,sny,snz   ! surface normals

   f = eface1(iside)
   e = gllel (eg)

   if (f.eq.1.or.f.eq.2) then      ! "r face"
      snx = unx(iy,iz,iside,e)                 ! Note:  iy,iz
      sny = uny(iy,iz,iside,e)
      snz = unz(iy,iz,iside,e)
   else if (f.eq.3.or.f.eq.4)  then ! "s face"
      snx = unx(ix,iz,iside,e)                 !        ix,iz
      sny = uny(ix,iz,iside,e)
      snz = unz(ix,iz,iside,e)
   else if (f.eq.5.or.f.eq.6)  then ! "t face"
      snx = unx(ix,iy,iside,e)                 !        ix,iy
      sny = uny(ix,iy,iside,e)
      snz = unz(ix,iy,iside,e)
   end if

   ux=0.0
   uy=0.0
   uz=0.0
   temp=0.0

   return
   end

This example will load a list of field files (filenames are read from a file) into the solver using the ``load_fld()`` function. After the data is loaded, the user is free to compute other postprocessing quantities. At the end the results are dumped onto a regular (uniform) mesh by a subsequent call to ``prepost()``.

Note: The regular grid data (field files) cannot be used as a restart file (uniform->GLL interpolation is unstable)!

.. code-block:: fortran

   subroutine userchk
   include 'SIZE'
   include 'TOTAL'
   include 'RESTART'

   character*80 filename(9999)

   ntot   = nx1*ny1*nz1*nelv

   ifreguo = .true.   ! dump on regular (uniform) grid instead of GLL
   nrg     = 16       ! dimension of regular grid (nrg**ndim)

   ! read file-list
   if (nid.eq.0) then
      open(unit=199,file='file.list',form='formatted',status='old')
      read(199,*) nfiles
      read(199,'(A80)') (filename(i),i=1,nfiles)
      close(199)
   end if
   call bcast(nfiles,isize)
   call bcast(filename,nfiles*80)

   do i = 1,nfiles
      call load_fld(filename(i))

      ! do something
      ! note: make sure you save the result into arrays which are
      !       dumped by prepost() e.g. T(nx1,ny1,nz1,nelt,ldimt)
      ! ...

      ! dump results into file
      call prepost(.true.,'his')
   end do

   ! we're done
   call exitt

--------------
History Points    
--------------

Assuming a case named ``blah``, a list of monitor points can be defined in file ``blah.his`` to evaluate velocity, 
temperature, pressure and passive scalars. Results will be appended to this file each time subroutine ``hpts()`` 
is called. Depending on the numnber of monitoring points you may need to increase parameter ``lhis`` in SIZE.
Usage example:

- setup an ASCII file called ``blah.his``, e.g.:

  .. code-block:: none

     3 !number of monitoring points
     1.1 -1.2 1.0
     . . .
     x y z

- add ``call hpts()`` to ``userchk()``

--------------------------
Grid-to-Grid Interpolation
--------------------------

To interpolate an existing field file (e.g. ``base.fld``) onto a new mesh do the following:

- set ``lpart`` in SIZE to a large value (e.g. 100,000 or larger) depending on your memory footprint
- compile Nek with MPI/IO support
- set ``NSTEPS=0`` in the ``.rea`` file (post-processing mode)
- run Nek using the new geometry (e.g. ``new_geom.f0000``)
- run Nek using the old geometry and add this code snippet to ``userchk()``

  .. code-block:: fortran

     character*132  newfld, oldfld, newgfld
     data newfld, oldfld, newgfld /'new0.f0001','base.fld','new_geom.f0000'/
     call g2gi(newfld, oldfld, newgfld) ! grid2grid interpolation
     call exitt()

----------------------------
Lagrangian Particle Tracking
----------------------------

The interpolation tool can be used for Lagrangian particle tracking (the particles are the interpolation points).

Workflow: Set initial particle positions (e.g. reading a file ``particle.pos0``) ``x_part <- x_pos0``

LOOP

- compute field quantities
- interpolate field quantities for all particles using ``intpts()``
- dump/store particle data
- advect particles using ``particle_advect()``

END LOOP

.. code-block:: fortranfixed

         subroutine particle_advect(rtl,mr,npart,dt_p)
   c
   c     Advance particle position in time using 4th-order Adams-Bashford.
   c     U[x\_ i(t)] for a given x\_ i(t) will be evaluated by spectral interpolation.
   c     Note: The particle timestep dt_p has be constant!
   c
         include 'SIZE'
         include 'TOTAL'

         real rtl(mr,1)

         real vell(ldim,3,lpart)  ! lagged velocities
         save vell

         integer icalld
         save    icalld
         data    icalld /0/

         if(npart.gt.lpart) then
           write(6,*) 'ABORT: npart>lpart - increase lpart in SIZE. ',nid
           call exitt
         end if

         ! compute AB coefficients (for constant timestep)
         if (icalld.eq.0) then
            call rzero(vell,3*ldim*npart) ! k = 1
            c0 = 1.
            c1 = 0.
            c2 = 0.
            c3 = 0.
            icalld = 1
         else if (icalld.eq.1) then        ! k = 2
            c0 = 1.5
            c1 = -.5
            c2 = 0.
            c3 = 0.
            icalld = 2
         else if (icalld.eq.2) then        ! k = 3
            c0 = 23.
            c1 = -16.
            c2 = 5.
            c0 = c0/12.
            c1 = c1/12.
            c2 = c2/12.
            c3 = 0.
            icalld = 3
         else                             ! k = 4
            c0 = 55.
            c1 = -59.
            c2 = 37.
            c3 = -9.
            c0 = c0/24.
            c1 = c1/24.
            c2 = c2/24.
            c3 = c3/24.
         end if

         ! compute new position x[t(n+1)]
         do i=1,npart
            do k=1,ndim
               vv = rtl(1+2*ndim+k,i)
               rtl(1+k,i) =  rtl(1+k,i) +
        &                    dt_p*(
        &                    + c0*vv
        &                    + c1*vell(k,1,i)
        &                    + c2*vell(k,2,i)
        &                    + c3*vell(k,3,i)
        &                    )
               ! store velocity history
               vell(k,3,i) = vell(k,2,i)
               vell(k,2,i) = vell(k,1,i)
               vell(k,1,i) = vv
            end do
         end do

         return
         end




