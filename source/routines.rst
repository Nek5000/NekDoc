-----------
Commonly used Subroutines
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


--------------
History Points    
--------------

Assuming a case named ``foo``, a list of monitor points can be defined in file ``foo.his`` to evaluate velocity, 
temperature, pressure and passive scalars. Results will be appended to this file each time subroutine ``hpts()`` 
is called. Depending on the numnber of monitoring points you may need to increase parameter ``lhis`` in SIZE.
Usage example:

- setup an ASCII file called ``foo.his``, e.g.:

  .. code-block:: none

     3 !number of monitoring points
     1.1 -1.2 1.0
     . . .
     x y z

- add ``call hpts()`` to ``userchk()``

--------------------------
Grid-to-Grid Interpolation
--------------------------

To restart from an existing field file onto a new mesh you can call the generic field file
reader interpolation subroutine in userchk. Note that selection of specific fields to read is not currently
supported. That means all fields included in base.fld will be overwritten. Usage example:

   .. code-block:: none

      if (istep.eq.0) call gfldr('foo.f00001')

