===========================
Performance and Memory Considerations
===========================

...................
Performance Tips
...................

- Design your mesh for a polynomial order N=7 (``lx1=8`` in SIZE)
- Turn on dealiasing only if needed and try to minimize the polynomial order used for overintegration (e.g. ``lxd=10`` in SIZE)
- Ensure that you have enough elements (say >20) per MPI process
- Use AMG instead of XXT coarse grid solver
- Make sure ``usrchk()`` does not contain time consuming operations getting called in the time loop
- Enable tuned MxM implementation for your platform (see ``makenek`` options)
- Make sure residual projection for pressure is turned on and experiment with velocity projection to see if it reduced the solver time
- Tune your pressure/velocity solver tolerances (typically tolerances around 1e-5/1e-8 are acceptable)
- Try to maximize the timestep by switch in 2nd order in time and OIFS (target a Courant number of between 2-4) 
- Use binary input files e.g. ``.re2`` and ``.ma2`` to minimize solver initialization time 

...................
Memory Requirements
...................

The memory footprint of a Nek5000 run depends on many factors and is printed to
screen whenever Nek5000 exits. What follows is a rough a-priori estimate::

  lx1*ly1*lz1*lelt * 3000byte + lelg * 12byte + MPI

The memory allocated by MPI will depend heavily on the total number of
ranks and the considered MPI implementation. For large ranks counts (say > 100'000) it's
easily 50-100MB.

Note, the output of GNU`s SIZE utility is inaccurate as it does not
take into account the dynamic memory alloation of MPI, gslib, CVODE, etc. 
