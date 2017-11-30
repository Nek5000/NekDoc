===========================
Performance Considerations
===========================

- Design your mesh for a polynomial order N=7 (``lx1=8`` in SIZE)
- Turn on dealiasing only if needed and try to minimize the polynomial order used for the fine grid (e.g. ``lxd=10`` in SIZE)
- Ensure that you have at least 20 elements per MPI process (the more the better)
- Use AMG instead of XXT coarse grid solver (sweet spot depends on number of processors and elements)
- Make sure ``usrchk()`` does not contain time consuming operations getting called in the time loop
- Enable tuned MxM implementation for your platform (see ``makenek`` options)
- Make sure residual projection for pressure is turned on and experiment with velocity to see if it reduced the solver time
- Tune your pressure/velocity solver tolerances (typically 1e-4 for pressure and 1e-8 velocity are acceptable)
- Try to maximize the timestep by switch in 2nd order in time and OIFS (target a Courant number of 4 or 2) 
- ``.re2`` and ``.ma2`` are strongly recommended 
