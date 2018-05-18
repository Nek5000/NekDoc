.. _faq:

==============
FAQ
==============

--------------
General Info
--------------

**How can I properly reference Nek5000?**

   If you use our software, please cite the following:

::

  NEK5000 Version X.Y. Release date. Argonne National Laboratory, Illinois. 
  Available: https://nek5000.mcs.anl.gov.

**What is the license for Nek5000?**

   Nek5000 is licensed under BSD.  
   For more information see the ``LICENSE`` file included with the distribution in the root level directory.

**Where can I find support or ask questions about Nek5000?**

   For any questions not answered in this document, you can email the Nek5000 community at large via the `mailing list <https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users>`_.

----------------------------------
Installing, Compiling, and Running
----------------------------------

**Which platforms are supported by Nek5000?**

   All posix compliant operations system. 
   Compiling and running on Microsoft Windows 10 has been demonstrated using CMake and MinGW, although this is **not** supported.

**Which compilers work?**

   Currently GNU, Intel, PGI and IBM compilers are supported.

**When trying to recompile after increasing the number of elements I get the error below. Why does this happen and how can I avoid it?**

   ``relocation truncated to fit: R_X86_64 against symbol 'foo_' defined in COMMON section in obj/bar.o``

   This happens when the resultant executable requires more memory allocation than available, usually because ``lelt`` is too large.  
   The best way to avoid this is to increase the minimum number of MPI ranks, ``lpmin`` in ``SIZE``.  
   Alternatively, you can add ``-mcmodel=medium`` to the ``FFLAGS`` variable in ``makenek``.

-------------------
Computational Speed
-------------------

**How many elements should I have per processor?**

  This is highly dependent on the computer architecture you are using and ranges anywhere from a couple 100 in a workstation environment down to under 10 on a Blue Gene.
  This will also depend on the polynomial order, with a higher order running more efficiently on more cores.
  If possible, we recommend running a short scaling test by recording the amount of time required to run a set numer of time steps for various MPI ranks to find your particular strong-scaling limit.

.. **What is "projection" and should I use it?**

---------------------------
Problem Setup
---------------------------

**What is the difference between Pn/Pn and Pn/Pn-2?**

   In brief, in the Pn/Pn formulation the pressure equation is solved on the same mesh as the velocity equations (:math:`P_n`), whereas in Pn/Pn-2 it is solved on a GLL mesh of lower order (:math:`P_{n-2}`). 
   Typically Pn/Pn-2 is more stable and slightly faster, but will yield a discontinuous pressure solution while Pn/Pn is able to better predict the wall shear.

**What polynomial order should I use?**

   This depends what you are after. The sweet spot is typically :math:`N=7` (lx1=8)

**How do I specify/change the polynomial order?**

   Change ``lx1`` in the SIZE file.

**How do I specify/change the solver runtime parameters?**

   See ``.par`` or ``.rea`` file.

**Why is ``userbc`` only called for certain element faces?**

   ``userbc`` is ONLY called for element boundary conditions specified with a lower-case letter, e.g. 'v', 't', or 'o' but NOT 'W', 'E', or 'O'.  Note that this implies it is not necesarily called on all MPI ranks.

**How do I solve for a scalar?**

   Nek5000 supports solving up to 99 additional scalars.  
   To solve an additional scalar equation, increase ``ldimt`` in the ``SIZE`` file to accomodate the additional scalar and specify the appropriate parameter in the .:ref:`case_files_par` file.  

---------------------------
Physical Models
---------------------------

**What turbulence models are available in Nek5000?**

   LES with explicit filtering (spectral damping) or high-pass filtering is available. 
   Other turbulence models are available through user file implementation including: dynamic Smagorinsky (turbChannel example), k-ω, k-ω SST, etc. (contact the developers for more information).

-------------------
Pre-Processing
-------------------

**How can I generate a mesh for use with Nek5000?**

   Please see the section on :ref:`mesh_gen`.

**What element types are supported?**

   Conformal curved quadrilateral/hexahedral elements.

---------------
Post-Processing
---------------

**The local coordinate axes of my elements are not aligned with the global coordinate system, is this normal?**

   Yes, there is no guarantee that the elements are generated with any particular orientation (except if you use genbox).

**Where are my solution files and how do I visualize them?**

   By default Nek5000 outputs solution files in the binary ``0.f%05d`` format.  These can be read by both VisIt and ParaView in conjunction with a meta-data file.  For more information see :ref:`qstart_vis`.

**I have calculated additional fields from my solution (e.g. vorticity), how do I visualize them?**

   Using the ``.par`` file, define an additional scalar and include ``solver = none``, for example:

.. code-block:: none

   [SCALAR01] #vorticity
   solver = none

..

   Then store the calculated field in ``t(1,1,1,1,iscal+1)`` where ``iscal`` is your passive scalar index (in this example 1).

**How do I obtain values of variables at a specific point?**

  The simplest way is through the use of history points. See the section on the .:ref:`case_files_his` file.

