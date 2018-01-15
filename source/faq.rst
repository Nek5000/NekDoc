.. _faq:

==============
FAQ
==============

.. topic:: How can I properly reference Nek5000?

   If you use our software, please cite the following:

::

  NEK5000 Version X.Y. Release date. Argonne National Laboratory, Illinois. 
  Available: https://nek5000.mcs.anl.gov.

.. topic:: What is the license for Nek5000?

   to do

.. topic:: Where can I find support or ask questions about Nek5000?

   For any questions not answered in this document, you can email the Nek5000 community at large via the `mailing list <https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users>`_.

.. topic:: Which platforms are supported by Nek5000?

   All posix compliant operations system. Microsoft Windows is not supported.

.. topic:: Which compilers work?

   Currently GNU, Intel, PGI and IBM compilers are supported.

.. topic:: What is the difference between Pn/Pn and Pn/Pn-2?

   In brief, in the Pn/Pn formulation the pressure equation is solved on the same mesh as the velocity equations (:math:`P_n`), whereas in Pn/Pn-2 it is solved on a GLL mesh of lower order (:math:`P_{n-2}`). 
   Typically Pn/Pn-2 is more stable and slightly faster, but will yield a discontinuous pressure solution while Pn/Pn is able to better predict the wall shear.

.. topic:: What turbulence models are available in Nek5000?

    Only implicit LES (filtering or relaxtion term) is available.  

.. topic:: How do I specify/change the polynomial order?

   Change ``lx1`` in the SIZE file.

.. topic:: What element types are supported?

   Conformal curved quadrilateral/hexahedral elements.

.. topic:: How do I specify/change the solver runtime parameters?

   See ``.par`` or ``.rea`` file.

.. topic:: What polynomial order should I use?

    This depends what you are after. The sweet spot is typically :math:`N=7` (lx1=8)

.. topic:: How can I generate a mesh for use with Nek5000?

   Please see the section on `Mesh Generation and Conversion <https://nek5000.github.io/NekDoc/geometry.html>`_.

.. topic:: How do I solve for a scalar?

   Nek5000 supports solving up to 99 additional scalars.  To solve an additional scalar equation, increase ``ldimt`` in the ``SIZE`` file to accomodate the additional scalar and specify the appropriate parameter in the ``.par`` file (`for details <https://nek5000.github.io/NekDoc/user_files.html#par>`_).  

.. topic:: Where are my solution files and how do I visualize them?

   By default Nek5000 outputs solution files in the binary ``0.f%05d`` format.  These can be read by both VisIt and ParaView in conjunction with a meta-data file.  For more information see `here <https://nek5000.github.io/NekDoc/quickstart.html#visualization>`_.

.. topic:: How do I obtain values of variables at a specific point?

   to do

.. topic:: Why is ``userbc`` only getting called for certain element faces?

   ``userbc`` is ONLY called for element boundary conditions specified with a lower-case letter, e.g. 'v', 't', or 'o' but NOT 'W', 'E', or 'O'.  Note that this implies it is not necesarily called on all MPI ranks.

.. topic:: The local coordinate axes of my elements are not aligned with the global coordinate system, is this normal?

   Yes, there is no guarantee that the elements are generated with any particular orientation (except if you use genbox).
