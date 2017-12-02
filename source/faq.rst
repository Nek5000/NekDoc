.. _faq:

==============
FAQ
==============

**How can I properly reference Nek5000?**

todo

**What is the license for Nek5000?**

todo

**Where can I find support or ask questions about Nek5000?**

**Which platforms are supported by Nek5000?**

All posix compliant operations system. Microsoft Windows is not supported.

**Which compilers work?**

Currently GNU, Intel, PGI and IBM compilers are supported.

**What is the difference between Pn/Pn and Pn/Pn2?**

**What turbulence models are available in Nek5000?**

Only implicit LES (filtering or relaxtion term) is available.  

**How do I specify/change the polynomial order?**

Change ``lx1`` in the SIZE file.

**What element types are supported?**

Conformal curved quadrilateral/hexahedral elements.

**How do I specify/change the solver runtime parameters?**

See ``.par`` or ``.rea`` file.

**What polynomial order should I use?**

This depends what you are after. The sweet spot is typically :math:`N=7` (lx1=8)

**How can I generate a mesh for use with Nek5000?**

**How do I solve for a scalar?

**Where are my solution files and how do I visualize them?

**How do I obtain values of variables at a specific point?
