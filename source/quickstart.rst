.. _quickstart:

==============
Quickstart
==============

--------------------
Downloading the code
--------------------
All release tarballs can be found `here <https://github.com/Nek5000/Nek5000/releases>`_.
We do **not** recommend using the master branch on `GitHub <https://github.com/Nek5000/Nek5000>`_ 
in a production environment! 
 
If you use our software, please cite the following:

::

  NEK5000 Version X.Y. Release date. Argonne National Laboratory, Illinois. 
  Available: https://nek5000.mcs.anl.gov.

-------------------
Directory structure
-------------------

Here’s a brief description of each top-level directory:

.. topic:: /core

   contains the Nek5000 application sources.

.. topic:: /bin

   contains scripts for running nek5000 and manipulating its output.

.. topic:: /tools

   contains the sources for the pre- and post-processing tools which are stand-alone.

.. topic:: /short-tests

   contains light-weight regression tests for verification.
 
.. topic:: /run

   consistent place for users to place their problem cases.

.. topic:: /examples

   contains example problems.

.. topic:: /doc

   contains the user documentation.
 
.. topic:: /3rd_party

   its purpose it to provide a consistent place for 3rd party code.

---------------------
Case files
---------------------


.. topic::  SIZE

   contains some hardwired runtime parameters to dimension static arrays.

.. topic::  foo.par

   contains runtime parameters.

.. topic::  foo.re2

   contains mesh and boundary data.

.. topic::  foo.ma2

   contains partioning data.

.. topic::  foo.usr

   contains user specific code to initialize solver, set source terms and boundary conditions or to manipulate solver internals.

.. topic::  foo.his

   contains probing points.
 
.. topic::  foo.f00000

   contains checkpoint data.

.. topic::  foo.nek5000

   contains metadata for VisIt or ParaView.

.. topic::  foo.rea (legacy)

   contains runtime parameters and mesh in ASCII. Replaced by .par and .re2 file.

.. topic::  foo.map (legacy)

   contains partioning data in ASCII.

Note: The old legacy files (.rea & .map) are recommended for debugging purposes only.

-------------------
Scripts
-------------------

Let’s walk through some useful batch scripts:

- ``makenek <case>`` compiles your case
- ``nek/nekb <case>`` runs a serial job in foreground or background
- ``nekmpi/nekbmpi <case> <number of ranks>`` runs a parallel job
- ``neknek <case1> <cas2> <ranks 1> <ranks 2>`` runs Nek5000 with two overlapping component grids 
- ``visnek <case>`` creates metadata file required by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ and `ParaView <https://www.paraview.org/>`_. 
- ``mvn <old name> <new name>`` renames all case files
- ``cpn <old name> <new name>`` copies all case files

----------------------------------
Running your very first simulation
----------------------------------

Hold your horses, this needs less than 5 min.

::

  cd ~
  tar -xvzf Nek5000_X.Y.tar.gz
  export PATH=$HOME/Nek5000/bin:$PATH
  cd ~/Nek5000/tools; ./maketools genmap
  cd ~/Nek5000/run
  cp -r ../examples/eddy .
  cd eddy
  genmap                       # run partioner, on input type eddy_uv 
  makenek eddy_uv              # build case, edit script to change settings
  nekbmpi eddy_uv 2            # run Nek5000 on 2 ranks in the background
  tail logfile                 # view solver output
  visnek eddy_uv; visit -o eddy_uv.nek5000 # requires a VisIt/Paraview installation

-------------------
Meshing
-------------------

Nek5000 is mainly a solver. However, simple box type meshes can be generated with the ``genbox`` tool. For more complex meshes please consider using ``PRENEK`` and the meshing tools ``nekmerge`` and ``n2to3``. We provide mesh converters like ``exo2nek`` and ``msh2nek`` which are quite handy if you want to use your favorite mesh generator. Also check our 
`Bazaar <https://github.com/Nek5000/NekBazaar>`_ for 3rd party meshing tools.

-------------------
Visualization
-------------------
Nek5000 output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or 
`ParaView <https://www.paraview.org/>`_. This requires using ``visnek`` to generate a metadata file.  
There is also an build-in X-Window based postprocessor called ``POSTNEK`` located in tools.

-------------------
Troubleshooting
-------------------
If you run into problems compiling, installing, or running Nek5000, first check the User’s Guide. 
If you are not able to find a solution to your problem there, please send a message 
to the User’s Group `mailing list <https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users>`_.

-------------------
Reporting Bugs
-------------------
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the `Issues <https://github.com/Nek5000/Nek5000/issues>`_ feature on GitHub. 
However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble 
installing the code or getting your model to run properly, you should first send a message to the User’s Group mailing list. 
If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code,
you may create an Issue on GitHub.

-------------------
Contributing
-------------------
Our project is hosted on `GitHub <https://github.com/Nek5000>`_. Here are the most important things you need to know:

- follow the usual “fork-and-pull” Git workflow
- all development happens on the master branch
- anything in master is always deployable
- upcoming releases get their own tags out of master

If you are planning a large contribution, we encourage you to discuss the concept here on GitHub and interact with us frequently to ensure that your effort is well-directed.
