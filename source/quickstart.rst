.. _qstart:

==============
Quickstart
==============

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

   contains example problems. Note: NOT included in the master branch on the GitHub repo.

.. topic:: /doc

   contains the user documentation. Note: NOT included in the master branch on the GitHub repo.
 
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
Begin by downloading the latest release tarball from `here <https://github.com/Nek5000/Nek5000/releases>`_.
Then follow the instructions below

::

  cd ~
  tar -xvzf Nek5000_X.Y.tar.gz
  export PATH=$HOME/Nek5000/bin:$PATH
  cd ~/Nek5000/tools; ./maketools genmap
  cd ~/Nek5000/run
  cp -r ../examples/eddy_uv .
  cd eddy_uv
  genmap                       # run partioner, on input type eddy_uv 
  makenek eddy_uv              # build case, edit script to change settings
  nekbmpi eddy_uv 2            # run Nek5000 on 2 ranks in the background
  tail logfile                 # view solver output
  visnek eddy_uv; visit -o eddy_uv.nek5000 # requires a VisIt/Paraview installation

Note that this will not work if you clone the master branch from GitHub, as the ``examples`` folder is NOT included.

.. _qstart_meshing:

-------------------
Meshing
-------------------

Nek5000 is mainly a solver. However, simple box type meshes can be generated with the ``genbox`` tool. For more complex meshes please consider using ``PRENEK`` and the meshing tools ``nekmerge`` and ``n2to3``. We provide mesh converters like ``exo2nek`` and ``msh2nek`` which are quite handy if you want to use your favorite mesh generator. Also check our 
`Bazaar <https://github.com/Nek5000/NekBazaar>`_ for 3rd party meshing tools.

.. _qstart_vis:

-------------------
Visualization
-------------------
Nek5000 output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or 
`ParaView <https://www.paraview.org/>`_. This requires using ``visnek`` to generate a metadata file.  
There is also an build-in X-Window based postprocessor called ``POSTNEK`` located in tools.


