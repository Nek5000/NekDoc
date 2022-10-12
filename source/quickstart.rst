.. _qstart:

==============
Quickstart
==============

-------------------
Before You Begin...
-------------------

It is recommended that the user have basic familiarity with Linux or another Unix-based OS.
The instructions provided in the quickstart guide and the tutorials use basic bash commands and assume the user has this knowledge.
At an absolute minimum, to successfully build and run *Nek5000* will require compatible C and Fortran compilers, such as ``gcc`` and ``gfortran``, and *GNU Make*.
To successfully build and run *Nek5000* in parallel will additionally require a compatible MPI wrapper, such as *OpenMPI* or *MPICH*.
Some of the tools and advanced features will have additional dependencies, such as *CMake*.

On Ubuntu, the following packages will provide the necessary compilers and wrappers to run *Nek5000* in parallel. 

.. code-block:: console

  $ sudo apt-get update
  $ sudo apt install build-essential  #provides gcc
  $ sudo apt install gfortran         #provides gfortran
  $ sudo apt install libopenmpi-dev   #provides openmpi

Additionally, to install all of the NekTools, the following packages are necessary

.. code-block:: console

  $ sudo apt install cmake            #provides cmake
  $ sudo apt install libx11-dev       #provides X11
  $ sudo apt install libxt-dev        #provides X11

:Note: These package lists also work for Ubuntu on the Windows Subsystem for Linux

-------------------
Directory structure
-------------------

Here’s a brief description of each top-level directory:

.. topic:: /core

   Contains the Nek5000 application sources.

.. topic:: /bin

   Contains scripts for running nek5000 and manipulating its output, and binaries for the tools. This directory should be added to your environment PATH (see :ref:`below<sec:PATH>`).

.. topic:: /tools

   Contains the sources for the pre- and post-processing tools which are stand-alone.

.. topic:: /short-tests

   Contains light-weight regression tests for verification.
 
.. topic:: /run

   A place for users to keep their problem cases. Note that many HPC systems recommend keeping source code and output on separate file systems, in which case this directory should not be used. Consult your system administrator for best practices.

.. topic:: /examples

   Contains example problems. Note that this directory is NOT included in the master branch on the GitHub repo. The *NekExamples* repository can be found `here <https://github.com/Nek5000/NekExamples>`__.

.. topic:: /3rd_party

   Contains third party software not part of the *Nek5000* core, e.g. *gslib*, *HYPRE*, and *CVODE*.

.. _sec:PATH:

--------------------
Setting up your PATH
--------------------

We recommend adding the ``bin`` directory to your shell's execution PATH.
In the ``bash`` shell, this can be done temporarily (only for your active session) with the command

.. code-block:: none

   $ export PATH+=:$HOME/Nek5000/bin

To do this more permanently, this line can be added to your ``.bashrc`` file in your ``$HOME`` directory.
This will require you to restart your current session, i.e. log out and log back in, to become active.
You can check your current execution PATH with

.. code-block:: none

  $ echo $PATH

This will produce a colon-separated list of the directories searched by Linux for the commands typed into the command line.
The ``Nek5000/bin`` entry is likely either the first or last value in this list, depending on your environment.
For more information on the execution PATH in Linux, see `here <https://opensource.com/article/17/6/set-path-linux>`__ (warning: links to a 3rd party website).

---------------------
Case files
---------------------


.. topic::  SIZE

   Contains some hardwired runtime parameters to dimension static arrays.

.. topic::  foo.par

   Contains runtime parameters.

.. topic::  foo.re2

   Contains mesh and boundary data.

.. topic::  foo.ma2

   Contains partioning data.

.. topic::  foo.usr

   Contains user specific code to initialize solver, set source terms and boundary conditions or to manipulate solver internals.

.. topic::  foo.his

   Contains probing points.
 
.. topic::  foo.f00000

   Contains checkpoint data.

.. topic::  foo.nek5000

   Contains metadata for VisIt or ParaView.

.. topic::  foo.rea (legacy)

   Contains runtime parameters and mesh in ASCII. Replaced by .par and .re2 file.

.. topic::  foo.map (legacy)

   Contains partioning data in ASCII.

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

.. code-block:: console

  $ cd ~
  $ tar -xvzf Nek5000_X.Y.tar.gz
  $ export PATH=$HOME/Nek5000/bin:$PATH
  $ cd ~/Nek5000/tools
  $ ./maketools genmap
  $ cd ~/Nek5000/run
  $ cp -r ../examples/eddy_uv .
  $ cd eddy_uv
  $ genmap                       # run partioner, on input type eddy_uv 
  $ makenek eddy_uv              # build case, edit script to change settings
  $ nekbmpi eddy_uv 2            # run Nek5000 on 2 ranks in the background
  $ tail logfile                 # view solver output
  $ visnek eddy_uv               # produces the eddy_uv.nek5000 file

As the case runs, it will generate multiple ``eddy_uv0.fXXXXX`` files.
These are the restart checkpoint and visualization data files.
The metadata file, ``eddy_uv.nek5000``, can be opened with either VisIt or ParaView, which will look for the data files in the same directory as the ``eddy_uv.nek5000`` file.

Note that this will not work if you clone the master branch from GitHub, as the ``examples`` folder is NOT included.
To obtain the examples using git, clone the ``Nek5000/NekExamples.git`` repository.

.. _qstart_meshing:

-------------------
Meshing
-------------------

*Nek5000* is mainly a solver. 
However, simple box type meshes can be generated with the ``genbox`` tool. 
For more complex meshes please consider using *preNek* and the meshing tools ``nekmerge`` and ``n2to3``. 
We provide mesh converters like ``exo2nek`` and ``gmsh2nek`` which are quite handy if you want to use your favorite mesh generator. 

.. _qstart_vis:

-------------------
Visualization
-------------------
*Nek5000* output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `ParaView <https://www.paraview.org/>`_. 
This requires using ``visnek`` to generate a metadata file.  
There is also an built-in X-Window based postprocessor called ``postnek`` located in tools.


