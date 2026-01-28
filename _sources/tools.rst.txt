.. _tools:

=========
Nek Tools
=========

Nek5000* includes several tools to help with case setup.  To build or clean
them, go to the ``tools/`` directory and run ``./maketools``.  The resulting
executables will be placed in ``bin/``.

.. code-block:: console

  $ maketools --help
  Usage: maketools [clean|all|core|tool(s)]
  Examples:
    ./maketool all                    # build all tools
    ./maketool core                   # build commonly used tools
    ./maketool n2to3 reatore2         # only build listed tools
    ./maketool clean n2to3 reatore2   # only clean listed tools
    ./maketool clean                  # clean all tools

For details on each tool, see the pages below.

.. toctree::
  :glob:
  :maxdepth: 1

  exo2nek - converts EXODUS II format meshes to re2 <tools/exo2nek>
  genbox - generates simple 2D and 3D meshes <tools/genbox>
  gmsh2nek - converts GMSH meshes to re2  <tools/gmsh2nek>
  n2to3 - extrudes 2D meshes <tools/n2to3>
  reatore2 - converts ASCII (rea) meshes to binary (.re2) <tools/reatore2>

..

 * cgns2nek - converts cgns mesh files to re2
 * gencon - generates connectivity for large meshes
 * genmap - generates partitioning information
 * preNek - native *Nek5000* mesh generating tool

.. * nekmerge 
 * postNek
