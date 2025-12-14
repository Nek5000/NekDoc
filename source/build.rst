=================
Build Options
=================

The shell script ``makenek`` is designed to assist the compilation process of *Nek5000*. The script will create a ``makefile`` based on the user settings section in ``makenek``. The GNU gmake utility is used to build *Nek5000*.

-----------------------
Basic ``makenek`` Usage
-----------------------

Before using ``makenek``, ensure that ``NEK_SOURCE_ROOT`` inside the script
points to the *Nek5000* source tree. Alternatively, you can set the environment
variable in your shell, for example:

.. code-block:: bash

   export NEK_SOURCE_ROOT=/home/Nek5000

If you run ``makenek`` without any arguments, it will search the current
directory for a single ``.usr`` file and use that as the case name:

.. code-block:: bash

   ./makenek

In this mode, ``makenek``

- creates a ``makefile`` based on the selected ``.usr`` file and
  configuration options, and
- immediately builds the executable ``./nek5000``.

Once the ``makefile`` has been created, you may interact directly with
``gmake`` (e.g., ``make -j8``). However, compiling via ``makenek`` is still
recommended because it ensures that the ``.usr`` file is always recompiled,
including cases where ``gmake`` fails to detect changes in dependent files such
as included headers.

If there is more than one ``.usr`` file in the directory, you must specify
the case name explicitly (without the ``.usr`` extension). For example, to
build a case defined in ``eddy_uv.usr``:

.. code-block:: bash

   ./makenek eddy_uv

---------------
Compile Options
---------------

Available configurations options:

.. _tab:bdms:

.. csv-table:: Compiler options
   :header: name,values,default,description
   :widths: 12,7,12,20

   PPLIST, string, , "list of pre-processor symbols (CVODE, ...)"
   MPI, "1, 0", 1, use MPI (needed for a multiprocessor computation)

   FC, string, optional, Fortran compiler (mpif77)
   CC, string, optional, C compiler (mpicc)
   FCLAGS, string, optional, optional Fortan compilation flags
   CCLAGS, string, optional, optional C compilation flags
   Nek5000_SOURCE_ROOT, string, optional, path of *Nek5000* source
   USR, string, optional, object list of additional files to compile make intructions (``makefile_usr.inc`` required)
   USR_LFLAGS, string, optional, optional linking flags
   PROFILING, "1, 0", 1, enable internal timers for performance statistics
   VISIT, "1, 0", 0, Toggles Visit in situ. See Visit_in_situ for details
   VISIT_INSTALL, string, VISIT in situ, Path to VISIT install path. See Visit_in_situ for details.
   VISIT_STOP, "true, false", false, "When running VISIT in situ, simulation stops after step 1 to connect VISIT."

.. _build_pplist:

-----------------
Preprocessor List
-----------------

The ``PPLIST`` field can be used to activate several features at compilation time.
A list of possible options is shown in the table below.
For a comprehensive list of options available with your version of *Nek5000*, set ``PPLIST="?"`` and execute ``makenek``.

.. _tab:PPLIST:

.. csv-table:: PPLIST options
   :header: Symbol, Description

   CVODE, compile with CVODE support for scalars
   DPROCMAP, use distributed processor mapping
   HYPRE, enable HYPRE support for AMG preconditioner (requires Cmake)
   NOMPIIO, deactivate MPI-IO support
   PARRSB, use online RSB partitioner
   VENDOR_BLAS, use VENDOR BLAS/LAPACK
   XSMM, use libxsmm for mxm

The ``HYPRE`` symbol simplifies the use of the AMG preconditioner and is necessary to use any of the HYPRE pressure preconditioners as specified in the ``.par`` file.
This will download and install `HYPRE <https://github.com/hypre-space/hypre>`_ V2.15.1.
An AMG preconditioner is necessary to run large (:math:`E>350,000`) meshes.

The ``PARRSB`` symbol compiles Nek with the new online mesh partitioner.
This can be used along with the tool ``gencon`` as an alternative to ``genmap``.
This is recommended for large (:math:`E>150,000`) meshes.

The following options have been deprecated and are no longer available or necessary as of the latest version availble on github:

.. _tab:PPLIST_dep:

.. csv-table:: legacy PPLIST options
   :header: Symbol, Description

   BGQ, use Blue Gene Q optimized mxm (supported through V19)
   EXTBAR, add underscore to exit call (supported through V17)
   NEKNEK, activate overlapping mesh solver (supported through V17)
   CMTNEK, activate discontinuous Galerkin (supported through V17)

The NekNek capability has been fully integrated into *Nek5000* and a preprocessor symbol is no longer necesary to use this feature.

.. _build_compflags:

--------------
Compiler Flags
--------------

In addition to these preprocessor items, the user can add compilation and linking flags.
``FFLAGS`` allows the user to add Fortran compilation flags while ``CCFAGS`` allows the user to
add C compilation flags.
These will be compiler dependent and the user is encouraged to consult the manual of the compiler if specific options are needed/desired.

A commonly used flag is ``-mcmodel`` which allows for arrays of size larger than 2GB.
This option  tells the compiler to use a specific memory model to generate code and store data.
It can affect code size and performance.
If your program has global and static data with a total size smaller than 2GB, ``-mcmodel=small`` is sufficient.
Global and static data larger than 2GB requires ``-mcmodel=medium`` or ``-mcmodel=large``.

.. Another useful flag is related to implicit typesetting.
.. Nek5000 relies often on implicit typesetting as default in the example cases.
.. This means in practice that if the user defines a new variable in the user file and forgets to define its type explicitly then variable beginning with a character from I to N, its type is ``INTEGER``.
.. Otherwise, it is ``REAL``.
.. To avoid confusion the user not accustomed to implicit typesetting may use the warning flag ``-Wimplicit``.
.. This flag warns whenever a variable, array, or function is implicitly declared and has an effect similar to using the ``IMPLICIT NONE`` statement in every program unit.

------------------
Cleaning the Build
------------------

To remove object files and the compiled executable using standard ``make``,
run:

.. code-block:: bash

   make clean

Alternatively, you can ask ``makenek`` to clean the build environment:

.. code-block:: bash

   ./makenek clean

which will prompt you whether to also remove third-party dependencies:

.. code-block:: text

   using source path: /home/Nek5000
   do you want to cleanup all 3rd party dependencies too? [N] y

For a more thorough cleanup, including removal of downloaded/expanded
third-party tar files, use:

.. code-block:: bash

   ./makenek distclean

which asks a similar question and removes almost everything that was
built or downloaded.
