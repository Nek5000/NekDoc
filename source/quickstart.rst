.. _quickstart:

==========
Quickstart
==========

This chapter provides a quick overview to using Nek5000 for some basic flow problems distributed
with Nek5000.
           
Nek5000 runs under UNIX or any UNIX-like OS such as Linux, Mac OS/X, AIX, BG, Cray etc.  The source
is maintained in a legacy SVN repository (last updated June 2, 2016) and up-to-date Git
repositories.  Both repositories are available for download, but we recommended that new users work
with the Git repositories.  

.. _quickstart_git:

--------------------
The Git Repositories
--------------------

.. _quickstart_git_source:

___________________________
Downloading the Source Code
___________________________

.. highlight:: bash

The Nek5000 source code is hosted on GitHub at https://github.com/Nek5000/Nek5000.  You can also
download the source code as a `.zip archive
<https://github.com/Nek5000/Nek5000/archive/master.zip>`_ or a `.tar.gz archive
<https://github.com/Nek5000/Nek5000/archive/master.tar.gz>`_.  You can also download a working copy
of the repository itself.  The following commands will clone the repository in a directory named
``$HOME/Nek5000`` (where ``$HOME`` denotes your home directory)::

  cd ~
  git clone https://github.com/Nek5000/Nek5000.git -b master

You can also clone the repository in any directory of your choice::

  cd ~/my-repos/ 
  git clone https://github.com/Nek5000/Nek5000.git -b master

.. _quickstart_git_examples:

________________________________
Downloading the Example Problems
________________________________

.. highlight:: bash

The Nek5000 example problems are hosted in a second repository on GitHub,
https://github.com/Nek5000/NekExamples.  Likewise, you can download the examples as a `.zip archive
<https://github.com/Nek5000/NekExample/archive/master.zip>`_, a `.tar.gz archive
<https://github.com/Nek5000/NekExamples/archive/master.tar.gz>`_.  If you'd like to use Git, you can
download a working copy of the repository to your home directory::

  cd ~
  git clone https://github.com/Nek5000/NekExamples.git -b master

or any directory of your choice::

  cd ~/my-repos/
  git clone https://github.com/Nek5000/Nek5000.git -b master

.. _quickstart_git_tools:

_________________
Tools and Scripts
_________________

.. highlight:: bash

The ``Nek5000/tools`` directory contains programs for pre- and post-processing tasks, such as
generating meshes from geometry descriptions.  The following commands will build the tools in
``Nek5000/bin``::

  cd Nek5000/tools
  ./maketools all

By default, the ``maketools`` script will use gfortran and gcc as the Fortran and C compilers; you
can specify a different set of compilers by editing ``maketools``.  You can also build the tools
in another directory of your choice by providing a second argument to ``maketools``; for example::

  cd Nek5000/tools
  ./maketools all ~/my-tools-bin/

Besides the compiled tools, there are several convenience scripts in ``Nek5000/bin`` that allow you
to easily set up and execute runs of Nek5000.  These scripts are executable as-is and do not need
to be compiled. 

We recommend that you add ``Nek5000/bin`` to the ``$PATH`` variable in your shell's environment.  This
will allow you execute the tools and scripts without typing the full path to the executables.  For
example, In the bash shell, if you have cloned Nek5000 in ``$HOME/Nek5000`` and compiled tools
in ``$HOME/Nek5000/bin``, you can edit your executable path with::

  export PATH="$HOME/Nek5000/bin:$PATH"

In the following examples, we will assume that the tools and scripts have been added to ``$PATH``.


----------------
A Worked Example
----------------

.. highlight:: bash

As a first example, we consider the eddy problem presented by Walsh
[Walsh1992]_.  To get started, execute the following commands for the Git
repositories::

  cd $HOME/Nek5000/run
  cp -r $HOME/Nek5000/short_tests/eddy_par .
  cd eddy_par
  cp $HOME/Nek5000/bin/makenek .

Or, to run using the legacy ``.rea`` setup, execute the commands as::

  cd $HOME/Nek5000/run
  cp -r $HOME/Nek5000/short_tests/eddy .
  cd eddy
  cp $HOME/Nek5000/bin/makenek .

_________________
Modifying makenek
_________________

Makenek allows for several customization options. 
If the source is not installed in the standard directory option::

  SOURCE_ROOT="$HOME/Nek5000"  

the user can modify the variable ``SOURCE_ROOT``. 
The Fortran compiler is defined with the variable ``FC``, while the  C compiler is defined with the
 ``CC`` variable. They are set at common choices for most systems but may require to be changed 
for some specialized architectures (Blue GeneQ architectures). 

The ``PPLIST`` field can be used to activate several features at compilation time. 
A list a possible options is below: 

============= ====================================================== 
                            PPLIST Options      
--------------------------------------------------------------------   
 Symbol         Description           
============= ======================================================   
 NOMPIIO       deactivate MPI-IO support
 BGQ           use BGQ optimized mxm
 XSMM          use libxsmm for mxm
 CVODE         compile with CVODE support for scalars
 VENDOR_BLAS   use VENDOR BLAS/LAPACK
 EXTBAR        add underscore to exit call (for BGQ)
 NEKNEK        activate overlapping mesh solver (experimental)
 CMTNEK        activate DG compressible-flow solver (experimental)
============= ====================================================== 

In addition to these preprocessor items, the user can add compilation and linking flags. 
``FFLAGS`` allows the user to add Fortran compilation flags while ``CCFAGS`` allows the user to 
add C compilation flags. These will be compiler dependent and the user is encouraged to consult 
the manual of the compiler if specific options are needed/desired. 
A commonly used flag is ``-mcmodel`` which allows for arrays of size larger than 2GB. This option 
tells the compiler to use a specific memory model to generate code and store data. It can affect
code size and performance. If your program has global and static data with a total size smaller than
2GB, ``-mcmodel=small`` is sufficient. Global and static data larger than 2GB requires
``-mcmodel=medium`` or ``-mcmodel=large``. 
Another useful flag is related to implicit typesetting. Nek5000 relies often on 
implicit typesetting as default in the example cases. This means in practice that if
the user defines a new variable in the user file and forgets to define its type explicitly then
variable beginning with a character from I to N, its type is ``INTEGER``. Otherwise, it is ``REAL``.  
To avoid confusion the user not accustomed to implicit typesetting may use the warning flag 
``-Wimplicit``. This flag warns whenever a variable, array, or function is implicitly declared. 
Has an effect similar to using the ``IMPLICIT NONE`` statement in every program unit.

At linking time ``USR_LFLAGS`` can be used to add linking flags, while ``USR``
can be used to link additional files.

Nek5000 is typically run with MPI. If you do not have MPI installed on your system, uncomment 
the ``MPI=0`` flag, and change the Fortran and C compilers according to what is available  i
on your system.  In a similar manner profiling can be turned off by uncommenting ``PROFILING=0``.

Nek5000 has also some experimental features for in-situ visulization. In case you are interested 
in these features, they can be activated by uncommenting ``VISIT=1`` and follow instructions in 
makenek.  We reccomend to consult the developers before attempting to use this feature.  

_____________
Compiling Nek
_____________

.. highlight:: bash

If you have MPI installed on your system or have made the prescribed changes to makenek, the eddy
problem can be compiled as follows::

  makenek eddy_uv

If all works properly, upon comilation the executable ``nek5000`` will be generated using
``eddy_uv.usr`` to provide user-supplied initial conditions and analysis.  Note that if you
encountered a problem during a prior attempt to build the code you should type::

  makenek clean
  makenek eddy_uv

________________________
Running a case in serial
________________________

.. highlight:: bash

Once compilation is successful, start the simulation by typing::

  nekb eddy_uv

which runs the executable in the background (``nekb``, as opposed to ``nek``, which will run in the
foreground).  

__________________________
Running a case in parallel
__________________________

.. highlight:: bash

If you are running on a multi-processor machine that supports MPI, you can also run this case via::

  nekbmpi eddy_uv 4

which would run on 4 processors.    If you are running on a system that supports queuing for batch
jobs (e.g., pbs), then the following would be a typical job submission command::

  nekpbs eddy_uv 4

In most cases, however, the details of the ``nekpbs`` script would need to be modified to accommodate an
individual's user account, the desired runtime and perhaps the particular queue.   Note that the
scripts ``nek``, ``nekb``, ``nekmpi``, ``nekbmpi``, etc. perform some essential file manipulations prior to
executing ``nek5000``, so it is important to use them rather than invoking ``nek5000`` directly.

_______________________
Checking console output
_______________________

.. highlight:: bash

To check the error for this case, type::

  grep -i err eddy_uv.log | tail

or equivalently::

  grep -i err logfile | tail

where, because of the ``nekb`` script, ``logfile`` is linked to the ``.log`` file of the given
simulation.  If the run has completed, the above ``grep`` command should yield lines like

.. code-block:: none

  1000  1.000000E-01  6.759103E-05  2.764445E+00  2.764444E+00  1.000000E+00  X err
  1000  1.000000E-01  7.842019E-05  1.818632E+00  1.818628E+00  3.000000E-01  Y err

which gives for the :math:`x`- and :math:`y`-velocity components the step number, the physical time,
the maxiumum error, the maximum exact and computed values and the mean (bulk) values.

 
A common command to check on the progress of a simulation is::

  grep tep logfile | tail

which typically produces lines such as::

  Step    996, t= 9.9600000E-02, DT= 1.0000000E-04, C=  0.015 4.6555E+01 3.7611E-02

indicating, respectively, the step number, the physical time, the
timestep size, the Courant (or CFL) number, the cumulative wall clock time (in seconds)
and the wall-clock time for the most recent step.   Generally, one would 
adjust :math:`\Delta t` to have a CFL of :math:`\sim0.5`.
  
.. 
.. 
.. See Section \ref{sec:timestepping} for a comprehensive discussion of timestep selection.

----------------------------
Viewing the First 2D Example
----------------------------

.. \section{A Worked Example}
.. 

.. highlight:: bash

The preferred mode for data visualization and analysis with Nek5000 is
to use VisIt.  For a quick
peek at the data, however, we list a few commands for the native Nek5000 
postprocessor.   Assuming that the ``maketools`` script has been executed
and that ``/bin`` is in the execution path, then typing:: 

  postx 

in the working directory should open a new window with a sidebar menu.
With the cursor focus in this window (move the cursor to the window and
left click), hit ``return`` on the keyboard accept the default session name and click **plot** with the left mouse button.  This should bring up
a color plot of the pressure distribution for the first output file
from the simulation (here, ``eddy_uv.fld01``), which contains the
geometry, velocity, and pressure.  

Alternatively one can use the script *visnek*, to be found in ``/scripts``. It is sufficent to run:: 

  visnek eddy_uv

*(or the name of your session)* to obatain a file named ``eddy_uv.nek5000`` which can be recognized in `VisIt. <https://wci.llnl.gov/simulation/computer-codes/visit/>`_

.. 
.. 
.. \begin{comment}
.. To see the vorticity at the final time, load the last output file,
.. {\tt eddy\_uv.fld12}, by clicking/typing the following in the postx window:
.. \begin{tabular}{r l l l}
..   & {\bf click} \hspace*{1in} &{\bf type} \hspace*{1in} & {\bf comment} \\ \hline
.. 1.& SET TIME         & 12 & load fld12 \\
.. 2.& SET QUANTITY \\
.. 3.& VORTICITY \\
.. 4.& PLOT 
.. \end{tabular}
.. \end{comment}
.. 

**Plotting the error:**
For this case, the error has been written to 
``eddy_uv.fld11`` by making a call to ``outpost()`` in the ``userchk()``
routine in ``eddy_uv.usr``.  The error in the velocity components
is stored in the velocity-field locations and can be viewed with 
postx, or VisIt as before.

.. 
.. \begin{comment}
.. through the following sequence: 
.. \begin{tabular}{r l l l}
..   & {\bf click} \hspace*{1in} &{\bf type} \hspace*{1in} & {\bf comment} \\ \hline
.. 1.& SET TIME         & 11 & load fld11 \\
.. 2.& SET QUANTITY \\
.. 3.& VELOCITY \\
.. 4.& MAGNITUDE \\
.. 5.& PLOT  \\
.. \end{tabular}
.. \end{comment}
.. 
.. \subsection{Modifying the First Example}
.. 

---------------------------
Modifying the First Example
---------------------------

_______________________
Using the ``.par`` file
_______________________

.. highlight:: bash

A common step in the Nek5000 workflow is to rerun with a higher
polynomial order.   Typically, one runs a relatively low-order case
(e.g., ``lx1`` =5) for one or two flow-through times and then uses
the result as an initial condition for a higher-order run
(e.g., ``lx1`` =8).  We illustrate the procedure with the 
``eddy_uv`` example.
.. 
Assuming that the contents of ``Nek5000/bin``
are in the execution path, begin by typing::

  cpn eddy_uv eddy_new

which will copy the requisite ``eddy_uv`` case files
to ``eddy_new``.  
Next, edit ``SIZE`` and change the two lines defining
``lx1`` and ``lxd`` from::

       parameter (lx1=8)                ! p-order (avoid uneven and values <6)
       parameter (lxd=12)               ! p-order for over-integration (dealiasing)

to::

       parameter (lx1=12)                ! p-order (avoid uneven and values <6)
       parameter (lxd=18)               ! p-order for over-integration (dealiasing)

Then recompile the source by typing::

  makenek eddy_new

Next, edit ``eddy_new.par`` and add the following line in the ``GENERAL``

.. code-block:: none

     startFrom = eddy_uv.fld12

which tells Nek5000 to use the contents of ``eddy_uv.fld12``
as the initial condition for ``eddy_new``.
The simulation is started in the usual way::

  nekb eddy_new

after which the command::

  grep err logfile | tail

will show a much smaller error (:math:`\sim 10^{-9}`) than the ``lx1=8``
case. 

Note that one normally would not use a restart file for the *eddy*
problem, which is really designed as a convergence study.  The purpose here, however, was two-fold, namely,
to illustrate a change of order and its impact on the error, and to
demonstrate the frequently-used restart procedure. However for a higher order timestepping scheme an accurate restart would require a number of field files of the same size (+1) as the order of the multistep scheme


_______________________
Using the ``.rea`` file
_______________________

.. highlight:: bash

Modifying the legacy ``.rea`` and ``SIZE`` file formats is done in a similar way to above. In the lgeacy ``SIZE`` file, change::

       parameter (lx1=8,ly1=lx1,lz1=1,lelt=300,lelv=lelt)
       parameter (lxd=12,lyd=lxd,lzd=1)

to::

       parameter (lx1=12,ly1=lx1,lz1=1,lelt=300,lelv=lelt)
       parameter (lxd=18,lyd=lxd,lzd=1)

Then recompile the source by typing::
  
  makenek eddy_new


Next, edit the legacy ``.rea`` format, edit ``eddy_new.rea`` and change the line

.. code-block:: none
 
             0 PRESOLVE/RESTART OPTIONS  *****

(found roughly 33 lines from the bottom of the file) to

.. code-block:: none

             1 PRESOLVE/RESTART OPTIONS  *****
  eddy_uv.fld12

which tells Nek5000 to use the contents of ``eddy_uv.fld12``
as the initial condition for ``eddy_new``.

The simulation can now be run in a similar way to the above ``.par`` section.

_________________
The ``.par`` File
_________________

In the future of Nek5000, the ``.rea`` file format will be replaced by the new ``.par`` file format. This is because the new format is simpler to read (as can be seen in the *eddy* problem). Here, we introduce the basics of the ``.par`` file as used in this example case.

In the ``short_tests/eddy_par`` directory, there is a file labeled ``eddy_uv.par``, which is broken into four sections: ``GENERAL``, ``PROBLEMTPYE``, ``PRESSURE``, and ``VELOCITY``. Each of these sections contains necessary and optional parameters for Nek5000 to run that the user inputs, similar to the legacy ``.rea`` format. Analogs to each of the parameters can be found in ``eddy_uv.rea``, along with descriptions. As well, detailed analog information can be found in the ``core/readat_new.f`` file (where each numbered parameter is matched with its new name in the ``.par`` file).

.. code-block:: none

    #
    # nek parameter file
    #
    [GENERAL]
    #startFrom = restart.fld
    stopAt = numSteps #endTime
    numSteps = 1000

    dt = -1E-04
    timeStepper = bdf #char #steady
    tOrder = 3

    writeControl = timeStep #runTime
    writeInterval = 100

    dealiasing = yes

    maxCFL = 1

    userParam01 = 1 #u0 transational velocity (formerly param(96))
    userParam02 = 0.3 #v0 translational velocity (formerly param(97))

    [PROBLEMTYPE]
    #stressFormulation = yes

    [PRESSURE]
    preconditioner = schwarz #semg #amg
    residualTol = 1E-09
    residualProj = yes

    [VELOCITY]
    residualTol = 1E-12
    residualProj = yes
    density = 1
    viscosity = -20
