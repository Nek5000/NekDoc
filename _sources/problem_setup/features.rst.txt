===================
Additional Features
===================

*Nek5000* includes several features which make case setup, running, and data processing much easier.  
These capabilities are outlined here.

.. _features_his:

--------------
History Points    
--------------

Assuming a case named ``foo``, a list of monitor points can be defined in the file ``foo.his`` to evaluate velocity, temperature, pressure and passive scalars. 
Values for each scalar will be spectrally interpolated to each point and appended to this file each time the subroutine ``hpts()`` is called. 
The values printed to this file are controlled by the ``writeToFieldFile`` parameters for each quantity in the ``.par`` file.
Depending on the number of monitoring points, you may need to increase parameter ``lhis`` in SIZE.
Note that the monitoring points are assumed to be distributed among MPI ranks, so it may be necessary to increase ``lhis`` above the actual number of history points requested.

Usage example:

- setup an ASCII file called ``foo.his``, e.g.:

  .. code-block:: none

     3 !number of monitoring points
     1.1 -1.2 1.0
     . . .
     x y z

- add ``call hpts()`` to ``userchk()``

.. _features_gfldr:

--------------------------
Grid-to-Grid Interpolation
--------------------------

Nek includes the capability to transfer a solution from one mesh to an entirely different mesh.
This allows a user to restart from an existing field file with a new mesh. 
This is accomplished by calling the generic field reader which will spectrally interpolate a result file from a previous case.
If using the latest master branch from github, this can be done by specifying the ``int`` restart option in the ``.par`` file.
For example, to interpolate a result from case ``foo``:

.. code-block:: ini

   [GENERAL]
   startFrom = foo0.f00001 int

When invoked from the ``.par`` file, the ``int`` option is compatible with all other restart options, except for ``X``.

To use this feature in V19, add a call to the ``gfldr`` subroutine in ``userchk`` in the user routines file.

.. literalinclude:: g2g.txt
   :language: fortran
   :emphasize-lines: 9

Note that ``foo0.f00001`` must include the coordinates, i.e. it must have been created from a run with ``writeToFieldFile = yes`` in the ``[MESH]`` section of the ``.par`` file and that selection of specific fields to read is not currently supported.
That means all fields included in ``foo.f00001`` will be overwritten.

.. _features_restart:

---------------
Restart Options
---------------

By default, *Nek5000* will read all available variables from the restart file. 
Restart options can be added after the filename in the ``.par`` file to control which variables are loaded, to manually set the time, and in the latest github version, to specify an interpolated restart (see :ref:`features_gfldr`).
To control which variables are read from a restart file, add an additional string composed of the following:

.. csv-table:: Variables loaded with restart options
   :header: "Variable","Characters"
   :widths: 20, 30
 
   Coordinates,``X``
   Velocity,``U``
   Pressure,``P``
   Temperature,``T``
   Passive scalar 'i',``Si``
   reset time,``time=0.0``

For example, to load only velocity and passive scalar 2, and set the physical time to 5.0, use the following

.. code-block:: ini

   [GENERAL]
   startFrom = foo0.f00001 US2 time=5.0

:Feature:
   If a restart file contains coordinates, *Nek5000* will overwrite the coordinates generated from the ``.re2`` file. This behavior may or may not be desirable, use the restart options to control it!

.. _features_avg:

---------------
Averaging
---------------

When running a high fidelity case with DNS or LES turbulence models, it is often necessary to time-average the solution fields to extract meaningful quanties.
This may sometimes even be useful for a URANS case as well.
*Nek5000* includes a subroutine for calculating a running time-average of all the primitive variables, i.e. :math:`u`, :math:`v`, :math:`w`, :math:`p`, and :math:`T`, as well as the second order terms :math:`u^2`, :math:`v^2`, :math:`w^2`, :math:`uv`, :math:`uw`, and :math:`vw`.
Which can be used to reconstruct the Reynolds stresses.
To activate time-averaging, simply call ``avg_all`` in ``userchk``.

.. literalinclude:: avgall.txt
   :language: fortran
   :emphasize-lines: 9

Adding this call to ``userchk`` will output three additional files, ``avgfoo``, ``rmsfoo``, and ``rm2foo``, where "foo" is your case name.

:WARNING:
  Averaging files are written in double precision by default and can very quickly consume a large amount of disk space!

These files will be written at the same interval as the standard restart output.
When the files are written, the averaging restarts. 
The average files thus only contain averages over the window specified by ``writeInterval`` in the ``.par`` file.

The complete list of variables, including which file they are written to and the scalar position they occupy in that file are specified in the table below.
Additionally, the width of the time-window is recorded as the phyiscal time in each average file.

.. csv-table:: Variables included in ``avg_all``
   :header: "Variable","filename","scalar"
   :widths: 10, 30, 30

   :math:`\overline{u}`,avgfoo0.f00000,u-velocity
   :math:`\overline{v}`,avgfoo0.f00000,v-velocity
   :math:`\overline{w}`,avgfoo0.f00000,w-velocity
   :math:`\overline{p}`,avgfoo0.f00000,pressure
   :math:`\overline{T}`,avgfoo0.f00000,temperature
   :math:`\overline{\phi_i}`,avgfoo0.f00000,scalar i
   :math:`\overline{u^2}`,rmsfoo0.f00000,u-velocity
   :math:`\overline{v^2}`,rmsfoo0.f00000,v-velocity
   :math:`\overline{w^2}`,rmsfoo0.f00000,w-velocity
   :math:`\overline{p^2}`,rmsfoo0.f00000,pressure
   :math:`\overline{T^2}`,rmsfoo0.f00000,temperature
   :math:`\overline{\phi_i^2}`,rmsfoo0.f00000,scalar i
   :math:`\overline{uv}`,rm2foo0.f00000,u-velocity
   :math:`\overline{vw}`,rm2foo0.f00000,v-velocity
   :math:`\overline{uw}`,rm2foo0.f00000,w-velocity
   :math:`\overline{p^2}`,rm2foo0.f00000,pressure
   :math:`\overline{T^2}`,rm2foo0.f00000,temperature

The averaging files can then be reloaded into *Nek5000* as a standard restart file for post processing.
The files contain enough information to reconstruct Reynolds stresses considering that for a sufficiently large time-window at statistically steady state:

.. math::

   \overline{u'u'}=\overline{u^2}-\overline{u}^2

Note that by default, ``avg_all`` does NOT output enough information to reconstruct the turbulent heat fluxes.
Currently, custom user code is necessary to accomplish this.

.. _features_filt:

------------------
Filtering
------------------

TODO...

.. _features_obj:

------------------
Objects
------------------
TODO...


