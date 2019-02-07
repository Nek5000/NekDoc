===================
Additional Features
===================

Nek5000 includes several features which make case setup, running, and data processing much easier.  
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

--------------------------
Grid-to-Grid Interpolation
--------------------------
.. _features_gfldr:

Nek includes the capability to transfer a solution from one mesh to an entirely different mesh.
This allows a user to restart from an existing field file with a new mesh. 
This is accomplished by calling the generic field reader which will spectrally interpolate a result file from a previous case.
To use this feature, add a call to the ``gfldr`` subroutine in ``userchk`` in the user routines file.
For example, to interpolate a result from case ``foo``:

.. code-block:: fortran
  
  if(istep.eq.0) call gfldr("foo0.f00001")

Note that ``foo0.f00001`` must include the coordinates, i.e. it must have been created from a run with ``writeToFieldFile = yes`` in the ``[MESH]`` section of the ``.par`` file and that selection of specific fields to read is not currently supported.
That means all fields included in ``foo.f00001`` will be overwritten.

---------------
Averaging
---------------
.. _features_avg:

.. call ``avg_all`` in userchk

TODO...

------------------
Filtering
------------------
.. _features_filt:

TODO...

------------------
Objects
------------------
.. _features_obj:

TODO...


