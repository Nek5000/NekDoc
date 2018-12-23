.. _case_files_his:

--------------
History Points    
--------------

Assuming a case named ``foo``, a list of monitor points can be defined in file ``foo.his`` to evaluate velocity, temperature, pressure and passive scalars. 
Values for each scalar will be spectrally interpolated to each point and appended to this file each time the subroutine ``hpts()`` is called. 
Depending on the number of monitoring points you may need to increase parameter ``lhis`` in SIZE.
Usage example:

- setup an ASCII file called ``foo.his``, e.g.:

  .. code-block:: none

     3 !number of monitoring points
     1.1 -1.2 1.0
     . . .
     x y z

- add ``call hpts()`` to ``userchk()``

