--------
reatore2
--------

The ``reatore2`` tool is a simple mesh conversion utility that will convert a legacy mesh file (``.rea``) into a binary mesh file (``.re2``).
The ``.rea`` format is a combined parameter and mesh format.
Running ``reatore2`` will generate both an ``.re2`` file and a new ``.rea`` file. 
The new ``.rea`` file will contain the same parameter settings as the previous ``.rea`` file, but will not include any element information.
Instead, the number of elements listed in the new ``.rea`` file will be negative, indicating that the mesh information is contained in a matching ``.re2`` file.
