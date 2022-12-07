.. _tutorial_rans:

------------
RANS Channel
------------

This tutorial describes the essential setup details for RANS wall-resolved simulation, illustrated through 
a 2D channel case. The :math:`k-\tau` RANS model is employed for this tutorial, which is the recommended RANS
model in Nek5000. Other RANS models, including the regularized :math:`k-\omega`, are also available, listed 
:ref:`here <intro_ktau>`. 

..........................
Before You Begin
..........................

It is highly recommended that new users familiarize themselves with the basic Nek5000 simulation
setup files and procedures outlined in the :ref:`fdlf` and :ref:`perhill` tutorials before proceeding.

..............................
Mesh and Boundary Conditions
..............................
 
The mesh is generated with ``genbox`` using the following script

.. literalinclude:: rans/wallResolved/chan_WR.box
   :language: none 
   
It creates an infinite 2D half-channel of non-dimensional width :math:`1`. The streamwise 
(:math:`x`) direction has 3 elements with periodic boundary conditions, while the spanwise 
(:math:`y`) direction has 7 elements with symmetry boundary condition specified at the bottom face 
and wall on the top face. In addition to velocity and temperature, RANS simulation 
requires two additional scalar fields including turbulent kinetic energy, :math:`k`, and :math:`\tau`
stored at the 3rd and 4th field index respectively. Note that the 1st and 2nd index are reserved for
velocity and temperature, respectively. In wall-resolved simulations Dirichlet boundary conditions are
specified at the wall, denoted by ``W`` for velocity and ``t`` for :math:`k` and :math:`\tau`.


..............................
Control Parameters (.par file)
..............................

Details of the structure of parameter file can be found :ref:`here <case_files_par>`. For RANS simulations it
is critical to include ``[SCALAR01]`` and ``[SCALAR02]`` cards which correspond to the :math:`k` and 
:math:`\tau` fields respectively. In addition, it is essential to also include the ``[PROBLEMTYPE]`` card
and enable ``variableProperties`` and ``stressFormulation``.

For this particular tutorial, the simulation is :ref:`non-dimensionalized <intro_ns_nondim>` and flow
properties are ``density=1.0`` and ``viscosity=-125000``, the latter being the Reynolds number of the flow.
It is strongly recommended to run RANS simulations in non-dimensional form. ``density`` and ``diffusivity`` 
for ``SCALAR01`` and ``SCALAR02`` should be assigned identical values as ``density`` and ``viscosity`` for
velocity field, respectively (Note: property values for scalars are internally replaced with `velocity` 
properties).  Temperature field is not solved in this tutorial, but can be turned on by removing 
``solver=none`` entry under the ``[TEMPERATURE]`` card.

.. literalinclude:: rans/wallResolved/chan_WR.par
   :language: ini 

..............................
User Routines (.usr file)
..............................

This section describes the essential code snippets required in ``.usr`` file for RANS simulations.
Other details of all the subroutines can be found :ref:`here <case_files_usr>`. 

Foremost, it is essential to include the following header at the beginning of the ``.usr`` file.

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 2-3
      
Files in the above relative locations in the Nek5000 repo load the essential RANS subroutines.
RANS initialization is done through the ``rans_init`` subroutine call from ``usrdat2``. The 
required code snippet is shown below.

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 158-195
   
``ifld_k`` and ``ifld_omega`` variables specify the field index location of the transport variables of the 
two-equation RANS model. The specific RANS model used is identified by the ``m_id`` variable. All available
RANS models are annotated in the above code. ``m_id=4`` invokes the recommended :math:`k-\tau` model. 
``ifcoeffs=.false.`` selects the standard RANS coefficients outlined in [Wilcox2008]_. Advanced users have
the option of specifying their own coefficients which have to populated in the ``coeffs`` array, with 
``ifcoeffs=.true.``. Final parameter to be aware of is the wall distance function code ``w_id``. The recommended
value for it is ``w_id=2`` which provides a smooth distance, populated in the ``wd`` array. A cheaper option
is available through ``w_id=1`` which is recommended for periodic domains. The user also has the option of 
specifying their own wall distance array by setting ``w_id=0`` which will require ``wd`` array to
be populated with user computed wall distance before the ``rans_init`` call. 

Diffusion coefficients for all fields in RANS simulation runs must be modified to include eddy viscosity.
This is done by the following inclusions in the ``uservp`` subroutine

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 5-36
	
As above, eddy viscosity, ``mu_t``, is added to momentum diffusion coefficient and :math:`k` and :math:`\tau`
diffusion coefficients are modified as described in :eq:`ktau`. 

Source terms in the :math:`k-\tau` transport equations :eq:`ktau` are added with the following inputs in 
``userq`` subroutine

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 58-85
	
Implicit source terms contributions are specified to the ``avol`` variable, while remaining terms are 
in the ``qvol`` variable. 

For wall-resolved :math:`k-\tau` RANS, Dirichlet boundary conditions for velocity, :math:`k` and 
:math:`\tau` are all zero. ``userbc``, therefore, is simply

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 87-111
	
Initial conditions are specified in ``useric`` routine. For RANS simulations, a positive initial value for
the :math:`k` and :math:`\tau` fields is recommended. Following is used for the channel simulation,

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 113-135
	
Flow is driven through the channel by the application of streamwise mass flux, specified in ``usrdat``

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 146-156

..............................
SIZE file
..............................

The ``SIZE`` file can be copied from the available template as decribed in the :ref:`fdlf` tutorial. The user
needs to ensure that the auxiliary fields specified in the SIZE file is at minimum ``ldimt=3`` for RANS
simulations. Other details on the contents of the ``SIZE`` file can be found :ref:`here<case_files_SIZE>`.


..............................
Compilation
..............................

All required case files for RANS wall-resolved channel simulation can be downloaded using the links below:

 * :download:`chan_WR.usr <rans/wallResolved/chan_WR.usr>`
 * :download:`chan_WR.par <rans/wallResolved/chan_WR.par>`
 * :download:`chan_WR.box <rans/wallResolved/chan_WR.box>`
 * :download:`SIZE <rans/wallResolved/SIZE>`
 
..............................
Results
..............................

For reference, the results obtained from the :math:`k-\tau` RANS wall-resolved simulation are shown below and
compared with results from low-Re :math:`k-\tau` (``m_id=5``), regularized :math:`k-\omega` (``m_id=0``) and low-Re regularized 
:math:`k-\omega` (``m_id=1``) models.

.. _fig:streamwise_vel:

.. figure:: rans/U1.png
   :align: center
   :figclass: align-center

   Normalized stream-wise velocity from different RANS models.
   
.. _fig:chan_tke:

.. figure:: rans/k1.png
   :align: center
   :figclass: align-center

   Normalized TKE from different RANS models.
   
RANS models may be simply switched by using the appropriate ``m_id`` in ``usrdat2``. The low-Re versions of the model should 
be employed only for capturing the near wall TKE peak. It does not have an effect on the velocity profile. 
