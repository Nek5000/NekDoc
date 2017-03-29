.. _user_files:

==========
User Files
==========

Each simulation is defined by three files: the .rea file, the .usr file, and the SIZE file.  In
addition, there is a derived .map file that is generated from the .rea file by running *genmap*,
which will determine how the elements will be split across processors in the case of a parallel
run.  SIZE controls (at compile time) the polynomial degree used in the simulation, as well as the
space dimension :math:`d=2` or :math:`3`.

The SESSION.NAME file provides the name and path of the .rea file and the path to it.  It does not
however need to correspond to a .usr file of an identical name. This allows for different test
cases (.usr files) that use the same geometry and boundary conditions (.rea files).

This chapter provides an introduction to the basic files required to set up a Nek5000 simulation.

.. _user_files_session:

------------
SESSION File
------------

To run NEK5000, each simulation must have a SESSION.NAME file. This file is read in by the code and
gives the path to the relevant files describing the structure and parameters of the simulation. The
SESSION.NAME file is a file that contains the name of the simulation and the full path to
supporting files. For example, to run the eddy example from the repository, the SESSION.NAME file
would look like:

  eddy_uv
  /home/user_name/Nek5000/short_tests/eddy/ 

.. _user_files_usr:

----------------------
Case Setup File (.usr)
----------------------

.....................
Contents of .usr File
.....................


The most important interface to Nek5000 is the set of Fortran subroutines that are contained in the
*.usr* file.  This file allows direct access to all runtime variables.  Here, the user may
specify spatially varying properties (e.g., viscosity), volumetric heating sources, body forces,
and so forth.  One can also specify arbitrary initial and boundary conditions through the routines
``useric`` and ``userbc``.  The routine ``userchk`` allows the user to interrogate the
solution at the end of each timestep for diagnostic purposes.   The *.usr* files provided in
the *Nek5000/short_tests/* directory illustrate several of the more common analysis tools.  For
instance, there are utilities for computing the time average of :math:`u`, :math:`u^2`, etc. so that one
can analyze mean and rms distributions with the postprocessor.  There are routines for computing
the vorticity or the scalar :math:`\lambda_2` for vortex identification, and so forth.

.....................
Routines in .usr File
.....................



.. The routine ``uservp`` specifies the variable properties of the governing equations.  This
.. routine is called once per processor, and once per discrete point therein. 
.. 
.. \begin{tabular}{ l|l|l|l }
..    \hline
..    Equation & {\tt utrans} & {\tt udiff} & {\tt ifield} \\ \hline \hline
..    Momentum Eq.\ref{eq:ns_momentum} & \(\rho\) & \(\mu\) & 1 \\ 
..    Energy Eq.\ref{eq:energy} & \(\rho c_p\) & \(k\) & 2\\ 
..    Passive scalar Eq.\ref{eq:pass_scal} &\((\rho c_p)_i\) & \(k_i\)& i-1\\
..    \hline
.. \end{tabular}
.. 
.. 
.. \begin{lstlisting}
..       subroutine uservp (ix,iy,iz,eg)
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
..       
..       integer iel
..       iel = gllel(eg)
.. 
..       udiff =0.
..       utrans=0.
..       
..       return
..       end
.. \end{lstlisting}
.. 
.. The routine {\tt userdat} is called right after the geometry is loaded into NEK5000 and prior to
.. the distribution of the GLL points. This routine is called once per processor but for all the data
.. on that processor. At this stage the elements can be modified as long as the topology is preserved.
.. It is also possible to alter the type of boundary condition that is initially attributed in the
.. {\tt .rea} file, as illustrated below (the array {\tt cbc(face,iel,field}) contains the boundary
.. conditions per face and field of each element). Note the spacing allocated to each BC string is of
.. three units.
.. 
.. \begin{lstlisting}
..       subroutine usrdat
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
..       integer iel,f
.. 
..       do iel=1,nelt  !  Force flux BCs
..       do f=1,2*ndim
..          if (cbc(f,iel,1).eq.'W  ') cbc(f,iel,2) = 'f  ' ! flux BC for temperature
..       enddo
..       enddo
..    
..       return
..       end
.. \end{lstlisting}
.. 
.. The routine {\tt usrdat2} is called after the GLL points were distributed and allows at this point only for affine transformations of the geometry.
.. \begin{lstlisting}
..       subroutine usrdat2
..       include 'SIZE'
..       include 'TOTAL'
.. 
..       return
..       end
.. \end{lstlisting}
.. 
.. The routine {\tt userf} is called once for each point and provides the force term in Eq.\ref{eq:ns_momentum}. Not that according to the dimensionalization in Eq.\ref{eq:ns_momentum} the force term \(\vect f\) is in fact multiplied by the density \(\rho\).
.. \begin{lstlisting}
..       subroutine userf  (ix,iy,iz,eg)
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
.. 
..       ffx = 0.0
..       ffy = 0.0
..       ffz = 0.0
.. 
..       return
..       end
.. \end{lstlisting}
.. 
.. Similarly to {\tt userf} the routine {\tt userq} provides the force term in Eq.\ref{eq:energy} and the subsequent passive scalar equations according to Eq.\ref{eq:pass_scal}.
.. \begin{lstlisting}	
..       subroutine userq  (ix,iy,iz,eg)
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
..       
..       qvol   = 0.
.. 
..       return
..       end
..       \end{lstlisting}
..       
..       The boundary conditions are assigned in {\tt userbc} for both the fluid, temperature and all other scalars. An extensive list of such possible boundary conditions is available in Section.~\ref{sec:boundary}. 
..       \begin{lstlisting}
..       subroutine userbc (ix,iy,iz,iside,ieg)
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
.. 
..       ux=0.0
..       uy=0.0
..       uz=0.0
..       temp=0.0
..       flux = 1.0
..       
..       return
..       end
.. \end{lstlisting}
.. 
.. Initial conditions are attributed in {\tt useric} similarly to the boundary conditions
.. \begin{lstlisting}
..       subroutine useric (ix,iy,iz,ieg)
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
..    
..       uy=0.0
..       ux=0.0
..       uz=1.0
.. 
..       return
..       end
..       
.. \end{lstlisting}
.. The routine {\tt userchk} is called once per processor after each timestep (and once after the initialization is finished). This is the section where the solution can be interrogated and subsequent changes can be made.
.. \begin{lstlisting}
..       subroutine userchk
..       include 'SIZE'
..       include 'TOTAL'
..       include 'NEKUSE'
.. 
..       call outpost(vx,vy,vz,pr,t,'ext')
..            
..       return
..       end
..       \end{lstlisting}
..       
.. The routine {\tt usrdat3} is not widely used, however it shares the same properties with {\tt usrdat2}.
.. \begin{lstlisting}
..       subroutine usrdat3
..       include 'SIZE'
..       include 'TOTAL'
.. c
..       return
..       end
.. \end{lstlisting}
.. 
.. Nek5000 can solve the dimensional or non-dimensional equations by setting the following parameters
.. 
.. \begin{table}
.. 
.. \begin{tabular}{ l|l| }
..    \hline
..    Dimensional parameters & Non-dimensional parameters\\ \hline \hline
.. {\tt p1}=\(\rho\)      &      {\tt p1}=1\\
.. {\tt p2}=\(\nu\)       &      {\tt p2}=1/Re (-Re)\\
.. {\tt p7}=\(\rho C_p\)  &      {\tt p7}=1\\
.. {\tt p8}=\(k\)         &      {\tt p8}=1/Pe (-Pe)\\
..    \hline
.. \end{tabular}
.. \end{table}
.. 
.. alternatively the variable properties can be set in the USERVP routine.

 


------------------------
Problem-Size File (SIZE)
------------------------

-----------------------------------
Geometry and Parameters File (.rea)
-----------------------------------

-----------
Data Layout
-----------
