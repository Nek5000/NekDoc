==============
Theory
==============

.. _intro_comput_approach:

----------------------
Computational Approach
----------------------

The spatial discretization is based on the spectral element method (SEM) [Patera1984]_, which is a
high-order weighted residual technique similar to the finite element method.   In the SEM, the
solution and data are represented in terms of :math:`N` th-order tensor-product polynomials within each
of :math:`E` deformable hexahedral (brick) elements. Typical discretizations involve
:math:`E`\=100--10,000 elements of order :math:`N`\=8--16 (corresponding to 512--4096 points per
element).  Vectorization and cache efficiency derive from the local lexicographical ordering within
each macro-element and from the fact that the action of discrete operators, which nominally have
:math:`O(EN^6)` nonzeros, can be evaluated in only :math:`O(EN^4)` work and :math:`O(EN^3)` storage
through the use of tensor-product-sum factorization [Orszag1980]_.   The SEM exhibits very little
numerical dispersion and dissipation, which can be important, for example, in stability
calculations, for long time integrations, and for high Reynolds number flows. We refer to
[Denville2002]_ for more details.

Nek5000 solves the unsteady incompressible two-dimensional, axisymmetric, or three-dimensional
Stokes or Navier-Stokes equations with forced or natural convection heat transfer in both
stationary (fixed) or time-dependent geometry. It also solves the compressible Navier-Stokes in the
Low Mach regime, the magnetohydrodynamic equation (MHD).  The solution variables are the fluid
velocity :math:`\mathbf u=(u_{x},u_{y},u_{z})`, the pressure :math:`p`, the temperature :math:`T`.
All of the above field variables are functions of space :math:`{\bf x}=(x,y,z)` and time :math:`t`
in domains :math:`\Omega_f` and/or :math:`\Omega_s` defined in :numref:`fig-walls`.
Additionally Nek5000 can handle conjugate heat transfer problems.

.. _fig-walls:

.. figure:: figs/walls.png
    :align: center
    :figclass: align-center
    :alt: domains

    Computational domain showing respective fluid and solid subdomains, :math:`\Omega_f` and
    :math:`\Omega_s`.  The shared boundaries are denoted :math:`\partial\Omega_f=\partial\Omega_s`
    and the solid boundary which is not shared by fluid is :math:`\overline{\partial\Omega_s}`,
    while the fluid boundary not shared by solid :math:`\overline{\partial\Omega_f}`.

.. _intro_ns:

--------------------------------------
Incompressible Navier-Stokes Equations
--------------------------------------

The governing equations of flow motion in dimensional form are

.. math::
    :label: ns_momentum

    \rho\left(\frac{\partial\mathbf u}{\partial t} +\mathbf u \cdot \nabla \mathbf u\right) = - \nabla p + \nabla \cdot \boldsymbol\tau + \rho {\bf f} \,\, , \text{in } \Omega_f , \quad \text{  (Momentum)  } 

where :math:`\boldsymbol\tau=\mu[\nabla \mathbf u+\nabla \mathbf u^{T}]` and :math:`\mathbf f` is a user defined acceleration.

.. math::
    :label: ns_cont

    \nabla \cdot \mathbf u =0 \,\, , \text{in } \Omega_f, \quad \text{  (Continuity)  }   

If the fluid viscosity is constant in the entire domain the viscous stress tensor can be contracted
:math:`\nabla\cdot\boldsymbol\tau=\mu\Delta \mathbf u`, therefore one may solve the Navier--Stokes equations
in either the stress formulation, or no stress

- Variable viscosity requires the full stress tensor :math:`\nabla \cdot \boldsymbol\tau=\nabla \cdot
  \mu[\nabla \mathbf u+\nabla \mathbf u^{T}]`, and we shall refer to this as the stress formulation
- Constant viscosity leads to a simpler stress tensor :math:`\nabla \cdot \boldsymbol\tau=\mu\Delta \mathbf u`,
  which we refer to as the 'no stress' formulation

.. _intro_ns_nondim:

-----------------------------
Non-Dimensional Navier-Stokes
-----------------------------

Let us introduce the following non-dimensional variables :math:`\mathbf x^*\ = \frac{\mathbf x}{L}`,
:math:`\mathbf u^*\ = \frac{u}{U}`, :math:`t^*\ = \frac{tU}{L}`, and :math:`\mathbf f^* =\frac{\mathbf f L}{U^2}`.  For the pressure scale we have
two options:

- Convective effects are dominant i.e. high velocity flows :math:`p^* = \frac{p}{\rho U^2}`
- Viscous effects are dominant i.e. creeping flows (Stokes flow) :math:`p^* = \frac{p L}{\mu U}`

For highly convective flows we choose the first scaling of the pressure and obtain the
non-dimensional Navier-Stokes:

.. math::
    :label: NS_nondim

    \frac{\partial \mathbf{u^*}}{\partial t^*} + \mathbf{u^*} \cdot \nabla \mathbf{u^*}\ = -\nabla p^* + \frac{1}{Re} \nabla\cdot \boldsymbol\tau^* + \mathbf f^*.

where :math:`\boldsymbol\tau^*=[\nabla \mathbf u^*+\nabla \mathbf u^{*T}]` and :math:`\mathbf f^*` is the dimensionless user defined forcing function, e.g. gravity.

The non-dimensional number here is the Reynolds number :math:`Re=\frac{\rho U L}{\mu}`.

.. _intro_energy:

---------------
Energy Equation
---------------

In addition to the fluid flow, Nek5000 computes automatically the energy equation

.. math::
    :label: energy

    \rho c_{p} \left( \frac{\partial T}{\partial t} + \mathbf u \cdot \nabla T \right) =
       \nabla \cdot (k \nabla T) + q_{vol}\,\, ,\text{in } \Omega_f\cup \Omega_s  \text{  (Energy)  } 

.. _intro_energy_nondim:

------------------------------------------------
Non-Dimensional Energy / Passive Scalar Equation
------------------------------------------------

A similar non-dimensionalization as for the flow equations using the non-dimensional variables
:math:`\mathbf x^*\ = \frac{\mathbf x}{L}`,  :math:`\mathbf u^*\ = \frac{u}{U}`, :math:`t^*\ =
\frac{t}{L/U}`, :math:`T=\frac{T^*-T_0}{\delta T}` leads to

.. math::
    :label: energy_nondim

    \frac{\partial T^*}{\partial t^*} + \mathbf u^* \cdot \nabla T^* =
      \frac{1}{Pe} \nabla \cdot \nabla T^* + q_{vol}\,\, ,\text{in } \Omega_f\cup \Omega_s  \text{  (Energy)  } 

where :math:`Pe=LU/\alpha`, with :math:`\alpha=k/\rho c_p`.

.. _intro_pass_scal:

---------------
Passive Scalars
---------------

We can additionally solve a convection-diffusion equation for each passive scalar :math:`\phi_i`,
:math:`i = 1,2,\ldots` in :math:`\Omega_f \cup \Omega_s`

.. math::
    :label: pass_scal

    (\rho c_{p})_i \left( \frac{\partial \phi_{i}}{\partial t} + \mathbf u \cdot \nabla \phi_{i} \right) =
    \nabla \cdot (k_i \nabla \phi_{i}) + (q_{vol})_i.

The terminology and restrictions of the temperature equations are retained for the passive scalars,
so that it is the responsibility of the user to convert the notation of the passive scalar
parameters to their thermal analogues.  For example, in the context of mass transfer, the user
should recognize that the values specified for temperature and heat flux will represent
concentration and mass flux, respectively.  Any combination of these equation characteristics is
permissible with the following restrictions. First, the equation must be set to unsteady if it is
time-dependent or if there is any type of advection. For these cases, the steady-state (if it
exists) is found as stable evolution of the initial-value-problem. Secondly, the stress formulation
must be selected if the geometry is time-dependent. In addition, stress formulation must be
employed if there are traction boundary conditions applied on any fluid boundary, or if any mixed
velocity/traction boundaries, such as symmetry and outflow/n, are not aligned with either one of
the Cartesian :math:`x,y` or :math:`z` axes.  Other capabilities of Nek5000 are the linearized
Navier-Stokes for flow stability, magnetohydrodynamic flows etc.

.. _intro_ns_stokes:

---------------
Unsteady Stokes
---------------

In the case of flows dominated by viscous effects Nek5000 can solve the reduced Stokes equations

.. math::
    :label: ns_momentum_stokes

    \rho\left(\frac{\partial \mathbf u}{\partial t} \right) = - \nabla p + \nabla \cdot \boldsymbol\tau + \rho {\bf f} \,\, , \text{in } \Omega_f \text{  (Momentum)  }

where :math:`\boldsymbol\tau=\mu[\nabla \mathbf u+\nabla \mathbf u^{T}]` and

.. math::
    :label: ns_cont_stokes

    \nabla \cdot \mathbf u =0 \,\, , \text{in } \Omega_f  \text{  (Continuity)  } 

Also here we can distinguish between the stress and non-stress formulation according to whether the
viscosity is variable or not. The non-dimensional form of these equations can be obtained using the
viscous scaling of the pressure.

.. _intro_ns_steady_stokes:

-------------
Steady Stokes
-------------

If there is no time-dependence, then Nek5000 can further reduce to

.. math::
    :label: ns_momentum_steady_stokes

    - \nabla p + \nabla \cdot \boldsymbol\tau + \rho {\bf f}=0 \,\, , \text{in } \Omega_f \text{  (Momentum)  }

where :math:`\boldsymbol\tau=\mu[\nabla \mathbf u+\nabla {\mathbf u}^{T}]` and

.. math::
    :label: ns_cont_steady_stokes

    \nabla \cdot \mathbf u =0 \,\, , \text{in } \Omega_f  \text{  (Continuity)  } 

.. _intro_linear_eq:

--------------------
Linearized Equations
--------------------

In addition to the basic evolution equations described above, Nek5000 provides support for the
evolution of small perturbations about a base state by solving the *linearized equations*

.. math::
    :label: pertu

    \rho\left(\frac{\partial \mathbf u_i'}{\partial t} + \mathbf u \cdot \nabla {\mathbf u_i}^{'} + \mathbf u_i' \cdot \nabla \mathbf u \right) &=
    - \nabla p_i' + \mu \nabla^2 \mathbf u_i'\\
    \nabla \cdot \mathbf u_i' &= 0 \nonumber

for multiple perturbation fields :math:`i=1,2,\dots` subject to different initial
conditions and (typically) homogeneous boundary conditions.  

These solutions can be evolved concurrently with the base fields :math:`(\mathbf u,p,T)`.  There is
also support for computing perturbation solutions to the Boussinesq equations for natural
convection.  Calculations such as these can be used to estimate Lyapunov exponents of chaotic
flows, etc.

.. _intro_steady_conduct:

-----------------
Steady Conduction
-----------------

The energy equation :eq:`energy` in which the advection term :math:`\mathbf u \cdot \nabla T` and the
transient term :math:`\partial T/\partial t` are zero. In essence this represents a Poisson equation.

.. _intro_low_mach:

----------------------
Low-Mach Navier-Stokes
----------------------

The compressible Navier-Stokes differ mathematically from the incompressible ones mainly in the
divergence constraint :math:`\nabla \cdot \mathbf u\neq 0`. 
In this case the system of equations is not closed and an additional equation of state (EOS) is required to connect the state variables, e.g. :math:`\rho=f(p,T)`. 
Nek5000 includes the ability to solve the low-Mach approximation of the compressible Navier-Stokes, :math:`\rho\approx f(T)`. 
The low-Mach approximation decouples the pressure from the velocity leading to a system of equations which can be solved numerically in a similar fashion as the incompressible Navier-Stokes.

The low-Mach equations are 

.. math::
    :label: lowmach

    \rho\left(\frac{\partial \mathbf u}{\partial t}+ \mathbf u\cdot\nabla\mathbf u\right)&=-\nabla p+\nabla \cdot\boldsymbol\tau+\rho\mathbf f\ \\
    \nabla \cdot \mathbf u &= -\frac{1}{\rho}\frac{\mathrm d \rho}{\mathrm d T}\left(\frac{\partial T}{\partial t}+ \mathbf u\cdot\nabla T\right) \\ 
    \rho c_p\left(\frac{\partial T}{\partial t}+ \mathbf u\cdot\nabla T\right)&=-\nabla \cdot k \nabla T + q_{vol}

where :math:`\boldsymbol\tau=\mu[\nabla \mathbf u+\nabla \mathbf u^{T}-\frac{2}{3}\nabla \cdot
\mathbf u \mathbf I]`.

.. The implementation of the equation of state for the low-Mach formulation is for the moment hard-coded to be the ideal gas equation of state :math:`p=\rho R T`. 

This allows for both variable density and variable viscosity. 
The system is solved by substituting :math:`\rho\approx f(T)` into the continuity equation and obtaining a so-called thermal divergence.

.. _intro_mhd:

----------------------------
Incompressible MHD Equations
----------------------------

Magnetohydrodynamics is based on the idea that magnetic fields can induce currents in a moving
conductive fluid, which in turn creates forces on the fluid and changing the magnetic field itself.
The set of equations which describe MHD are a combination of the Navier-Stokes equations of fluid
dynamics and Maxwell's equations of electromagnetism. These differential equations have to be
solved simultaneously, and Nek5000 has an implementation for the incompressible MHD.

Consider a fluid of velocity :math:`\mathbf u` subject to a magnetic field :math:`\mathbf B` then
the incompressible MHD equations are

.. math::
    :label: mhd

    \rho\left(\frac{\partial\mathbf u}{\partial t} + \mathbf u \cdot \nabla \mathbf u\right) &= - \nabla p + \mu \Delta \mathbf u + \mathbf B\cdot \nabla \mathbf B \ ,\\ 
    \nabla \cdot \mathbf u &= 0\\ \nonumber
    \frac{\partial \mathbf B}{\partial t} + \mathbf u \cdot \nabla \mathbf B &= - \nabla q + \eta \Delta \mathbf B + \mathbf B\cdot \nabla \mathbf u \ ,\\ 
    \nabla \cdot \mathbf B &= 0 

where :math:`\rho` is the density :math:`\mu` the viscosity, :math:`\eta` resistivity, and pressure :math:`p`.

The total magnetic field can be split into two parts: :math:`\mathbf{B} = \mathbf{B_0} +
\mathbf{b}` (mean + fluctuations). The above equations become in terms of Elsässer variables
(:math:`\mathbf{z}^{\pm} =  \mathbf{u} \pm \mathbf{b}`) 

.. math::

  \frac{\partial {\mathbf{z}^{\pm}}}{\partial t}\mp\left(\mathbf {B}_0\cdot{\mathbf \nabla}\right){\mathbf z^{\pm}} + \left({\mathbf z^{\mp}}\cdot{\mathbf \nabla}\right){\mathbf z^{\pm}} = -{\mathbf \nabla}p 
  + \nu_+ \nabla^2 \mathbf{z}^{\pm} + \nu_- \nabla^2 \mathbf{z}^{\mp} 

where :math:`\nu_\pm = \nu \pm \eta`.

The important non-dimensional parameters for MHD are :math:`Re = U L /\nu` and the magnetic Re :math:`Re_M = U L /\eta`.

-----------------------------------
Arbitrary Lagrangian-Eulerian (ALE)
-----------------------------------

We consider unsteady incompressible flow in a domain with moving boundaries:

.. math::
    :label: mhd1

    \frac{\partial\mathbf u}{\partial t} = -\nabla p +\frac{1}{Re}\nabla\cdot(\nabla + \nabla^T)\mathbf u  + NL,\\
    \nabla \cdot \mathbf u = 0 

Here, :math:`NL` represents the quadratic nonlinearities from the convective term.

Our free-surface hydrodynamic formulation is based upon the arbitrary Lagrangian-Eulerian (ALE)
formulation described in [Ho1989]_.  Here, the domain :math:`\Omega(t)` is also an unknown.  As
with the velocity, the geometry :math:`\mathbf x` is represented by high-order polynomials.  For
viscous free-surface flows, the rapid convergence of the high-order surface approximation to the
physically smooth solution minimizes surface-tension-induced stresses arising from non-physical
cusps at the element interfaces, where only :math:`C^0` continuity is enforced.  The geometric
deformation is specified by a mesh velocity :math:`\mathbf w := \dot{\mathbf x}` that is
essentially arbitrary, provided that :math:`\mathbf w` satisfies the kinematic condition
:math:`\mathbf w \cdot \hat{\mathbf n}|^{}_{\Gamma} = \mathbf u \cdot \hat{\mathbf
n}|^{}_{\Gamma}`, where :math:`\hat{\mathbf n}` is the unit normal at the free surface
:math:`\Gamma(x,y,t)`.  The ALE formulation provides a very accurate description of the free
surface and is appropriate in situations where wave-breaking does not occur.

To highlight the key aspects of the ALE formulation, we introduce the weighted residual formulation
of Eq. :eq:`mhd1`: *Find* :math:`(\mathbf u,p) \in X^N \times Y^N` *such that:*

.. math::
    :label: wrt1

    \frac{\mathrm d}{\mathrm d t}(\mathbf v,\mathbf u) = (\nabla \cdot \mathbf v,p) - \frac{2}{Re}(\nabla \mathbf v,\mathbf S)
    +(\mathbf v,N\!L) + c(\mathbf v,\mathbf w,\mathbf u),
    \qquad
    (\nabla \cdot \mathbf u,q) = 0,

for all test functions :math:`(\mathbf v,q) \in X^N \times Y^N`.  Here :math:`(X^N,Y^N)` are the
compatible velocity-pressure approximation spaces introduced in [Maday1989]_, :math:`(.,.)` denotes
the inner-product :math:`(\mathbf f,\mathbf g) := \int_{\Omega(t)} \mathbf f \cdot \mathbf g \,dV`,
and :math:`\mathbf S` is the stress tensor :math:`S_{ij}^{} := \frac{1}{2}( \frac{\partial
u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} )`.  For simplicity, we have neglected the
surface tension term.  A new term in Eq.  :eq:`wrt1` is the trilinear form involving the mesh
velocity

.. math::
    :label: trilin

    c(\mathbf v,\mathbf w,\mathbf u) :=
    \int_{\Omega(t)}^{}
    \sum_{i=1}^3 
    \sum_{j=1}^3 v_i^{} \frac{\partial w_j^{} u_i^{}}{\partial x_j^{}} \,dV,

which derives from the Reynolds transport theorem when the time derivative is moved outside the
bilinear form :math:`(\mathbf v,\mathbf u_t^{})`.  The advantage of Eq. :eq:`wrt1` is that it
greatly simplifies the time differencing and avoids grid-to-grid interpolation as the domain
evolves in time.  With the time derivative outside of the integral, each bilinear or trilinear form
involves functions at a specific time, :math:`t^{n-q}`, integrated over :math:`\Omega(t^{n-q})`.
For example, with a second-order backward-difference/extrapolation scheme, the discrete form of
Eq. :eq:`wrt1` is

.. math::
    :label: bdk

    \frac{1}{2 \Delta t}\left[ 
     3 (\mathbf v^n,\mathbf u^n)^n
    -4 (\mathbf v^{n-1},\mathbf u^{n-1})^{n-1}
     + (\mathbf v^{n-2},\mathbf u^{n-2})^{n-2} \right]
    = L^n (\mathbf u) + 
    2 \widetilde{N\!L}^{n-1}
    - \widetilde{N\!L}^{n-2}.

Here, :math:`L^n(\mathbf u)` accounts for all *linear* terms in Eq. :eq:`wrt1`, including the
pressure and divergence-free constraint, which are evaluated implicitly (i.e., at time level
:math:`t^n`, on :math:`\Omega(t^n)`), and :math:`\widetilde{N\!L}^{n-q}` accounts for all 
*nonlinear* terms, including the mesh motion term :eq:`trilin`, at time-level :math:`t^{n-q}`.
The superscript on the inner-products :math:`(.,.)^{n-q}` indicates integration over
:math:`\Omega(t^{n-q})`.  The overall time advancement is as follows.  The mesh position
:math:`\mathbf x^n \in \Omega(t^n)` is computed explicitly using :math:`\mathbf w^{n-1}` and
:math:`\mathbf w^{n-2}`; the new mass, stiffness, and gradient operators involving integrals and
derivatives on :math:`\Omega(t^n)` are computed;  the extrapolated right-hand-side terms are
evaluated; and the implicit linear system is solved for :math:`\mathbf u^n`.   Note that it is only
the *operators* that are updated, not the *matrices*.  Matrices are never formed in Nek5000
and because of this, the overhead for the moving domain formulation is very low.
