##############################################
Governing equations
##############################################

--------------------------------------------------------------------------------

**********************
Multicomponent mixture
**********************

Let :math:`\mathbf{v}` be the mass-averaged velocity of the mixture, defined as:

.. math::

  \mathbf{v} = \sum_{i = 1}^{n} Y_i \mathbf{u}_i

where:

- :math:`Y_i` is the mass fraction of species :math:`i`
- :math:`\mathbf{u}_i` is the velocity of species :math:`i`
- :math:`n` is the number of species (components) in the mixture.

At a given point in space and time, transport of physical quantities in a
multicomponent mixture can be described by the following set of governing
equations written in the conservative (strong) form:

Continuity equation
========================

.. math::

  \frac{\partial \rho}{\partial t} = - \nabla \cdot \rho \mathbf{v}

where:

- :math:`\rho` is the mixture density

Species mass conservation equation
=====================================

.. math::

  \frac{\partial \rho Y_i}{\partial t} = - \nabla \cdot \rho Y_i \mathbf{v} - \nabla \cdot \mathbf{j}_i + \omega_i

where:

- :math:`\mathbf{j}_i` is the mass diffusive flux of species :math:`i` relative to a mass-averaged velocity
- :math:`\omega_i` is the net mass production rate of species :math:`i`

Momentum equation
=====================================

.. math::

  \frac{\partial \rho \mathbf{v}}{\partial t} = - \nabla \cdot \rho \mathbf{v} \mathbf{v} - \nabla \cdot \pmb{\tau} - \nabla \cdot p \mathbf{I} + \rho \sum_{i=1}^{n} Y_i \mathbf{f}_i

where:

- :math:`\pmb{\tau}` is the viscous momentum flux tensor
- :math:`p` is the pressure
- :math:`\mathbf{I}` is the identity matrix
- :math:`\mathbf{f}_i` is the net acceleration from body forces acting on species :math:`i`

Total internal energy equation
=====================================

.. math::

  \frac{\partial \rho e_0}{\partial t} = - \nabla \cdot \rho e_0 \mathbf{v} - \nabla \cdot \mathbf{q} - \nabla \cdot \pmb{\tau} \cdot \mathbf{v} - \nabla \cdot p \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{n}_i

where:

- :math:`e_0` is the mixture total internal energy
- :math:`\mathbf{q}` is the heat flux
- :math:`\mathbf{n}_i` is the total mass flux of species :math:`i`

Internal energy equation
=====================================

.. math::

  \frac{\partial \rho e}{\partial t} = - \nabla \cdot \rho e \mathbf{v} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} - p \nabla \cdot \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

where:

- :math:`e` is the mixture internal energy

Enthalpy equation
=====================================

.. math::

  \frac{\partial \rho h}{\partial t} = - \nabla \cdot \rho h \mathbf{v} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} + \frac{Dp}{Dt} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

where:

- :math:`h` is the mixture enthalpy

Temperature equation
=====================================

.. math::

  \rho c_p \frac{DT}{D t} = - \nabla \cdot \mathbf{q} + \alpha T \frac{Dp}{Dt} - \pmb{\tau} : \nabla \mathbf{v} + \sum_{i=1}^{n} \big( h_i (\nabla \cdot \mathbf{j}_i - \omega_i) + \mathbf{f}_i \cdot \mathbf{j}_i \big)

where:

- :math:`c_p` is the mixture isobaric specific heat capacity
- :math:`T` is the mixture temperature
- :math:`\alpha` is the coefficient of thermal expansion
- :math:`h_i` is the enthalpy of species :math:`i`

Entropy equation
=====================================

.. math::

  \frac{\partial \rho s}{\partial t} = - \nabla \cdot \rho s \mathbf{v} - \nabla \Big( \frac{1}{T} \big( \mathbf{q} - \sum_{i=1}^{n} \tilde{\mu}_i \mathbf{j}_i \big) \Big) + \mathbf{q} \cdot \nabla \Big( \frac{1}{T} \Big) - \sum_{i=1}^{n} \mathbf{j}_i \cdot \nabla \Big( \frac{\tilde{\mu}_i}{T} \Big) - \frac{1}{T} \pmb{\tau} : \nabla \mathbf{v} + \frac{1}{T} \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i - \frac{1}{T} \sum_{i=1}^{n} \tilde{\mu}_i \omega_i

where:

- :math:`\tilde{\mu}_i` is the chemical potential of species :math:`i`
