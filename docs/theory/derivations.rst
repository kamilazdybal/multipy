##############################################
Derivations
##############################################

This section documents few less commonly presented derivations.

.. contents:: Table of contents
   :depth: 5

--------------------------------------------------------------------------------

********************************************
Governing equations
********************************************

Summing over all :math:`n` species mass conservation equations yields the continuity equation
================================================================================================

If we sum over all :math:`n` species mass conservation equations, we should
arrive at the continuity equation for the entire mixture. Below, we show that
for the Eulerian differential and integral forms of the species mass conservation
equation.

The Eulerian **differential** form of the species mass conservation equation
--------------------------------------------------------------------------------

We begin with:

.. math::

  \frac{\partial \rho Y_i}{\partial t} = - \nabla \cdot \mathbf{n}_i + \omega_i

Summing over all :math:`n` equations:

.. math::

  \sum_{i=1}^{n} \frac{\partial \rho Y_i}{\partial t} = - \sum_{i=1}^{n} \nabla \cdot \mathbf{n}_i + \sum_{i=1}^{n} \omega_i

Since :math:`\sum_{i=1}^{n} \omega_i = 0`:

.. math::

  \sum_{i=1}^{n} \frac{\partial \rho Y_i}{\partial t} = - \sum_{i=1}^{n} \nabla \cdot \mathbf{n}_i

Which we can also write as:

.. math::

  \frac{\partial}{\partial t} \sum_{i=1}^{n} \rho Y_i = - \nabla \cdot \sum_{i=1}^{n}  \mathbf{n}_i

The density :math:`\rho` can now be taken out of the summation. Since
:math:`\sum_{i=1}^{n} Y_i = 1` and :math:`\sum_{i=1}^{n} \mathbf{n}_i = \mathbf{n}_t`:

.. math::

  \frac{\partial \rho}{\partial t} = - \nabla \cdot \mathbf{n}_t

which is the continuity equation in the differential form that we were after.
In addition, we may express the total fluxes as: :math:`\mathbf{n}_t = \rho \mathbf{v}` which gives:

.. math::

  \frac{\partial \rho}{\partial t} = - \nabla \cdot \rho \mathbf{v}

The Eulerian **integral** form of the species mass conservation equation
--------------------------------------------------------------------------------

We begin with:

.. math::

  \int_{V(t)} \omega_i dV = \int_{V(t)} \frac{\partial \rho Y_i}{\partial t} dV + \int_{S(t)} \mathbf{n}_i \cdot \mathbf{a} dS

Summing over :math:`n` equations:

.. math::

  \sum_{i=1}^{n} \int_{V(t)} \omega_i dV = \sum_{i=1}^{n} \int_{V(t)} \frac{\partial \rho Y_i}{\partial t} dV + \sum_{i=1}^{n} \int_{S(t)} \mathbf{n}_i \cdot \mathbf{a} dS

which we can write as:

.. math::

  \int_{V(t)} \sum_{i=1}^{n} \omega_i dV = \int_{V(t)}  \sum_{i=1}^{n} \frac{\partial \rho Y_i}{\partial t} dV + \int_{S(t)}  \sum_{i=1}^{n} \mathbf{n}_i \cdot \mathbf{a} dS

We can thus use the same reasoning for terms :math:`\sum_{i=1}^{n} \omega_i`
and :math:`\sum_{i=1}^{n} \frac{\partial \rho Y_i}{\partial t}` as previously to get:

.. math::

  0 = \int_{V(t)} \frac{\partial \rho}{\partial t} dV + \int_{S(t)}  \sum_{i=1}^{n} \mathbf{n}_i \cdot \mathbf{a} dS

Finally, we need to think whether :math:`\sum_{i=1}^{n} (\mathbf{n}_i \cdot \mathbf{a}) = \mathbf{n}_t \cdot \mathbf{a}`?
Once we obtain a discretization of surface :math:`S(t)` into multiple :math:`dS`,
each :math:`dS` has got a fixed normal vector, :math:`\mathbf{a}`, associated with it.
Since the same Eulerian boundary encapsulates the total flux of species :math:`i`
in the above equation and the total flux of the mixture, we can conclude that
:math:`\mathbf{a}` is independent of :math:`i` for a given :math:`S(t)`
and pull it out of the summation. This gives us:

.. math::

  \sum_{i=1}^{n} (\mathbf{n}_i \cdot \mathbf{a}) = \mathbf{a} \cdot \sum_{i=1}^{n} \mathbf{n}_i = \mathbf{a} \cdot \mathbf{n}_t

which brings us to the continuity equation in the integral form:

.. math::

  \int_{V(t)} \frac{\partial \rho}{\partial t} dV  = - \int_{S(t)}   \mathbf{n}_t \cdot \mathbf{a} dS

We may at this point substitute :math:`\mathbf{n}_t = \rho \mathbf{v}` and use the divergence theorem to get:

.. math::

  \int_{V(t)} \frac{\partial \rho}{\partial t} dV  = - \int_{V(t)} \nabla \cdot \rho \mathbf{v} dS

--------------------------------------------------------------------------------

Total energy equation from the momentum equation
================================================================================

We begin with the total internal energy equation:

.. math::

  \frac{\partial \rho e_0}{\partial t} = - \nabla \cdot \rho e_0 \mathbf{v} - \nabla \cdot \mathbf{q} - \nabla \cdot \pmb{\tau} \cdot \mathbf{v} - \nabla \cdot p \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{n}_i

and with the momentum equation:

.. math::

  \frac{\partial \rho \mathbf{v}}{\partial t} = - \nabla \cdot \rho \mathbf{v} \mathbf{v} - \nabla \cdot \pmb{\tau} - \nabla p + \rho \sum_{i=1}^{n} Y_i \mathbf{f}_i

We multiply the momentum equation by :math:`\mathbf{v}` to obtain the kinetic energy equation:

.. math::

  \boxed{\mathbf{v} \cdot \frac{\partial \rho \mathbf{v}}{\partial t} = - \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} \mathbf{v} - \mathbf{v} \cdot \nabla \cdot \pmb{\tau} - \mathbf{v} \cdot \nabla p + \mathbf{v} \rho \sum_{i=1}^{n} Y_i \mathbf{f}_i}

We note that the total internal energy is equal to the internal energy and the kinetic energy:

.. math::

  e_0 = e + \frac{1}{2} \mathbf{v} \cdot \mathbf{v}

We first substitute the expression for :math:`e_0` in the total internal energy equation:

.. math::

  \boxed{\frac{\partial \rho (e + \frac{1}{2} \mathbf{v} \cdot \mathbf{v})}{\partial t} = - \nabla \cdot \rho (e + \frac{1}{2} \mathbf{v} \cdot \mathbf{v}) \mathbf{v} - \nabla \cdot \mathbf{q} - \nabla \cdot \pmb{\tau} \cdot \mathbf{v} - \nabla \cdot p \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{n}_i}

We then subtract the kinetic energy equation from the total internal energy equation, since :math:`e = e_0 - \frac{1}{2} \mathbf{v} \cdot \mathbf{v}`.

LHS
---

LHS after subtraction is:

.. math::

  \frac{\partial \rho (e + \frac{1}{2} \mathbf{v} \cdot \mathbf{v})}{\partial t} - \mathbf{v} \cdot \frac{\partial \rho \mathbf{v}}{\partial t}

Using the chain rule on all terms:

.. math::

  \frac{\partial \rho e}{\partial t} + \frac{1}{2} \rho \mathbf{v} \frac{\partial \mathbf{v}}{\partial t} + \frac{1}{2} \rho \mathbf{v} \frac{\partial \mathbf{v}}{\partial t} + \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t} - \rho \mathbf{v} \cdot \frac{\partial \mathbf{v}}{\partial t} - \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t}

This is equal to:

.. math::

  \frac{\partial \rho e}{\partial t} - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t}

RHS
---

We will split RHS after subtraction into few groups of terms.
The first group are terms involving density:

- :math:`- \nabla \cdot \rho (e + \frac{1}{2} \mathbf{v} \cdot \mathbf{v}) \mathbf{v} + \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} \mathbf{v} = - \nabla \cdot \rho e \mathbf{v} - \nabla \cdot \frac{1}{2} \rho \mathbf{v} \cdot \mathbf{v} \mathbf{v} + \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} \mathbf{v}`

The second group is the heat flux term:

- :math:`- \nabla \cdot \mathbf{q}`

The third group are terms involving the viscous molecular flux tensor:

- :math:`- \nabla \cdot \pmb{\tau} \cdot \mathbf{v} + \mathbf{v} \cdot \nabla \cdot \pmb{\tau} = - \pmb{\tau} : \nabla \mathbf{v} - \mathbf{v} \cdot \nabla \cdot \pmb{\tau} + \mathbf{v} \cdot \nabla \cdot \pmb{\tau} = - \pmb{\tau} : \nabla \mathbf{v}`

The fourth group are terms involving pressure:

- :math:`- \nabla \cdot p \mathbf{v} + \mathbf{v} \cdot \nabla p = - p \nabla \cdot \mathbf{v} - \mathbf{v} \cdot \nabla p + \mathbf{v} \cdot \nabla p = - p \nabla \cdot \mathbf{v}`

The fifth group are the terms involving body forces:

- :math:`\sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{n}_i - \mathbf{v} \rho \sum_{i=1}^{n} Y_i \mathbf{f}_i = \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{n}_i - \sum_{i=1}^{n} \rho Y_i \mathbf{v} \cdot \mathbf{f}_i = \sum_{i=1}^{n} \underbrace{(\mathbf{n}_i - \rho_i \mathbf{v})}_\text{diffusive flux of $i$} \cdot \mathbf{f}_i = \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i`

LHS with RHS
------------

We now put LHS and RHS together to get:

.. math::

  \frac{\partial \rho e}{\partial t} \underbrace{- \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t}} = - \nabla \cdot \rho e \mathbf{v} \underbrace{- \nabla \cdot \frac{1}{2} \rho \mathbf{v} \cdot \mathbf{v} \mathbf{v}} + \underbrace{\mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} \mathbf{v}} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} - p \nabla \cdot \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

We tackle the underbraced terms below:

.. math::

  - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t} = - \nabla \cdot \frac{1}{2} \rho \mathbf{v} \cdot \mathbf{v} \mathbf{v} + \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} \mathbf{v}

Applying the chain rule, we get:

.. math::

  - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t} = - \frac{1}{2} \rho \mathbf{v} \cdot \nabla \cdot \mathbf{v} \mathbf{v} - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} + \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} + \rho \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \mathbf{v}

Applying the chain rule one more time on the :math:`- \frac{1}{2} \rho \mathbf{v} \cdot \nabla \cdot \mathbf{v} \mathbf{v}` term we get :math:`- \rho \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \mathbf{v}`. This term cancels out with the :math:`\rho \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \mathbf{v}` term. We are thus left with:

.. math::

  - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t} = - \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} + \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} = \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v}

Collecting all terms on one side and factoring out :math:`\mathbf{v} \cdot \mathbf{v}` we get:

.. math::

  \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \frac{\partial \rho}{\partial t} + \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \cdot \nabla \cdot \rho \mathbf{v} = \frac{1}{2} \mathbf{v} \cdot \mathbf{v} \underbrace{\Big( \frac{\partial \rho}{\partial t} + \nabla \cdot \rho \mathbf{v} \Big)}_\text{equal to 0 from continuity} = 0

Finally, the internal energy equation is:

.. math::

  \boxed{\frac{\partial \rho e}{\partial t} = - \nabla \cdot \rho e \mathbf{v} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} - p \nabla \cdot \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i}

--------------------------------------------------------------------------------

Enthalpy equation from the internal energy equation
================================================================================

Now that we've obtained the internal energy equation, we can use relation:

.. math::

  h = e + \frac{p}{\rho}

to obtain the enthalpy equation. Multiplying the above by :math:`\rho`:

.. math::

  \rho h = \rho e + p

and differentiating with respect to time:

.. math::

  \frac{\partial \rho h}{\partial t} = \frac{\partial \rho e}{\partial t} + \frac{\partial p}{\partial t}

If we insert the two above relationships into the internal energy equation, we get:

.. math::

  \frac{\partial \rho h}{\partial t} \underbrace{- \frac{\partial p}{\partial t}} = \underbrace{- \nabla \cdot (\rho h - p) \mathbf{v}} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} \underbrace{- p \nabla \cdot \mathbf{v}} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

With some re-arrangements of the underbraced terms, this is equivalent to:

.. math::

  \frac{\partial \rho h}{\partial t} = \underbrace{- \nabla \cdot \rho h \mathbf{v} + \nabla \cdot p \mathbf{v} - p \nabla \cdot \mathbf{v} + \frac{\partial p}{\partial t}} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

Applying the chain rule on the :math:`\nabla \cdot p \mathbf{v}` term:

.. math::

  \frac{\partial \rho h}{\partial t} = - \nabla \cdot \rho h \mathbf{v} + p \nabla \cdot \mathbf{v} + \mathbf{v} \cdot \nabla p - p \nabla \cdot \mathbf{v} + \frac{\partial p}{\partial t} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i

We can cancel out the :math:`p \nabla \cdot \mathbf{v}` term and write :math:`\frac{\partial p}{\partial t} + \mathbf{v} \cdot \nabla p = \frac{D p}{D t}`.

The enthalpy equation thus becomes:

.. math::

  \boxed{\frac{\partial \rho h}{\partial t} = - \nabla \cdot \rho h \mathbf{v} + \frac{D p}{D t} - \nabla \cdot \mathbf{q} - \pmb{\tau} : \nabla \mathbf{v} + \sum_{i=1}^{n} \mathbf{f}_i \cdot \mathbf{j}_i}
