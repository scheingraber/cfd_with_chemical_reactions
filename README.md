## Group Project: Solver for Navier Stokes, taking heat into account, with transport of substances and chemical reactions.

* transport of temperature (w/ Boussinesq approx. for Navier-Stokes)
* transport of chemical substances
* chemical reactions of substances

Our project simulates continuous flows using the Navier-Stokes model.
In this project, the previously implemented model was extended to take
the influence of heat into account. This is performed using the
Boussinesq approximation, which neglects changes in density.
The temperature differences directly result in the buoyancy force.
The project also models transport (i.e., diffusion and convection) of
an arbitrary number of chemical substances.
These substances, together with their reaction coefficients, are defined
inside the configuration files for the different scenarios. The initial
distribution of these substances can be defined as a constant value over
the whole domain or with the aid of a PGM-file.

#### Example: Rayleigh–Bénard convection

  A shallow/wide container, filled with a fluid at ambient temperature
  (e.g. T=0), is heated from below (Dirichlet, e.g. T=5) and has a
  fixed ambient temperature at the top (Dirichlet, e.g. T=0). The
  lateral walls are isolated (Neumann).
  In addition, there is a small drop of another fluid inside the
  container.
  
  The initially produced temperature layering is unstable and soon the
  flow evolves in regular, contra-rotating convection cells. In the
  example provided, those patterns will be stable after some period of
  time.
  The added substance can be seen to move with the flow. It is being
  slowly transported from cell to cell, which suggests, that there is
  some exchange of fluid among them.
  
![Visualization](http://archive.scheingraber.net/files/cfd.gif)
