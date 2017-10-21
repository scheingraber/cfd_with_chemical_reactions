## Group Project: Solver for Navier Stokes with transport of chemical substances and chemical reactions written in C++.

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


[![Visualization](http://archive.scheingraber.net/files/cfd.gif)](https://github.com/ChrisPara/cfd_navier_stokes_with_chemical_reaction)
