# scalarLatticeQFTcpp
Lattice Monte Carlo for Scalar Boson QFT with quartic interaction in 2-3-4 space-time dimensions. Ising model like single cluster algorithm for the sign of the field at each location. An additional update for the magnitude of the field at each lattice site. Periodic boundary conditions are applied.

For more details on physics, purpose, algorithms and implementation read the pdf report:
BaraBodur_phy781ProjectReport.pdf

Requires CERN ROOT 6 libraries
For compiling use: GNUmakefile
make qftLattice

To run ./qftLattice INPUTS
INPUTS are described in detail in qftLattice.cc
