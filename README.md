# README #

### What is this repository for? ###

* BVP.jl contains implementations of BVP solvers for ODE, implemented in Julia. The primary source is Ascher, Mattheij, and Russell, *Numerical Solution of Boundary Value Problems for Ordinary Differential Equations*, SIAM, 1995.
* version 0.1

### How do I get set up? ###

* The main code repository is self-contained, requiring only [Julia](http://julialang.org). After cloning the package to your `~/.julia/v0.3` or `~/.julia/v0.4` directory, running
```julia
Pkg.test("BVP")
```
will run the extensive set of solver tests. Alternatively, run `julia test/runtests.jl` from the current directory.

### How do I use the library? ###

Documentation TBD. Perhaps more importantly, the details of the API is TBD, but the eventual goal is to have a unified, MATLAB-like interface to the solver collection, as well as special

### Who do I talk to? ###

* Jonathan Goldfarb <jgoldfar@gmail.com>
