# README #
[![Build Status](https://magnum.travis-ci.com/jgoldfar/BVP.jl.svg?token=zdgtXHoeQwBuQetRRHxV)](https://magnum.travis-ci.com/jgoldfar/BVP.jl)

### What is this repository for? ###

* BVP.jl contains implementations of BVP solvers for ODE, implemented in Julia. The primary source is Ascher, Mattheij, and Russell, *Numerical Solution of Boundary Value Problems for Ordinary Differential Equations*, SIAM, 1995.
* version 0.0.1

### How do I get set up? ###

* The main code repository is self-contained, requiring only [Julia](http://julialang.org). Run `Pkg.add("BVP")` to add the package to your Julia installation, and call `Pkg.test("BVP")`
to run the extensive set of solver tests. Alternatively, run `julia test/runtests.jl` from the current directory.

### How do I use the library? ###

Documentation TBD. Perhaps more importantly, the details of the API is TBD, but the eventual goal is to have a unified, MATLAB-like interface to the solver collection, as well as internal solvers implementing specific algorithms for specific problems.

### Who do I talk to? ###

* Jonathan Goldfarb <jgoldfar@gmail.com>
* *Note that the GitHub repo is available only as a convenience for use in Julia. Please submit issues to the main repo at [https://bitbucket.org/jgoldfar/bvp.jl](https://bitbucket.org/jgoldfar/bvp.jl)*
