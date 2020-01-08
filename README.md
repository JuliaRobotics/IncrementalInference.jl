# IncrementalInference.jl

> Click on badges to follow links

Stable | Dev | Coverage | Documentation
-------|-----|----------|--------------
[![Build Status](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=master)](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl) | [![Build Status](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=release/v0.8)](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl) | [![codecov.io](https://codecov.io/github/JuliaRobotics/IncrementalInference.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaRobotics/IncrementalInference.jl?branch=master) | [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://juliarobotics.github.io/Caesar.jl/latest/)

Optimization routines for incremental non-parametric and parametric solutions based on factor graphs and the Bayes (Junction) tree implemented in the [Julia language](http://www.julialang.org/) (and [JuliaPro](http://www.juliacomputing.com)).

<p align="center">
<a href="https://vimeo.com/190052649" target="_blank"><img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/mmfgbt.gif" alt="IMAGE ALT TEXT HERE" width="480" height="320" /></a>
</p>

This package furthermore forms a cardinal piece of the [Caesar.jl](https://github.com/JuliaRobotics/Caesar.jl) robotics toolkit, including 3D visualization and database interaction, which can serve as a base station for a robotic platform. A standalone [Robot Motion Estimate](https://github.com/JuliaRobotics/RoME.jl) package is also available.


Introduction
------------

This package implements [Multi-modal iSAM [1]](http://frc.ri.cmu.edu/~kaess/pub/Fourie16iros.pdf), a descendant of the iSAM2 [3] algorithm. The main algorithm is focused towards hybrid non-parametric and parametric inference over large factor graphs. Inference is performed via the Bayes tree (similar to Junction tree) where non-parametric and parametric solutions are based on belief propagation -- also known as the sum-product algorithm.  Immediate benefits such as branch recycling is carried over as well.  Also see [related research work here [2]](https://darchive.mblwhoilibrary.org/bitstream/handle/1912/9305/Fourie_thesis.pdf?sequence=1).

Installation
------------

Pre-install the following packages system wide packages[, and easily draw factor graph and Bayes tree]:
```bash
sudo apt-get install hdf5-tools
sudo apt-get install graphviz  # optional
```

Install the package from inside Julia
```julia
(v1.0) pkg> add IncrementalInference
```

Examples
========

This library is built as solver back-end which can be easily modified and extended. Specific emphasis is placed on allowing outside user defined constraint definitions to be used. The current major use case is through [RoME.jl](http://github.com/JuliaRobotics/RoME.jl) which introduces various sensor measurement and motion manifold functions for use in Robot Motion Estimate.

A few short example is available in the example folder here.

Contributors
============

D. Fourie, M. Kaess, J. Leonard, as well as long list of Contributors in the community. 

References
==========

See [references of interest here](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/)
