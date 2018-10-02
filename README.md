# IncrementalInference.jl

[![Build Status](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=master)](https://travis-ci.org/JuliaRobotics/IncrementalInference.jl)
[![codecov.io](https://codecov.io/github/JuliaRobotics/IncrementalInference.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaRobotics/IncrementalInference.jl?branch=master)

[![IncrementalInference](http://pkg.julialang.org/badges/IncrementalInference_0.6.svg)](http://pkg.julialang.org/?pkg=IncrementalInference&ver=0.6)
[![IncrementalInference](http://pkg.julialang.org/badges/IncrementalInference_0.7.svg)](http://pkg.julialang.org/?pkg=IncrementalInference&ver=0.7)
[![IncrementalInference](http://pkg.julialang.org/badges/IncrementalInference_1.0.svg)](http://pkg.julialang.org/?pkg=IncrementalInference&ver=1.0)

Optimization routines for incremental non-parametric and parametric solutions based on factor graphs and the Bayes (Junction) tree implemented in the [Julia language](http://www.julialang.org/) (and [JuliaPro](http://www.juliacomputing.com)).

<p align="center">
<a href="https://vimeo.com/190052649" target="_blank"><img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/mmfgbt.gif" alt="IMAGE ALT TEXT HERE" width="480" height="320" /></a>
</p>

This package furthermore forms a cardinal piece of the [Caesar.jl](https://github.com/JuliaRobotics/Caesar.jl) robotics toolkit, including 3D visualization and database interaction, which can serve as a base station for a robotic platform. A standalone [Robot Motion Estimate](https://github.com/JuliaRobotics/RoME.jl) package is also available.

Introduction
------------

This package implements [Multi-modal iSAM [1]](http://frc.ri.cmu.edu/~kaess/pub/Fourie16iros.pdf), a descendant of the iSAM2 [3] algorithm. The main algorithm is focused towards hybrid non-parametric and parametric inference over large factor graphs. Inference is performed via the Bayes tree (similar to Junction tree) where non-parametric and parametric solutions are based on belief propagation -- also known as the sum-product algorithm.  Immediate benefits such as branch recycling is carried over as well.  Also see [related research work here [2]](https://darchive.mblwhoilibrary.org/bitstream/handle/1912/9305/Fourie_thesis.pdf?sequence=1).


The animation below shows 50% confidence lines of marginal beliefs relating to 6DOF robot poses. The approximate beliefs are being inferred through a process of successive approximation. The black trace shows the initial belief, and red the final output belief estimate. Notice the mode cycling during the process, brought about by information from elsewhere in the graph. This animation illustrates the sum-product (belief propagation) process, during the upward pass on  Bayes tree from a real data.

<p align="center">
<img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/x60mcmc.gif" alt="successive approximation" width="480"/></img>
</p>

<!-- ![alt tag](https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/BayesTreeExample.png) -->

Comments, questions and issues welcome.

Installation
------------

Pre-install the following packages system wide packages[, and easily draw factor graph and Bayes tree]:
```bash
sudo apt-get install hdf5-tools
sudo apt-get install graphviz  # optional
```

Install the package from inside Julia
```julia
Pkg.add("IncrementalInference")
```

Basic example
=============

This library is built as solver back-end which can be easily modified and extended. Specific emphasis is placed on allowing outside user defined constraint definitions to be used. The current major use case is through [RoME.jl](http://github.com/JuliaRobotics/RoME.jl) which introduces various sensor measurement and motion manifold functions for use in Robot Motion Estimate.

A few short examples, such as the multi-modal 4 door robot example, is available in the example folder:

    examples/RobotfourDoor.jl

Here 4 simultaneous modes are considered producing multi-modal posterior beliefs in the continuous domain, final consensus output and ground truth belief are show below.

<p align="center">
<img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/4doors.png" alt="Four door final result" width="640"/></img>
</p>

DataBase interaction layer
==========================

The data layer of the solver can be swapped away from the default Julia based [Graphs.jl](http://www.github.com/JuliaArchive/Graphs.jl). For using the solver on a DataBase layer please see [Caesar.jl](http://www.github.com/JuliaRobotics/Caesar.jl) and associated [CloudGraphs](http://github.com/GearsAD/CloudGraphs.jl) project.

Contributors
============

D. Fourie, M. Kaess, J. Leonard

References
==========


[1]  Fourie, Dehann, et al. ["A Nonparametric Belief Solution to the Bayes Tree"](http://www.ri.cmu.edu/pub_files/2016/10/Fourie16iros.pdf) IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), (2016).

[2]  Fourie, Dehann, ["Multi-modal and Inertial Sensor Solutions for Navigation-type Factor Graphs"](https://darchive.mblwhoilibrary.org/bitstream/handle/1912/9305/Fourie_thesis.pdf?sequence=1), Joint Program with Massachusetts Institute of Technology and Woods Hole Oceanographic Institution, Cambridge, MA, USA, August 2017.

[3]  Kaess, Michael, et al. ["iSAM2: Incremental smoothing and mapping using the Bayes tree"](http://journals.sagepub.com/doi/abs/10.1177/0278364911430419) The International Journal of Robotics Research (2011): 0278364911430419.
