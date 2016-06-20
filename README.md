# IncrementalInference.jl

[![IncrementalInference](http://pkg.julialang.org/badges/IncrementalInference_0.4.svg)](http://pkg.julialang.org/?pkg=IncrementalInference&ver=0.4)  [![codecov.io](https://codecov.io/github/dehann/IncrementalInference.jl/coverage.svg?branch=master)](https://codecov.io/github/dehann/IncrementalInference.jl?branch=master)

[![IncrementalInference](http://pkg.julialang.org/badges/IncrementalInference_0.5.svg)](http://pkg.julialang.org/?pkg=IncrementalInference&ver=0.5)

Optimization routines for incremental non-parametric and parametric solutions based on factor graphs and the Bayes (Junction) tree implemented in the [Julia language](http://www.julialang.org/).

This code is still in development (and fixing bugs as we go along).

Introduction
------------

This work is an extension of the iSAM2 [1] work towards hybrid non-parametric [2] and parametric inference techniques. Inference is performed from a graphical model standpoint, where non-parametric and parametric solutions are based on belief propagation -- also known as the sum-product algorithm.


![alt tag](https://raw.githubusercontent.com/dehann/IncrementalInference.jl/master/doc/images/BayesTreeExample.png)

Installation
------------

    Pkg.add("IncrementalInference")

It will be useful to also install these two packages

    for p in ["Cairo"; "Fontconfig"]  Pkg.add(p) end

Basic example
=============

A multimodal 4 door robot example is available at:

    examples/RobotfourDoor.jl

Which should produce maginal beliefs over all variables in the factor graphs as shown below

![alt tag](https://raw.githubusercontent.com/dehann/IncrementalInference.jl/master/doc/images/4doors.png)

References
==========

    [1]  Kaess, Michael, et al. "iSAM2: Incremental smoothing and mapping using the Bayes tree." The International Journal of Robotics Research (2011): 0278364911430419.
    [2]  Fourie, Dehann, et al. "A Nonparametric Belief Solution to the Bayes Tree." Submitted for review to IEEE IROS (2016).
