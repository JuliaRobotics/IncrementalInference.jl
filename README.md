# IncrementalInference.jl
Optimization routines for incremental non-parametric and parametric solutions based on factor graphs and the Bayes (Junction) tree implemented in the ![Julia language](http://julialang.org/).

Introduction
------------

This work is an extension of the iSAM2 [1] work towards hybrid non-parametric [2] and parametric inference techniques. Inference is perfomed from a graphical model standpoint, where non-parametric and parametric solutions are based on belief propagation -- also known as the sum-product algorithm.

This code is still in development.

![alt tag](https://raw.githubusercontent.com/dehann/IncrementalInference.jl/master/doc/images/BayesTreeExample.png)

Installation
------------

Currently an unregistered package, so use:

    Pkg.clone("https://github.com/dehann/IncrementalInference.jl.git")
    
It will be useful to also install these two packages

    for p in ["Cairo"; "Fontconfig"]  Pkg.add(p) end

Basic example
=============

A multimodal 4 door robot example is available at:

    examples/RobotfourDoor.jl

References
==========

    [1]  Kaess, Michael, et al. "iSAM2: Incremental smoothing and mapping using the Bayes tree." The International Journal of Robotics Research (2011): 0278364911430419.
    [2]  Fourie, Dehann, et al. "A Nonparametric Belief Solution to the Bayes Tree." Submitted for review to IEEE IROS (2016).
