# IncrementalInference.jl

> Click on badges to follow links:

| Stable Release | Dev branch | Coverage | Documentation |
|----------------|------------|----------|---------------|
| [![iif-ci-stb][iif-ci-stb-img]][iif-ci-stb-url] <br> [![version][iif-ver-img]][iif-rel-url] | [![iif-ci-dev-img]][iif-ci-dev-url] <br> [![iif-commits-url]][contributors-url] <br> [![issues-time]][issues-url] | [![doi-img]][doi-url] <br> [![iif-cov-img]][iif-cov-url] <br> [![issues-open]][issues-url] | [![cjl-slack-badge]][cjl-slack] <br> [![caesar-docs]][cjl-docs-url] <br> [![dfg-docs]][dfg-docs-url] |


[iif-deps-img]: https://juliahub.com/docs/IncrementalInference/deps.svg
[iif-deps-jlh]: https://juliahub.com/ui/Packages/IncrementalInference/NrVw2??page=2
[doi-img]: https://zenodo.org/badge/DOI/10.5281/zenodo.5146221.svg
[doi-url]: https://doi.org/10.5281/zenodo.5146221

<!-- replicated in Caesar.jl README -->
[iif-ci-dev-img]: https://github.com/JuliaRobotics/IncrementalInference.jl/actions/workflows/ci.yml/badge.svg
[iif-ci-dev-url]: https://github.com/JuliaRobotics/IncrementalInference.jl/actions/workflows/ci.yml
[iif-ci-stb-img]: https://github.com/JuliaRobotics/IncrementalInference.jl/actions/workflows/ci.yml/badge.svg?branch=release%2Fv0.26
[iif-ci-stb-url]: https://github.com/JuliaRobotics/IncrementalInference.jl/actions/workflows/ci.yml
[iif-ver-img]: https://juliahub.com/docs/IncrementalInference/version.svg
[iif-rel-url]: https://github.com/JuliaRobotics/IncrementalInference.jl/releases
[iif-milestones]: https://github.com/JuliaRobotics/IncrementalInference.jl/milestones
[iif-cov-img]: https://codecov.io/github/JuliaRobotics/IncrementalInference.jl/coverage.svg?branch=master
[iif-cov-url]: https://codecov.io/github/JuliaRobotics/IncrementalInference.jl?branch=master

[iif-commits-url]: https://img.shields.io/github/commit-activity/y/JuliaRobotics/IncrementalInference.jl.svg?color=dark-green
[contributors-url]: https://github.com/JuliaRobotics/IncrementalInference.jl/graphs/contributors
[issues-time]: https://isitmaintained.com/badge/resolution/JuliaRobotics/IncrementalInference.jl.svg
[issues-open]: https://isitmaintained.com/badge/open/JuliaRobotics/IncrementalInference.jl.svg
[issues-url]: https://github.com/JuliaRobotics/IncrementalInference.jl/issues

[caesar-docs]: https://img.shields.io/badge/CaesarDocs-latest-blue.svg
[cjl-docs-url]: http://juliarobotics.github.io/Caesar.jl/latest/
[dfg-docs]: https://img.shields.io/badge/DFGDocs-latest-blue.svg
[dfg-docs-url]: https://juliarobotics.org/DistributedFactorGraphs.jl/latest/

[cjl-slack-badge]: https://img.shields.io/badge/Caesarjl-Slack-green.svg?style=popout
[cjl-slack]: https://join.slack.com/t/caesarjl/shared_invite/zt-ucs06bwg-y2tEbddwX1vR18MASnOLsw

Optimization routines for incremental non-parametric and parametric solutions based on factor graphs and the Bayes (Junction) tree implemented in the [Julia language](http://www.julialang.org/) (and [JuliaPro](http://www.juliacomputing.com)).


Introduction
============

This package implements a few different non-Gaussian factor graph inference algorithms, primarily 
- Multi-Modal iSAM (MM-iSAM) ([see references](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/)) which does hybrid non-parametric and parametric inference/state-estimation over large factor graphs.  
- Batch Parametric (akin to conventional "non-linear least squares"),
- Max-mixtures parametric,
- Other multiparametric and non-Gaussian algorithms are in the works and will be announced in due course.

Fundamentally, inference is performed via the Bayes (junction) tree where Chapman-Kolmogorov transit integral solutions are based on marginal-joint belief estimation (a sum-product / belief-propagation approximation algorithm).  Many benefits such as clique recycling are also available.  See the common Caesar.jl documenation for more details.  [![caesar-docs]][cjl-docs-url]

This package forms a cardinal piece of the [Caesar.jl](https://github.com/JuliaRobotics/Caesar.jl) robotics toolkit, including 3D visualization and database interaction, which can serve as a base station for a robotic platform. A standalone Robot Motion Estimate ([RoME.jl](https://github.com/JuliaRobotics/RoME.jl)) package extends the available variables, factors, and utilities for use in robotic navigation.  [![iif-deps-img]][iif-deps-jlh]  

Note, that IncrementalInference.jl **does not** have to be used with RoME.jl / Caesar.jl -- IncrementalInference.jl only implements the algebraic inference operations against mathematical abstractions such as Manifolds.jl. 

Furthermore, please contact info@navability.io for more formal support on this package, [NavAbility.io](https://www.navability.io). 

Installation
------------

Install the package from inside Julia
```julia
(v1.6) pkg> add IncrementalInference
```

Pre-install the following packages system wide packages[, and easily draw factor graph and Bayes tree]:
```bash
sudo apt-get install hdf5-tools
sudo apt-get install graphviz xdot # optional
```

Examples
========

See the common Caesar.jl documenation for more details [![caesar-docs]][cjl-docs-url].  Further examples can be found in the examples and test folders.

Cite and Contributors
=====================

We are grateful for many, many contributions within the Julia package ecosystem -- see the [`Project.toml`](https://github.com/JuliaRobotics/Caesar.jl/blob/master/Project.toml) files for a far reaching list of upstream packages and contributions.

Consider citing our work using the common reference at [Caesar.jl Citation with IncrementalInference.jl DOI](https://github.com/JuliaRobotics/Caesar.jl#contributors)

Get Involved, and Code of Conduct
---------------------------------

This project adheres to the [JuliaRobotics code of conduct](https://github.com/JuliaRobotics/administration/blob/master/code_of_conduct.md), and we invite contributions or comments from the community.  Use the slack channel, Julia Discourse, or Github issues to get in touch.

References
==========

See [references of interest here](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/)
