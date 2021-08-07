# IncrementalInference.jl

> Click on badges to follow links:

Stable v0.23 | Stable v0.24 | Dev | Coverage | Documentation
--------------|-------------|-------------|-----|---------
[![build-0-23]][travis-url] | [![build-0-24]][travis-url] | [![build-master]][travis-url] <br> [![commits-url]][contributors-url] <br> [![issues-time]][issues-url] | [![doi-img]][doi-url] <br> [![codecov-io]][codecov-url] <br> [![issues-open]][issues-url] | [![cjl-slack-badge]][cjl-slack] <br> [![caesar-docs]][cjl-docs-url] <br> [![dfg-docs]][dfg-docs-url]



[doi-img]: https://zenodo.org/badge/55802838.svg
[doi-url]: https://zenodo.org/badge/latestdoi/55802838

[travis-url]: https://travis-ci.org/JuliaRobotics/IncrementalInference.jl
[build-master]: https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=master
[build-0-23]: https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=release/v0.23
[build-0-24]: https://travis-ci.org/JuliaRobotics/IncrementalInference.jl.svg?branch=release/v0.24

[codecov-io]: https://codecov.io/github/JuliaRobotics/IncrementalInference.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JuliaRobotics/IncrementalInference.jl?branch=master
[commits-url]: https://img.shields.io/github/commit-activity/y/JuliaRobotics/IncrementalInference.jl.svg?color=dark-green
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

This package implements Multi-Modal iSAM (MM-iSAM) algorithm ([see references](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/)) which does hybrid non-parametric and parametric inference/state-estimation over large factor graphs.  Inference is performed via the Bayes (junction) tree where non-parametric and parametric solutions are based on marginal-joint belief propagation (an approximate sum-product algorithm).  Many benefits such as clique recycling is available.  See the common Caesar.jl documenation for more details.  Note, that IncrementalInference.jl **does not** have to be used with Caesar.jl, and aims to only implement the inference operations against mathematical abstractions such as Manifolds.jl. 

This package forms a cardinal piece of the [Caesar.jl](https://github.com/JuliaRobotics/Caesar.jl) robotics toolkit, including 3D visualization and database interaction, which can serve as a base station for a robotic platform. A standalone [Robot Motion Estimate](https://github.com/JuliaRobotics/RoME.jl) package extends the available variables, factors, and utilities for use in robotic navigation.

Please contact info@navability.io for further support on this package, [NavAbility](https://www.navability.io). 

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

This library is built as a back-end solver which is closer to the mathetical operations can be easily modified and extended for a variety of uses. Specific emphasis is placed on allowing outside user defined variables and factor definitions to be used. The current major use case is through [RoME.jl](http://github.com/JuliaRobotics/RoME.jl) and [Caesar.jl](http:///www.github.com/JuliaRobotics/Caesar.jl) which introduces various sensor measurement and motion manifold functions for use in Robot Motion Estimate (a.k.a SLAM).  See these and related packages for documentation and examples.

<p align="center">
<a href="https://vimeo.com/190052649" target="_blank"><img src="https://raw.githubusercontent.com/JuliaRobotics/IncrementalInference.jl/master/doc/images/mmfgbt.gif" alt="IMAGE ALT TEXT HERE" width="480" height="320" /></a>
</p>

Contributors
============

We are grateful for many, many contributions within the Julia package ecosystem -- see the [`Project.toml`](https://github.com/JuliaRobotics/Caesar.jl/blob/master/Project.toml) files for a far reaching list of upstream packages and contributions.

Consider citing our work using the common reference at [Caesar.jl Citation with IncrementalInference.jl DOI](https://github.com/JuliaRobotics/Caesar.jl#contributors)

Get Involved, and Code of Conduct
=================================

This project adheres to the [JuliaRobotics code of conduct](https://github.com/JuliaRobotics/administration/blob/master/code_of_conduct.md), and we invite contributions or comments from the community.  Use the slack channel, Julia Discourse, or Github issues to get in touch.

References
==========

See [references of interest here](http://www.juliarobotics.org/Caesar.jl/latest/refs/literature/)
