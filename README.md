# Blockster

(c) Timofey Mukha

Essentially, Blockster is a port of the CFD mesh generator blockMesh which
comes together with OpenFOAM.
The goal is to have most of the capabilities of blockMesh available in
Blockster, and then build on top of that.

At this point the main motivation for this project is to learn Julia "by doing",
and see how Blockster will compare to blockMesh in terms of speed.


## blockMesh features not yet in Blockster

* Fast algorithm for block merging.

* mergeFacePairs()

* Support for mesh regions, cell-zones etc.

* Support for non-hexahedral cells, non-quad faces.

* Curvelinear edges.

* Adding a default patch.

* Support for naming vertices

* Support for defining geometry and projecting mesh onto it.

## Features unique to Blockster

* Ability to define variables and use them for defining vertex position,
number of cells, expansion ratios etc.
