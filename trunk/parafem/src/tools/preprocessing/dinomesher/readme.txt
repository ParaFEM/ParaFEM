This directory contains source code for the meshgeneration program developed by Peter Falkingham for use in dinosaur footprint simulations.

It has been only slighlty modified from the form used for personal use, though will hopefully be updated in the the future to include more functionality.

It's primary purpose is to generate a 3D cubic mesh composed of 20-node hexehedral elements.  Optionally, an indenter can be meshed on the surface in the shape of a picture defined by a PPM image file.

The main use of the program is to generate a cubic mesh of any dimensions and element resolution.  A scaling factor can be incorporated to produce a mesh where the central elements are smaller (higher resolution) than elements twoards the edges, alowing large meshes without reaching silly numbers of elements.

I have added limited functionality for generating meshes of 8-node hexahedral elements.  Using 8-node elements will result in no load file, and is incompatible with an indenter.  This is included only for users who wish to generate a cubic mesh of 8-node hex elements with the scaling algorithm, but wish to add their own loading conditions later.  I think at the moment you need a ppm file at the same resolution as the dense mesh - sorry.

I also have utilities for converting the meshes into Abaqus format.  These may be uploaded at a later date, or can be obtained by emailing peter.falkingham@manchester.ac.uk




Peter Falkingham 17/12/10
