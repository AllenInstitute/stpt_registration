# Overview

This repository contains code for methods to register serial two-photon tomography (STPT) images
from the Allen Mouse Connectivity Atlas as described in 
[Oh et. al., 2014](https://www.ncbi.nlm.nih.gov/pubmed/24695228/) and
[Kuan et. al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/25536338/).
These methods have also been used in the construction of average template of the 
[Allen Mouse Common Coordinate Framework](https://community.brain-map.org/t/allen-mouse-ccf-accessing-and-using-related-data-and-tools/359).

# STPT Registration Methods

The three main executables are:

## idpProjectionAlignmentModule
This module performs a "global" registration between the STPT images and a target reference.
The module first performs a rigid registration followed by an affine (12 parameters) registration. In both cases using 
mutual information as the similarity metric and employing a multi-resolution optimization strategy. Input to this module is an xml file
containing meta-data for the image-series including the location, size and resolution of each image file. 
The output of this module is a "transform" xml.
For each image, it contains the 2D transform to reconstruction the specimen. For STPT, this is simply the identity transform since
the data is already inherently 3D. 
For the image series, it contains the 3D transform that will align the STPT specimen to the target reference.

Source: [lims2_modules/common/register/Registration/idpProjectionAlignmentModule.cxx](lims2_modules/common/register/Registration/idpProjectionAlignmentModule.cxx)


## idpProjectionLocalAlignmentModule
Following the "global" registration step, this module performs a "local" deformable registration
to improve the fine scale alignment of the STPT images to the target reference. 
In this method, deformation was parameterized as 3D B-splines (order of 3) where the knots are placed on a regularly spaced grid.
A four level coarse-to-fine optimization strategy was used. 
A coarse grid implicitly only allows smooth/stiff deformations while a finer grid allows more elastic transformations.
The method starts with a coarse grid allowing large scale anatomy to be matched.
The result is then used to initialize the finer levels allowing local neighborhoods to be more accurately matched. 

At each level, cross-correlation is used as the similarity metric. To speed up computation,
the metric was only evaluated on a small portion of randomly selected voxels.
The optimization problem was formulated and solved as a graph labeling problem. 
Specifically, instead of considering all the possible continuous moves of each knot, 
the search space is discretized into a number of grids represented by labels. 
The registration was performed in both the forward and backward direction in each iteration and then composed to form the final transform, 
resulting in a symmetric and invertible deformation field. 

Input to this module are two xml files, the first containing image-series metadata and the second contains the affine transforms 
computed during the "global" registration step. The output of this module is the final B-spline grid and associated deformation field.

Source: [lims2_modules/common/register/Registration/idpProjectionLocalAlignmentModule.cxx](lims2_modules/common/register/Registration/idpProjectionLocalAlignmentModule.cxx)


## idpProjectionResampleVolumeRGBModule
Following the "local" deformable registration step, this module enable the warping
of the STPT images to create an output volume that is aligned to the target reference. Input to this module
are the image-series meta-data and transform xml files as well as the deformation field computed during "local" deformable registration.

Source: [lims2_modules/common/register/Registration/idpProjectionResampleVolumeRGBModule.cxx](lims2_modules/common/register/Registration/idpProjectionResampleVolumeRGBModule.cxx)

# The "model_directory"

Reference data and algorithm parameters needed by the registration methods is organized in a directory structure.
In the production pipeline, the modules access this directory at a standard default location. Each module has an optional "model_directory" input argument to configure the modules to look for that information in user specified location. 

A copy of the content of the "model_directory" used for registering STPT images is located [here](http://download.alleninstitute.org/publications/allen_mouse_brain_common_coordinate_framework/).


# Level of Support

These methods were implemented as part of an informatics data processing pipeline and as such the 
design is influenced by the need to communicate with the internal database, workflow and build management
systems and to operate efficiently with internal file formats (eg JPEG2000) and compute cluster. A high level of technical
expertise may be required to compile the code and all its dependent packages, 
format image files into the require format and the generation and parsing of highly structured and verbose input/output xml files.

We have not deployed these methods beyond the scope of our own STPT images and computing cluster.

We are not currently supporting this code, but simply releasing it to the community AS IS 
but are not able to provide any guarantees of support. 
The community is welcome to submit issues, but you should not expect an active response.

# Dependencies

These methods leverage several libraries. While many are open source, some require specific licensing agreements.
Please see the documentation for each library for further details.

* **ITK** https://itk.org/, Version: InsightToolkit-3.20.0

* **VTK** https://vtk.org/, Version: VTK-5.10.0

* **Kakadu** https://kakadusoftware.com/, Version: v6_4_1-00465M

* **Boost** https://www.boost.org/, Version: boost_1_54_0

* **TinyXML** http://www.grinninglizard.com/tinyxml/, Version: tinyxml.2.6.2

* **CMake** https://cmake.org/, Version: cmake-2.8


# Build Instructions

* Download each dependency and build according to the instructions for each package
* Update the following files with the location of each dependency
  - [lims2_modules/lib/ITKJP2Lib/build/configure.sh](lims2_modules/lib/ITKJP2Lib/build/configure.sh)
  - [lims2_modules/common/register/build/configure.sh](lims2_modules/common/register/build/configure.sh)
* Change directory to [lims2_modules/common/register/](lims2_modules/common/register/)
* Run command [build/build_projection.sh](lims2_modules/common/register/build/build_projection.sh)


# Example Runtimes

* Example global alignment job
  - Limits:	mem=6gb,vmem=6gb,walltime=05:00:00,neednodes=centos
  - Resources: cput=00:05:49,vmem=1265748kb,walltime=00:06:04,mem=1731724kb
  - Node: x86_64, Linux, 3.10.0-514.26.2.el7.x86_64, 32 CPUs, 3000 MHz

* Example local alignment job including image resampling
  - Limits: mem=20gb,vmem=20gb,walltime=05:00:00,neednodes=centos
  - Resources: cput=00:45:11,vmem=17240312kb,walltime=00:45:42,mem=14007624kb
  - Node: x86_64, Linux, 3.10.0-514.26.2.el7.x86_64, 24 CPUs, 2500 MHz

