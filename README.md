# Direct-Mesh-Booleans

<img width="11162" height="2338" alt="presentation_title_image" src="https://github.com/user-attachments/assets/6a5428c9-3b88-4c07-9cca-abbe1b1dc40b" />

This is the implementation of "Direct Mesh Booleans: A Step Towards Non-Restrictive Boolean Operations."

## Overview
This repository contains a library (DirectMeshBooleans.lib) for computing Boolean operations. Currently, Union, Intersection, and Difference are supported. The inputs may contain open boundaries, non-manifold, or self-intersecting geometry. Our method works with **triangle meshes only**.

Optionally, an executable example that loads input meshes and runs a Boolean operation can be built. 
## Installation
### Prerequisites
This project is using CMake > 3.24 ((https://cmake.org/).

This library directly depends on:
- Eigen (https://gitlab.com/libeigen/eigen)
- OpenMesh (https://www.graphics.rwth-aachen.de/software/openmesh/)
- FaRMA (https://github.com/gcherchi/FastAndRobustMeshArrangements), which has additional internal dependencies.

By default, all the necessary dependencies are automatically downloaded and built, using the CMake system. 

### Setup
Tested on Windows 11. 

Simply clone the repository and run CMake > 3.24. This will **automatically** download and build all the dependencies.

If needed, you can use ```DMB_PROVIDE_OWN_DEPENDENCIES``` option to provide your own version of the dependencies.

## Usage
By default, only the ```DirectMeshBooleans``` library is built.
Option ```DMB_BUILD_EXAMPLES``` creates a target for a simple example application ```DirectMeshBooleansExample```, which loads inputs and runs a Boolean operation.
Example application usage:
```
DirectMeshBooleansExample.exe [U|D|I] input1 input2 outputPath

```
* ```[U|D|I]``` - Specifies requested Boolean operation **U**nion, **D**ifference or **I**ntersection.
* ```input1``` - Path to the first input mesh file.  
* ```input2``` - Path to the second input mesh file.
* ```output``` - Path to the output mesh file.

Supported file formats for input and output files are: OBJ, OFF, and STL.
## Citate us
TBD
