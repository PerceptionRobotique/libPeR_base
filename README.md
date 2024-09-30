[![Github Releases](https://img.shields.io/github/release/PerceptionRobotique/libPeR_base.svg)](https://github.com/PerceptionRobotique/libPeR_base/releases)

# LibPeR: Perception library for vision-based state estimation of robots and scenes

**This is the UniphorM branch of the LibPeR library.**

This branch main features contain:

- The Delaunay triangulation of the 3D points of the icosahedron
- The Voronoi diagram of the 3D points of the icosahedron
- Above mentioned spherical grid are available as `.xml` files in the `data/` directory (automatically moved at install time)
- Up to 7 subdivisions of the icosahedron (images up to *570x285*!)

## Get started

### Installation

This library relies on the following dependencies:

- [ViSP](https://visp.inria.fr/) (tested with 3.5)
- [OpenCV](https://opencv.org/) (tested with 4.2)
- LibXML2

To compile the library, run:

```
mkdir build && cd build
cmake ..
make -j12
sudo make install
```

### Documentation

The documentation of the library can be compiled by running **from the build/ directory**:

```
make html-doc
```

The `index` file of the documentation is generated in `build/html/index.html`

## Projects using LibPeR

### [**VisualGyroscope**](https://github.com/PerceptionRobotique/VisualGyroscope)

### [**VisualServoing**](https://github.com/PerceptionRobotique/VisualServoing)

### [**ROS Direct Visual Servoing**](https://github.com/isri-aist/ros_direct_visual_servoing)

A series of short programs that highlights how **LibPeR** can be used to perform **Direct Visual Alignment**, using various methods of spherical image mapping.

### Omnidirectional image utility tools

- **[dual2dual](https://github.com/PerceptionRobotique/dual2dual)** that warps a dual fisheye image to a given pose and outputs it as a dual fisheye image
- **[equi2equi](https://github.com/PerceptionRobotique/equi2equi)** that warps an equirectangular omnidirectional image to a given pose and outputs it as an equirectangular image
- **[equi2omni](https://github.com/PerceptionRobotique/equi2omni)** that warps an equirectangular omnidirectional image to a given pose and outputs it as a dual fisheye image
- **[dualfisheye2equi](https://github.com/PerceptionRobotique/dualfisheye2equi)** that warps a dual fisheye image to a given pose and outputs it as an equirectangular image

## Credits

```
This software was developed at:
MIS - UPJV
33 rue Saint-Leu
80039 AMIENS CEDEX
France

and at
CNRS - AIST JRL (Joint Robotics Laboratory)
1-1-1 Umezono, Tsukuba, Ibaraki
Japan

This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.

Description:
Insight about how to set the project and build the program
Authors:
Guillaume CARON, Antoine ANDRE

```
