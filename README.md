Erosion thickness
=================

This software implements the method described in the following paper:

**Erosion Thickness on Medial Axes of 3D Shapes**.
ACM Transactions on Graphics, SIGGRAPH 2016. 
Yajie Yan (Washington University in St Louis), Kyle Sykes, Erin Chambers, David Letscher (St. Louis University), Tao Ju (Washington University in St Louis).

It computes a measure called "erosion thickness" on the medial axis of a 3d shape, and optionally a skeleton generated under the guidance of the measure.

Dev environment and dependencies
------------------------------------
This software is developed with the help of a set of 3rd party libraries, using Visual Studio 2015 on a Windows 10 machine. Specific version information about the dependencies and VS are as follows.

Dependencies (more can be found in `ET/third_party/version.txt`):
- qt: Qt5.7.0
- boost: boost-1.62
- cgal: CGAL-4.9 (compiled with boost-1.62)
- eigen: Eigen3.3.0-master-26667be4f70b
- gflags: gflags-master-8935ef4
- trimesh: trimesh2-2.12
- oglplus: oglplus-0.45
- glew: glew-2.0

Visual Studio 2015 version: 14.0.25123.0

New updated version supports Visual Studio 2017  

How to build (windows only)
---------------------------
1. All dependencies **except QT** is pre-shipped in the folder `third_party`. Therefore, the most effort amounts to installing QT of the right version (see above) on your machine. 
2. After installation, define the environment variable QTDIR that points to the root dir of your QT. This way VS can pick up the QT headers and libs. 
3. Finally, setup a symbolic link:
`C:\Libs\Qt\Qt5.7.0\5.7\msvc2015_64\bin` -> `bin` folder of your QT
as VS needs this path to call QT exes, e.g. moc, to properlly preprocess QT files.

After these steps are completed, the VS `.sln` file should build.

Platform other than Windows
---------------------------
The code is written with minimal win-specific functionality (as far as I can recall, only the memory profiling uses win api). Therefore, it's possible to compile the code on a Linux or Mac machine, although all dependencies need to be setup mannually.

Acknowledgement
---------------
- Xiaolong Zhang (xlzhang@cs.hku.hk): setting up a repo in bitbucket, and co-debugging.
- Dan Zeng (danzeng@wustl.edu): co-debugging, maintenance

Contact
-------
yajieyan@wustl.edu
