# Barycentric-OpenCL

Originally developed as part of my [thesis
work](https://github.com/DrJPepper/4d-multilayer-modeling), this library is an
OpenCL optimized implementation of the barycentric calculations needed to work
with parameterized triangle meshes, such as what can be generated using
[CGAL](https://doc.cgal.org/latest/Surface_mesh_parameterization/index.html).
The primary functions provided by this library are `barycentricLoopingR` which
takes one or more points in model space and returns their equivalent in
parameter space, and `barycentricLoopingF` which does the reverse.

## Features

* Automatic instantiation of parameterization from input triangle mesh
    * Support for open surface models, closed models via ray tracing, and Bezier
      patches
* Bi-directional conversion between model and parameter space
* OpenCL-based parallelization for an arbitrary number of points

## Installation

The code is written in C++ and the project is configured to use Cmake for
compilation.

### Dependencies

* Eigen
* CGAL
* FMT

### Compilation

A demo is included and compiled automatically by Cmake. The demo has an
additional dependency of toml++. To compile this demo standalone, run

    mkdir build
    cd build
    cmake ..

The demo expects a TOML config file for a triangle mesh model. Examples of these
can be found in
[https://github.com/DrJPepper/4d-modeling-example-inputs](https://github.com/DrJPepper/4d-modeling-example-inputs).
From the same directory where Barycentric-OpenCL was cloned, run

    git clone https://github.com/DrJPepper/4d-modeling-example-inputs
    ln -sr 4d-modeling-example-inputs barycentric-opencl/inputs
    cd barycentric-opencl/build
    ./demo/barycentric_ocl_demo ../input/configs/4dp/egg_chair.toml

which should output

    0.325464 0.412296
    14.44 3.84312  11.217

### Library Only

The following `CMakeLists.txt` code can be used to compile the library for use in a larger project

    FetchContent_Declare(
        BarycentricOpenCL
        GIT_REPOSITORY https://github.com/DrJPepper/barycentric-opencl
        GIT_TAG        main
    )
    set(DISABLE_DEMO ON CACHE INTERNAL "Disable BarycentricOCL demo")
    FetchContent_MakeAvailable(BarycentricOpenCL)
    target_link_libraries(${PROJECT_NAME} PUBLIC BarycentricOpenCL)

This code automatically downloads the library to your build directory, and
disables compilation of the demo.

## Usage

As it is currently configured, an input model must be provided to
the `Mesh` class, which then automatically initializes the requisite
parameterization. The `Mesh` constructor takes two arguments, the name of the input model file 

## TODO

* ☐ Improve memory management and resource cleanup
* ☐ Check available VRAM and determine max settings
* ☐ Allow user to request either max memory usage or max concurrent points
* ☐ Run optimization check to see at what point it becomes worth running on OpenCL vs CPU
* ☐ Auto switch between CPU, OCL and multiple iterations of OCL
* ☐ Decouple barycentric code further by allowing user to pass in their own vertices and faces instead of only working with input models
* ☐ Stop just returning 0,0/0,0,0 when barycentric calculations fail/error
* ☐ Fall back to CPU when a point fails on OCL and print a more detailed error message
<!--☑-->
