# Task-based-Lattice-Boltzmann
This repository contains the elaborations for my SimTech project.

## How to compile
This project uses CMake to create an executable file.
All parallel implementations use HPX (https://hpx-docs.stellar-group.org/latest/html/index.html) and thus require CMake to know the path to the HPX installation.
Every `CMakeLists.txt` file must include the line `find_package(HPX REQUIRED)`. 
This tells CMAKE that the project uses HPX and that the respective files need to be considered when building.
In order for CMAKE to find HPX, one must specify the `CMAKE_PREFIX_PATH` to the path leading to the file `HPXConfig.cmake`.
Assuming the `CMakeLists.txt` is located within the parent folder of the build folder, this is done by accessing CMAKE on the build directory with 
`cmake -DCMAKE_PREFIX_PATH=/path/to/HPXConfig.cmake ..`.
By default, this path looks something like this:
`~/Documents/spack/opt/spack/YOUR-LINUX-VERSION/YOUR-COMPILER-VERSION/HPX-VERSION-FOLLOWED-BY-GIBBERISH/lib/cmake/HPX`.

## General recommendations
If you want to use IntelliSense, I recommend making an addition to the `c_cpp_properties.json` file within the `.vscode` folder.
`"includePath"` usually contains `"${workspaceFolder}/**"` such that IntelliSense recursively searches through all files within the workspace folder.
Like this, HPX headers will not be found by IntelliSense.
To remedy this, set `"includePath"` to `["${workspaceFolder}/**,~/Documents/spack/opt/spack/YOUR-LINUX-VERSION/YOUR-COMPILER-VERSION/**]"`.
