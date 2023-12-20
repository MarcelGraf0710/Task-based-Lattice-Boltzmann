# Task-based-Lattice-Boltzmann
This repository contains the elaborations for my SimTech project.

## How to compile
This project uses CMake to create an executable file.
All parallel implementations will use HPX and thus require CMake to know the path to the HPX installation.
Specification of this path can be achieved on two ways:
- Add when accessing make, use