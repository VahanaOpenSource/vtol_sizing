#
#
# HYBRID DESIGN AND ROTORCRAFT ANALYSIS (HYDRA)
#
#
Dependencies
---------------

 - Fortran 90 compiler (tested with gfortran)  - brew install gcc, sudo apt-get install gfortran
 - Python (python3 on mac), f90wrap, numpy, scipy, matplotlib, tkinter packages, pyyaml
 - CMake (minimum version 3.0 required)
 - CMake does check for the dependencies, so the missing tools will be known

Steps to build and compile
---------------
 - mkdir cmake/build
 - cd cmake/build
 - FC=gfortran cmake ../../
 - make -j

Test case
---------------
 - Default inputs 'default.yaml' - contains constants used for component models 
 - User    inputs   'input.yaml' - user definition for mission, aircraft type and sizing loops
- User, theory manuals coming soon.
 
Usage
------
 - cd cases/Beta/
 - python3 xrun.py (python2 backward support is present, but needs a few extra steps on Mac)
 - outputs are in outputs/logs/log*.yaml (open with text editor)
 - Postprocessing scripts are in cases/Beta/Postprocessing/
