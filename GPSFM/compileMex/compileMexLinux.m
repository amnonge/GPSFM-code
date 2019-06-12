%This code compiles the mex BAceres
%This code assumes that ceres library: http://ceres-solver.org/
%is already installed on the computer.
%The library should be installed following the instructions here:
%http://ceres-solver.org/installation.html#linux
%(was tested on Ubuntu)

mex BAceres.cpp -I/usr/include/eigen3/  -lceres -lglog -llapack -lblas -largeArrayDims
