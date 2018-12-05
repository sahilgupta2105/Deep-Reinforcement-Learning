# 2D Physics Simulator for coupled Rigid-Fluid motion

This folder contatins the code for a 2D simulator which simulates the coupled motion of a fluid and rigid ball. The method is based on this [paper](https://www.cs.ubc.ca/labs/imager/tr/2007/Batty_VariationalFluids).

The code can be built using the following commands:
```
swig -c++ -python simulator_interface.i

python setup.py build_ext --inplace
```
The build requires Python 3.6+ and C++11. The reinforcement learning libraries are in Python. So, the C++ application is packaged as a Python module using an interface compiler, SWIG. The build also requires a reading/writing library called Partio, which can be installed from [here](https://github.com/wdas/partio). 

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library for linear algebra is used and is already provided with this code.

This code uses a Poisson Disk Sampler, which is originally available [here](https://github.com/thinks/poisson-disk-sampling).

The simulator supports two types of boundaries: rectangular and circular. 
